/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_PHASETABLE_HPP
#define DD_PACKAGE_PHASETABLE_HPP

#include "dd/Definitions.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace dd {
    template<std::size_t NBUCKET = 65537, std::size_t INITIAL_ALLOCATION_SIZE = 2048, std::size_t GROWTH_FACTOR = 2, std::size_t INITIAL_GC_LIMIT = 65536>
    class PhaseTable {
    public:
        struct Entry {
            fp       value{};
            Entry*   next{};
            RefCount refCount{};

            ///
            /// The sign of number is encoded in the least significant bit of its entry pointer
            /// If not handled properly, this causes misaligned access
            /// These routines allow to obtain safe pointers
            ///
            [[nodiscard]] static constexpr inline Entry* getAlignedPointer(const Entry* e) {
                return reinterpret_cast<Entry*>(reinterpret_cast<std::uintptr_t>(e) & (~3UL));
            }

            [[nodiscard]] static constexpr inline std::uint_fast8_t getQuadrant(const Entry* e) {
                return reinterpret_cast<std::uintptr_t>(e) & 3ULL;
            }

            [[nodiscard]] static constexpr inline Entry* setQuadrant(const Entry* e, std::uint_fast8_t newQuadrant) {
                return reinterpret_cast<Entry*>((reinterpret_cast<std::uintptr_t>(e) & ~3ULL) | newQuadrant);
            }

            [[nodiscard]] static inline fp val(const Entry* e) {
                switch (getQuadrant(e)) {
                    case 0: return e->value;
                    case 1: return getAlignedPointer(e)->value + 0.5;
                    case 2: return getAlignedPointer(e)->value + 1.0;
                    case 3: return getAlignedPointer(e)->value + 1.5;
                }
                throw std::runtime_error("This should be unreachable");
            }

            [[nodiscard]] static inline bool approximatelyEquals(const Entry* left, const Entry* right) {
                return std::fmod(std::abs(val(left) - val(right)), 2) < TOLERANCE;
            }

            [[nodiscard]] static inline bool approximatelyZero(const Entry* e) {
                return e == zeroPtr || std::fmod(std::abs(val(e)), 2) < TOLERANCE;
            }

            [[nodiscard]] static inline bool approximatelyOne(const Entry* e) {
                return e == onePtr || std::fmod(std::abs(val(e) - 1), 2) < TOLERANCE;
            }

            static void writeBinary(const Entry* e, std::ostream& os) {
                auto temp = val(e);
                os.write(reinterpret_cast<const char*>(&temp), sizeof(decltype(temp)));
            }
        };

        static inline Entry zero{0., nullptr, 1};

        static inline Entry* zeroPtr = &zero;
        static inline Entry* onePtr  = reinterpret_cast<Entry*>(reinterpret_cast<std::uintptr_t>(&zero) | 2u);

        PhaseTable():
            chunkID(0), allocationSize(INITIAL_ALLOCATION_SIZE), gcLimit(INITIAL_GC_LIMIT) {
            // allocate first chunk of numbers
            chunks.emplace_back(allocationSize);
            allocations += allocationSize;
            allocationSize *= GROWTH_FACTOR;
            chunkIt    = chunks[0].begin();
            chunkEndIt = chunks[0].end();

            // add 1/4 to the complex table and increase its ref count (so that it is not collected)
            lookup(0.25)->refCount++;
        }

        ~PhaseTable() = default;

        static fp tolerance() {
            return TOLERANCE;
        }

        static void setTolerance(fp tol) {
            TOLERANCE = tol;
        }

        static constexpr std::size_t MASK = NBUCKET - 1;

        // linear (clipped) hash function
        static constexpr std::size_t hash(const fp val) {
            assert(val >= 0);
            assert(val < 0.5);
            auto key = static_cast<std::size_t>(val * (2 * MASK) + static_cast<dd::fp>(0.5));
            return std::min<std::size_t>(key, MASK);
        }

        // access functions
        [[nodiscard]] std::size_t getCount() const { return count; }

        [[nodiscard]] std::size_t getPeakCount() const { return peakCount; }

        [[nodiscard]] std::size_t getAllocations() const { return allocations; }

        [[nodiscard]] std::size_t getGrowthFactor() const { return GROWTH_FACTOR; }

        [[nodiscard]] const auto& getTable() const { return table; }

        [[nodiscard]] bool availableEmpty() const { return available == nullptr; }

        Entry* lookupQuadrantized(const fp& val) {
            assert(val >= 0); // should be handled as special case before quadrantizing

            if (val <= TOLERANCE) {
                ++hits;
                return zeroPtr;
            }

            const std::size_t lowerKey = hash(val - TOLERANCE);
            const std::size_t upperKey = hash(val + TOLERANCE);

            if (upperKey == lowerKey) {
                ++findOrInserts;
                return findOrInsert(lowerKey, val);
            }

            // code below is to handle cases where the looked up value
            // could be in the lower or upper buckets and we have to go through them

            const std::size_t key = hash(val);

            Entry* p = find(table[key], val);
            if (p != nullptr) {
                return p;
            }

            // search in (potentially) lower bucket
            if (lowerKey != key) {
                ++lowerNeighbors;
                // buckets are sorted so we only have to look into the last entry of the lower bucket
                Entry* p_lower = tailTable[lowerKey];
                if (p_lower != nullptr && val - p_lower->value < TOLERANCE) {
                    return p_lower;
                }
            }

            // search in (potentially) higher bucket
            if (upperKey != key) {
                ++upperNeighbors;
                // buckets are sorted, we only have to look at the first element

                Entry* p_upper = table[upperKey];
                if (p_upper != nullptr && p_upper->value - val < TOLERANCE) {
                    return p_upper;
                }
            }

            // value was not found in the table -> get a new entry and add it to the central bucket
            Entry* entry = insert(key, val);
            return entry;
        }

        Entry* lookup(fp val) {
            assert(!std::isnan(val));
            ++lookups;
            //std::cout << val;
            if (std::abs(val) <= TOLERANCE) {
                ++hits;
                //std::cout << " is close to zero\n";
                return zeroPtr;
            }

            while (val < 0) {
                val += 2;
            }

            uint_fast8_t quadrant = 0;
            while (val >= 0.5) {
                quadrant++;
                val -= 0.5;
            }

            //std::cout << " --> " << val << " @ " << +quadrant << std::endl;
            auto entry = lookupQuadrantized(val);
            return Entry::setQuadrant(entry, quadrant % 4);
        }

        [[nodiscard]] Entry* getEntry() {
            // an entry is available on the stack
            if (!availableEmpty()) {
                Entry* entry = available;
                available    = entry->next;
                // returned entries could have a ref count != 0
                entry->refCount = 0;
                return entry;
            }

            // new chunk has to be allocated
            if (chunkIt == chunkEndIt) {
                chunks.emplace_back(allocationSize);
                allocations += allocationSize;
                allocationSize *= GROWTH_FACTOR;
                chunkID++;
                chunkIt    = chunks[chunkID].begin();
                chunkEndIt = chunks[chunkID].end();
            }

            auto entry = &(*chunkIt);
            ++chunkIt;
            return entry;
        }

        void returnEntry(Entry* entry) {
            entry->next = available;
            available   = entry;
        }

        // increment reference count for corresponding entry
        static void incRef(Entry* entry) {
            // get valid pointer
            auto entryPtr = Entry::getAlignedPointer(entry);

            if (entryPtr == nullptr)
                return;

            // important (static) numbers are never altered
            if (entryPtr != &zero) {
                if (entryPtr->refCount == std::numeric_limits<RefCount>::max()) {
                    std::clog << "[WARN] PhaseTable: MAXREFCNT " << entryPtr->refCount << " reached for " << entryPtr->value << ". Number will never be collected." << std::endl;
                    return;
                }

                // increase reference count
                entryPtr->refCount++;
            }
        }

        // decrement reference count for corresponding entry
        static void decRef(Entry* entry) {
            // get valid pointer
            auto entryPtr = Entry::getAlignedPointer(entry);

            if (entryPtr == nullptr)
                return;

            // important (static) numbers are never altered
            if (entryPtr != &zero) {
                if (entryPtr->refCount == std::numeric_limits<RefCount>::max()) {
                    return;
                }
                if (entryPtr->refCount == 0) {
                    throw std::runtime_error("In ComputeTable: RefCount of entry " + std::to_string(entryPtr->value) + " is zero before decrement");
                }

                // decrease reference count
                entryPtr->refCount--;
            }
        }

        [[nodiscard]] bool possiblyNeedsCollection() const { return count >= gcLimit; }

        std::size_t garbageCollect(bool force = false) {
            gcCalls++;
            // nothing to be done if garbage collection is not forced, and the limit has not been reached,
            // or the current count is minimal (the complex table always contains at least 0.5)
            if ((!force && count < gcLimit) || count <= 1)
                return 0;

            gcRuns++;
            std::size_t collected = 0;
            std::size_t remaining = 0;
            for (std::size_t key = 0; key < table.size(); ++key) {
                Entry* p     = table[key];
                Entry* lastp = nullptr;
                while (p != nullptr) {
                    if (p->refCount == 0) {
                        Entry* next = p->next;
                        if (lastp == nullptr) {
                            table[key] = next;
                        } else {
                            lastp->next = next;
                        }
                        returnEntry(p);
                        p = next;
                        collected++;
                    } else {
                        lastp = p;
                        p     = p->next;
                        remaining++;
                    }
                    tailTable[key] = lastp;
                }
            }
            // The garbage collection limit changes dynamically depending on the number of remaining (active) nodes.
            // If it were not changed, garbage collection would run through the complete table on each successive call
            // once the number of remaining entries reaches the garbage collection limit. It is increased whenever the
            // number of remaining entries is rather close to the garbage collection threshold and decreased if the
            // number of remaining entries is much lower than the current limit.
            if (remaining > gcLimit / 10 * 9) {
                gcLimit = remaining + INITIAL_GC_LIMIT;
            } else if (remaining < gcLimit / 128) {
                gcLimit /= 2;
            }
            count = remaining;
            return collected;
        }

        void clear() {
            // clear table buckets
            for (auto& bucket: table) {
                bucket = nullptr;
            }
            for (auto& entry: tailTable) {
                entry = nullptr;
            }

            // clear available stack
            available = nullptr;

            // release memory of all but the first chunk TODO: it could be desirable to keep the memory
            while (chunkID > 0) {
                chunks.pop_back();
                chunkID--;
            }
            // restore initial chunk setting
            chunkIt        = chunks[0].begin();
            chunkEndIt     = chunks[0].end();
            allocationSize = INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR;
            allocations    = INITIAL_ALLOCATION_SIZE;

            for (auto& entry: chunks[0]) {
                entry.refCount = 0;
            }

            count     = 0;
            peakCount = 0;

            collisions       = 0;
            insertCollisions = 0;
            hits             = 0;
            findOrInserts    = 0;
            lookups          = 0;
            inserts          = 0;
            lowerNeighbors   = 0;
            upperNeighbors   = 0;

            gcCalls = 0;
            gcRuns  = 0;
            gcLimit = INITIAL_GC_LIMIT;
        };

        void print() {
            for (std::size_t key = 0; key < table.size(); ++key) {
                auto p = table[key];
                if (p != nullptr)
                    std::cout << key << ": "
                              << "\n";

                while (p != nullptr) {
                    std::cout << "\t\t" << p->value << " " << reinterpret_cast<std::uintptr_t>(p) << " " << p->refCount
                              << "\n";
                    p = p->next;
                }

                if (table[key] != nullptr)
                    std::cout << "\n";
            }
        }

        [[nodiscard]] fp hitRatio() const { return static_cast<fp>(hits) / lookups; }

        [[nodiscard]] fp colRatio() const { return static_cast<fp>(collisions) / lookups; }

        std::map<std::string, std::size_t> getStatistics() {
            return {
                    {"hits", hits},
                    {"collisions", collisions},
                    {"lookups", lookups},
                    {"inserts", inserts},
                    {"insertCollisions", insertCollisions},
                    {"findOrInserts", findOrInserts},
                    {"upperNeighbors", upperNeighbors},
                    {"lowerNeighbors", lowerNeighbors},
                    {"gcCalls", gcCalls},
                    {"gcRuns", gcRuns},
            };
        }

        std::ostream& printStatistics(std::ostream& os = std::cout) {
            // clang-format off
            os << "hits: " << hits
               << ", collisions: " << collisions
               << ", looks: " << lookups
               << ", inserts: " << inserts
               << ", insertCollisions: " << insertCollisions
               << ", findOrInserts: " << findOrInserts
               << ", upperNeighbors: " << upperNeighbors
               << ", lowerNeighbors: " << lowerNeighbors
               << ", hitRatio: " << hitRatio()
               << ", colRatio: " << colRatio()
               << ", gc calls: " << gcCalls
               << ", gc runs: " << gcRuns
               << "\n";
            // clang-format on
            return os;
        }

        std::ostream& printBucketDistribution(std::ostream& os = std::cout) {
            for (auto bucket: table) {
                if (bucket == nullptr) {
                    os << "0\n";
                    continue;
                }
                std::size_t bucketCount = 0;
                while (bucket != nullptr) {
                    ++bucketCount;
                    bucket = bucket->next;
                }
                os << bucketCount << "\n";
            }
            os << std::endl;
            return os;
        }

    private:
        using Bucket = Entry*;
        using Table  = std::array<Bucket, NBUCKET>;

        Table table{};

        std::array<Entry*, NBUCKET> tailTable{};

        // table lookup statistics
        std::size_t collisions       = 0;
        std::size_t insertCollisions = 0;
        std::size_t hits             = 0;
        std::size_t findOrInserts    = 0;
        std::size_t lookups          = 0;
        std::size_t inserts          = 0;
        std::size_t lowerNeighbors   = 0;
        std::size_t upperNeighbors   = 0;

        // numerical tolerance to be used for floating point values
        static inline fp TOLERANCE = 1e-13;

        Entry*                                available{};
        std::vector<std::vector<Entry>>       chunks{};
        std::size_t                           chunkID;
        typename std::vector<Entry>::iterator chunkIt;
        typename std::vector<Entry>::iterator chunkEndIt;
        std::size_t                           allocationSize;

        std::size_t allocations = 0;
        std::size_t count       = 0;
        std::size_t peakCount   = 0;

        // garbage collection
        std::size_t gcCalls = 0;
        std::size_t gcRuns  = 0;
        std::size_t gcLimit = 100000;

        inline Entry* findOrInsert(const std::size_t key, const fp val) {
            [[maybe_unused]] const fp val_tol = val + TOLERANCE;

            Entry* curr = table[key];
            Entry* prev = nullptr;

            while (curr != nullptr && val_tol > curr->value) {
                if (std::abs(curr->value - val) < TOLERANCE) {
                    ++hits;
                    return curr;
                }
                ++collisions;
                prev = curr;
                curr = curr->next;
            }

            ++inserts;
            Entry* entry = getEntry();
            entry->value = val;

            if (prev == nullptr) {
                // table bucket is empty
                table[key] = entry;
            } else {
                prev->next = entry;
            }
            entry->next = curr;
            if (curr == nullptr) {
                tailTable[key] = entry;
            }
            count++;
            peakCount = std::max(peakCount, count);
            return entry;
        }

        /**
         * Inserts a value into the bucket indexed by key. This function assumes no element within TOLERANCE is
         * present in the bucket.
         * @param key index to the bucket
         * @param val value to be inserted
         * @return pointer to the inserted entry
         */
        inline Entry* insert(const std::size_t key, const fp val) {
            ++inserts;
            Entry* entry = getEntry();
            entry->value = val;

            Entry* curr = table[key];
            Entry* prev = nullptr;

            while (curr != nullptr && val < curr->value) {
                ++insertCollisions;
                prev = curr;
                curr = curr->next;
            }

            if (prev == nullptr) {
                // table bucket is empty
                table[key] = entry;
            } else {
                prev->next = entry;
            }
            entry->next = curr;
            if (curr == nullptr) {
                tailTable[key] = entry;
            }
            count++;
            peakCount = std::max(peakCount, count);
            return entry;
        }

        inline Entry* find(const Bucket& bucket, const fp& val) {
            Entry*   p       = bucket;
            const fp val_tol = val - TOLERANCE;
            while (p != nullptr && val_tol <= p->value) {
                if (p->value - val < TOLERANCE) {
                    ++hits;
                    return p;
                }
                ++collisions;
                p = p->next;
            }
            return nullptr;
        }
    };
} // namespace dd
#endif //DD_PACKAGE_PHASETABLE_HPP
