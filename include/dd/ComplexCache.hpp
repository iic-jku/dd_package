/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_COMPLEXCACHE_HPP
#define DD_PACKAGE_COMPLEXCACHE_HPP

#include "Complex.hpp"

#include <cassert>
#include <cstddef>
#include <vector>

namespace dd {

    template<std::size_t INITIAL_ALLOCATION_SIZE = 2048, std::size_t GROWTH_FACTOR = 2>
    class ComplexCache {
        using MagnitudeEntry = MagnitudeTable<>::Entry;
        using PhaseEntry     = PhaseTable<>::Entry;

    public:
        ComplexCache():
            chunkID(0), allocationSize(INITIAL_ALLOCATION_SIZE) {
            // allocate first chunk of cache entries
            magChunks.emplace_back(allocationSize);
            phaseChunks.emplace_back(allocationSize);
            allocations += allocationSize;
            allocationSize *= GROWTH_FACTOR;
            magChunkIt      = magChunks[0].begin();
            magChunkEndIt   = magChunks[0].end();
            phaseChunkIt    = phaseChunks[0].begin();
            phaseChunkEndIt = phaseChunks[0].end();
        }

        ~ComplexCache() = default;

        // access functions
        [[nodiscard]] std::size_t getCount() const { return count; }
        [[nodiscard]] std::size_t getPeakCount() const { return peakCount; }
        [[nodiscard]] std::size_t getAllocations() const { return allocations; }
        [[nodiscard]] std::size_t getGrowthFactor() const { return GROWTH_FACTOR; }

        [[nodiscard]] Complex getCachedComplex() {
            // an entry is available on the stack
            if (magAvailable != nullptr || phaseAvailable != nullptr) {
                assert(magAvailable != nullptr && phaseAvailable != nullptr);
                auto entry     = Complex{magAvailable, phaseAvailable};
                magAvailable   = magAvailable->next;
                phaseAvailable = phaseAvailable->next;
                count += 2;
                return entry;
            }

            // new chunk has to be allocated
            if (magChunkIt == magChunkEndIt || phaseChunkIt == phaseChunkEndIt) {
                assert(magChunkIt == magChunkEndIt && phaseChunkIt == phaseChunkEndIt);
                magChunks.emplace_back(allocationSize);
                phaseChunks.emplace_back(allocationSize);
                allocations += allocationSize;
                allocationSize *= GROWTH_FACTOR;
                chunkID++;
                magChunkIt      = magChunks[chunkID].begin();
                magChunkEndIt   = magChunks[chunkID].end();
                phaseChunkIt    = phaseChunks[chunkID].begin();
                phaseChunkEndIt = phaseChunks[chunkID].end();
            }

            Complex c{};
            c.mag = &(*magChunkIt);
            ++magChunkIt;
            c.phase = &(*phaseChunkIt);
            ++phaseChunkIt;
            count += 2;
            return c;
        }

        [[nodiscard]] Complex getTemporaryComplex() {
            // an entry is available on the stack
            if (magAvailable != nullptr || phaseAvailable != nullptr) {
                assert(magAvailable != nullptr && phaseAvailable != nullptr);
                return {magAvailable, phaseAvailable};
            }

            // new chunk has to be allocated
            if (magChunkIt == magChunkEndIt || phaseChunkIt == phaseChunkEndIt) {
                assert(magChunkIt == magChunkEndIt && phaseChunkIt == phaseChunkEndIt);
                magChunks.emplace_back(allocationSize);
                phaseChunks.emplace_back(allocationSize);
                allocations += allocationSize;
                allocationSize *= GROWTH_FACTOR;
                chunkID++;
                magChunkIt      = magChunks[chunkID].begin();
                magChunkEndIt   = magChunks[chunkID].end();
                phaseChunkIt    = phaseChunks[chunkID].begin();
                phaseChunkEndIt = phaseChunks[chunkID].end();
            }
            return {&(*magChunkIt), &(*(phaseChunkIt))};
        }

        void returnToCache(Complex& c) {
            assert(count >= 2);
            assert(c != Complex::zero);
            assert(c != Complex::one);
            assert(c.mag->refCount == 0);
            assert(c.phase->refCount == 0);

            c.phase->next  = phaseAvailable;
            phaseAvailable = c.phase;

            c.mag->next  = magAvailable;
            magAvailable = c.mag;
            count -= 2;
        }

        void clear() {
            // clear available stack
            magAvailable   = nullptr;
            phaseAvailable = nullptr;

            // release memory of all but the first chunk TODO: it could be desirable to keep the memory
            while (chunkID > 0) {
                magChunks.pop_back();
                phaseChunks.pop_back();
                chunkID--;
            }
            // restore initial chunk setting
            phaseChunkIt    = phaseChunks[0].begin();
            phaseChunkEndIt = phaseChunks[0].end();
            magChunkIt      = magChunks[0].begin();
            magChunkEndIt   = magChunks[0].end();
            allocationSize  = INITIAL_ALLOCATION_SIZE * GROWTH_FACTOR;
            allocations     = INITIAL_ALLOCATION_SIZE;

            count     = 0;
            peakCount = 0;
        };

    private:
        MagnitudeEntry*                                magAvailable{};
        PhaseEntry*                                    phaseAvailable{};
        std::vector<std::vector<MagnitudeEntry>>       magChunks{};
        std::vector<std::vector<PhaseEntry>>           phaseChunks{};
        std::size_t                                    chunkID;
        typename std::vector<MagnitudeEntry>::iterator magChunkIt;
        typename std::vector<MagnitudeEntry>::iterator magChunkEndIt;
        typename std::vector<PhaseEntry>::iterator     phaseChunkIt;
        typename std::vector<PhaseEntry>::iterator     phaseChunkEndIt;
        std::size_t                                    allocationSize;

        std::size_t allocations = 0;
        std::size_t count       = 0;
        std::size_t peakCount   = 0;
    };
} // namespace dd

#endif //DD_PACKAGE_COMPLEXCACHE_HPP
