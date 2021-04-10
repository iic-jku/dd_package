/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DDpackage_UNIQUETABLE_HPP
#define DDpackage_UNIQUETABLE_HPP

#include "Definitions.hpp"

#include <array>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

namespace dd {

    template<class Edge, std::size_t NBUCKET = 32768, std::size_t ALLOCATION_SIZE = 2000>
    class UniqueTable {
        using Node = std::remove_pointer_t<decltype(Edge::p)>;

    public:
        explicit UniqueTable(std::size_t nvars, std::size_t gcLimit = 250000, std::size_t gcIncrement = 0):
            nvars(nvars), gcInitialLimit(gcLimit), gcLimit(gcLimit), gcIncrement(gcIncrement) {}

        ~UniqueTable() {
            for (auto chunk: chunks) delete[] chunk;
        }

        static constexpr size_t MASK = NBUCKET - 1;

        void resize(std::size_t nq) {
            nvars = nq;
            tables.resize(nq);
            // TODO: if the new size is smaller than the old one we might have to release the unique table entries for the superfluous variables
            active.resize(nq);
            activeNodeCount = std::accumulate(active.begin(), active.end(), 0UL);
        }

        static std::size_t hash(const Node* p) {
            std::uintptr_t key = 0;
            for (std::size_t i = 0; i < p->e.size(); ++i) {
                key += ((reinterpret_cast<std::uintptr_t>(p->e[i].p) >> i) +
                        (reinterpret_cast<std::uintptr_t>(p->e[i].w.r) >> i) +
                        (reinterpret_cast<std::uintptr_t>(p->e[i].w.i) >> (i + 1))) &
                       MASK;
                key &= MASK;
            }
            return key;
        }

        // access functions
        [[nodiscard]] auto        getNodeCount() const { return nodeCount; }
        [[nodiscard]] auto        getPeakNodeCount() const { return peakNodeCount; }
        [[nodiscard]] auto        getAllocations() const { return allocations; }
        [[nodiscard]] const auto& getTables() const { return tables; }

        // lookup a node in the unique table for the appropriate variable; insert it, if it has not been found
        // NOTE: reference counting is to be adjusted by function invoking the table lookup and only normalized nodes shall be stored.
        Edge lookup(const Edge& e, bool keepNode = false) {
            // there are unique terminal nodes
            if (e.p->v == -1)
                return e;

            lookups++;
            auto key = hash(e.p);
            auto v   = e.p->v;

            // successors of a node shall either have successive variable numbers or be terminals
            for ([[maybe_unused]] const auto& edge: e.p->e)
                assert(edge.p->v == v - 1 || isTerminal(edge));

            auto p = tables[v][key];
            while (p != nullptr) {
                if (e.p->e == p->e) {
                    // Match found
                    if (e.p != p && !keepNode) {
                        // put node pointed to by e.p on available chain
                        returnNode(e.p);
                    }
                    hits++;

                    // variables should stay the same
                    assert(p->v == e.p->v);
                    // successors of a node shall either have successive variable numbers or be terminals
                    for ([[maybe_unused]] const auto& edge: e.p->e)
                        assert(edge.p->v == v - 1 || isTerminal(edge));

                    return {p, e.w};
                }
                collisions++;
                p = p->next;
            }

            // node was not found -> add it to front of available space chain
            e.p->next      = tables[v][key];
            tables[v][key] = e.p;
            nodeCount++;
            peakNodeCount = std::max(peakNodeCount, nodeCount);

            return e;
        }

        Node* getNode() {
            Node* r;
            if (available != nullptr) // get node from available chain if possible
            {
                r         = available;
                available = available->next;
            } else { // otherwise allocate new nodes
                r = new Node[ALLOCATION_SIZE];
                allocations += ALLOCATION_SIZE;
                chunks.push_back(r);
                Node* r2  = r + 1;
                available = r2;
                for (auto i = 0U; i < ALLOCATION_SIZE - 2; i++) {
                    r2->next = r2 + 1;
                    ++r2;
                }
                r2->next = nullptr;
            }
            r->next = nullptr;
            r->ref  = 0; // nodes in available space chain may have non-zero ref count
            return r;
        }

        void returnNode(Node* p) {
            p->next   = available;
            available = p;
        }

        // increment reference counter for node e points to
        // and recursively increment reference counter for
        // each child if this is the first reference
        void incRef(const Edge& e) {
            dd::ComplexNumbers::incRef(e.w);
            if (isTerminal(e))
                return;

            if (e.p->ref == std::numeric_limits<decltype(e.p->ref)>::max()) {
                std::clog << "[WARN] MAXREFCNT reached for p=" << reinterpret_cast<uintptr_t>(e.p) << ". Node will never be collected." << std::endl;
                return;
            }

            e.p->ref++;

            if (e.p->ref == 1) {
                for (const auto& edge: e.p->e) {
                    if (edge.p != nullptr) {
                        incRef(edge);
                    }
                }
                active[e.p->v]++;
                activeNodeCount++;
                maxActive = std::max(maxActive, activeNodeCount);
            }
        }

        // decrement reference counter for node e points to
        // and recursively decrement reference counter for
        // each child if this is the last reference
        void decRef(const Edge& e) {
            dd::ComplexNumbers::decRef(e.w);
            if (isTerminal(e)) return;
            if (e.p->ref == std::numeric_limits<decltype(e.p->ref)>::max()) return;

            if (e.p->ref == 0) {
                throw std::runtime_error("In decref: ref==0 before decref\n");
            }

            e.p->ref--;

            if (e.p->ref == 0) {
                for (const auto& edge: e.p->e) {
                    if (edge.p != nullptr) {
                        decRef(edge);
                    }
                }
                active[e.p->v]--;
                activeNodeCount--;
            }
        }

        std::size_t garbageCollect(bool force = false) {
            gcCalls++;
            if (!force && nodeCount < gcLimit)
                return 0;

            gcRuns++;
            std::size_t collected = 0;
            std::size_t remaining = 0;
            for (auto& table: tables) {
                for (auto& bucket: table) {
                    Node* lastp = nullptr;
                    Node* p     = bucket;
                    while (p != nullptr) {
                        if (p->ref == 0) {
                            if (p->v == -1) {
                                throw std::runtime_error("Tried to collect a terminal node.");
                            }
                            collected++;
                            auto nextp = p->next;
                            if (lastp == nullptr)
                                bucket = p->next;
                            else
                                lastp->next = p->next;
                            returnNode(p);
                            p = nextp;
                        } else {
                            lastp = p;
                            p     = p->next;
                            remaining++;
                        }
                    }
                }
            }
            gcLimit += gcIncrement;
            nodeCount = remaining;
            return collected;
        }

        void clear() {
            for (auto& table: tables) {
                for (auto& bucket: table) {
                    if (bucket == nullptr)
                        continue;

                    // successively return nodes to available space chain
                    auto current = bucket;
                    while (current != nullptr) {
                        bucket = current->next;
                        returnNode(current);
                        current = bucket;
                    }
                }
            }
            nodeCount     = 0;
            peakNodeCount = 0;

            collisions = 0;
            hits       = 0;
            lookups    = 0;

            std::fill(active.begin(), active.end(), 0);
            activeNodeCount = 0;
            maxActive       = 0;

            gcCalls = 0;
            gcRuns  = 0;
            gcLimit = gcInitialLimit;
        };

        void print() {
            Qubit q = nvars - 1;
            for (auto it = tables.rbegin(); it != tables.rend(); ++it) {
                auto& table = *it;
                std::cout << "\t" << static_cast<std::size_t>(q) << ":" << std::endl;
                for (size_t key = 0; key < table.size(); ++key) {
                    auto p = table[key];
                    if (p != nullptr)
                        std::cout << key << ": ";
                    while (p != nullptr) {
                        std::cout << "\t\t" << reinterpret_cast<uintptr_t>(p) << " " << p->ref << "\t";
                        p = p->next;
                    }
                    if (table[key] != nullptr)
                        std::cout << std::endl;
                }
                --q;
            }
        }

        void printActive() {
            std::cout << "#printActive: " << activeNodeCount << ", ";
            for (const auto& a: active)
                std::cout << a << " ";
            std::cout << std::endl;
        }

        [[nodiscard]] fp hitRatio() const { return static_cast<fp>(hits) / lookups; }
        [[nodiscard]] fp colRatio() const { return static_cast<fp>(collisions) / lookups; }

        [[nodiscard]] auto getActiveNodeCount() const {
            return activeNodeCount;
        }
        [[nodiscard]] auto getActiveNodeCount(Qubit var) { return active.at(var); }

        std::ostream& printStatistics(std::ostream& os = std::cout) {
            os << "hits: " << hits << ", collisions: " << collisions << ", looks: " << lookups << ", hitRatio: " << hitRatio() << ", colRatio: " << colRatio() << ", gc calls: " << gcCalls << ", gc runs: " << gcRuns << std::endl;
            return os;
        }

    private:
        using NodeBucket = std::array<Node*, NBUCKET>;

        // unique tables (one per input variable)
        std::size_t             nvars = 0;
        std::vector<NodeBucket> tables{nvars};

        Node*              available{};
        std::vector<Node*> chunks{};
        std::size_t        allocations   = 0;
        std::size_t        nodeCount     = 0;
        std::size_t        peakNodeCount = 0;

        // unique table lookup statistics
        std::size_t collisions = 0;
        std::size_t hits       = 0;
        std::size_t lookups    = 0;

        // (max) active nodes
        // number of active vector nodes for each variable
        std::vector<std::size_t> active{nvars};
        std::size_t              activeNodeCount = 0;
        std::size_t              maxActive       = 0;

        // garbage collection
        std::size_t gcCalls        = 0;
        std::size_t gcRuns         = 0;
        std::size_t gcInitialLimit = 250000;
        std::size_t gcLimit        = 250000;
        std::size_t gcIncrement    = 0;

        [[nodiscard]] static bool isTerminal(const Edge& e) { return e.p->v == -1; }
    };

} // namespace dd

#endif //DDpackage_UNIQUETABLE_HPP