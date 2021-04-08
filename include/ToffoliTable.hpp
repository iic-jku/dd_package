/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DDpackage_TOFFOLITABLE_HPP
#define DDpackage_TOFFOLITABLE_HPP

#include "Definitions.hpp"

#include <array>
#include <set>

namespace dd {
    template<class Edge, std::size_t NBUCKET = 2048>
    class ToffoliTable {
    public:
        ToffoliTable() = default;

        struct Entry {
            QubitCount        n = 0;
            std::set<Control> controls{};
            Qubit             target = 0;
            Edge              e;
        };

        static constexpr size_t MASK = NBUCKET - 1;

        // access functions
        [[nodiscard]] const auto& getTable() const { return table; }

        void insert(QubitCount n, const std::set<Control>& controls, Qubit target, const Edge& e) {
            const auto key = hash(controls, target);
            table[key]     = {.n = n, .controls = controls, .target = target, .e = e};
            ++count;
        }

        Edge lookup(QubitCount n, const std::set<Control>& controls, Qubit target) {
            lookups++;
            Edge        r{};
            const auto  key   = hash(controls, target);
            const auto& entry = table[key];
            if (entry.e.p == nullptr) return r;
            if (entry.n != n) return r;
            if (entry.target != target) return r;
            if (entry.controls != controls) return r;
            hits++;
            return entry.e;
        }

        static size_t hash(const std::set<Control>& controls, Qubit target) {
            auto key = static_cast<std::size_t>(std::make_unsigned<Qubit>::type(target));
            for (const auto& control: controls) {
                if (control.type == dd::Control::Type::pos) {
                    key *= 29u * control.qubit;
                } else {
                    key *= 71u * control.qubit;
                }
            }
            return key & MASK;
        }

        void clear() {
            if (count > 0) {
                for (auto& entry: table)
                    entry.e.p = nullptr;
                count = 0;
            }
            hits    = 0;
            lookups = 0;
        }

        [[nodiscard]] fp hitRatio() const { return static_cast<fp>(hits) / lookups; }
        std::ostream&    printStatistics(std::ostream& os = std::cout) {
            os << "hits: " << hits << ", looks: " << lookups << ", ratio: " << hitRatio() << std::endl;
            return os;
        }

    private:
        std::array<Entry, NBUCKET> table{};

        // compute table lookup statistics
        std::size_t hits    = 0;
        std::size_t lookups = 0;
        std::size_t count   = 0;
    };
} // namespace dd

#endif //DDpackage_TOFFOLITABLE_HPP
