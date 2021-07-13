/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_COMPLEXVALUE_HPP
#define DD_PACKAGE_COMPLEXVALUE_HPP

#include "Definitions.hpp"
#include "dd/tables/MagnitudeTable.hpp"
#include "dd/tables/PhaseTable.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

namespace dd {
    struct ComplexValue {
        fp mag, phase;

        [[nodiscard]] inline bool approximatelyEquals(const ComplexValue& c) const {
            return std::abs(mag - c.mag) < MagnitudeTable<>::tolerance() &&
                   std::abs(phase - c.phase) < PhaseTable<>::tolerance();
        }

        [[nodiscard]] inline bool approximatelyZero() const {
            return std::abs(mag) < MagnitudeTable<>::tolerance();
        }

        [[nodiscard]] inline bool approximatelyOne() const {
            return std::abs(mag - 1.) < MagnitudeTable<>::tolerance() &&
                   std::abs(phase) < PhaseTable<>::tolerance();
        }

        inline bool operator==(const ComplexValue& other) const {
            return mag == other.mag && phase == other.phase;
        }

        inline bool operator!=(const ComplexValue& other) const {
            return !operator==(other);
        }

        void readBinary(std::istream& is) {
            is.read(reinterpret_cast<char*>(&mag), sizeof(decltype(mag)));
            is.read(reinterpret_cast<char*>(&phase), sizeof(decltype(phase)));
        }

        void writeBinary(std::ostream& os) const {
            os.write(reinterpret_cast<const char*>(&mag), sizeof(decltype(mag)));
            os.write(reinterpret_cast<const char*>(&phase), sizeof(decltype(phase)));
        }

        void from_string(const std::string& mag_str, const std::string& phase_str) {
            mag   = mag_str.empty() ? 0. : std::stod(mag_str);
            phase = phase_str.empty() ? 0. : std::stod(phase_str);
        }

        static auto getLowestFraction(const double x, const std::uint64_t maxDenominator = std::uint64_t(1) << 10, const fp tol = dd::MagnitudeTable<>::tolerance()) {
            assert(x >= 0.);

            std::pair<std::uint64_t, std::uint64_t> lowerBound{0U, 1U};
            std::pair<std::uint64_t, std::uint64_t> upperBound{1U, 0U};

            while ((lowerBound.second <= maxDenominator) && (upperBound.second <= maxDenominator)) {
                auto num    = lowerBound.first + upperBound.first;
                auto den    = lowerBound.second + upperBound.second;
                auto median = static_cast<fp>(num) / static_cast<fp>(den);
                if (std::abs(x - median) < tol) {
                    if (den <= maxDenominator) {
                        return std::pair{num, den};
                    } else if (upperBound.second > lowerBound.second) {
                        return upperBound;
                    } else {
                        return lowerBound;
                    }
                } else if (x > median) {
                    lowerBound = {num, den};
                } else {
                    upperBound = {num, den};
                }
            }
            if (lowerBound.second > maxDenominator) {
                return upperBound;
            } else {
                return lowerBound;
            }
        }

        static void printFormatted(std::ostream& os, fp r, bool phase = false) {
            const auto tol = dd::MagnitudeTable<>::tolerance();
            if (std::abs(r) < tol) {
                return;
            }
            if (phase) {
                // special case treatment for +-i
                if (std::abs(r - 0.5) < tol) {
                    os << "i";
                    return;
                } else if (std::abs(r + 0.5) < tol) {
                    os << "-i";
                    return;
                }

                os << "ℯ(";
                if (std::signbit(r))
                    os << "-";
                os << "iπ";
            }

            auto absr = std::abs(r);

            if (std::abs(absr - 1.) < tol) { // +-1
                if (phase) {
                    os << ")";
                } else {
                    os << (std::signbit(r) ? "-" : "") << 1;
                }
                return;
            }

            auto fraction = getLowestFraction(absr);
            auto approx   = static_cast<fp>(fraction.first) / static_cast<fp>(fraction.second);
            auto error    = std::abs(absr - approx);

            if (error < tol) { // suitable fraction a/b found
                if (phase) {
                    if (fraction.second == 1U) {
                        os << " " << fraction.first << ")";
                    } else if (fraction.first == 1U) {
                        os << "/" << fraction.second << ")";
                    } else {
                        os << " " << fraction.first << "/" << fraction.second << ")";
                    }
                } else {
                    if (fraction.second == 1U) {
                        os << fraction.first;
                    } else {
                        os << fraction.first << "/" << fraction.second;
                    }
                }
                return;
            }

            auto abssqrt = absr / SQRT2_2;

            if (std::abs(abssqrt - 1.) < tol) { // +- 1/sqrt(2)
                if (phase) {
                    os << "/√2)";
                } else {
                    os << (std::signbit(r) ? "-" : "") << "1/√2";
                }
                return;
            }

            fraction = getLowestFraction(abssqrt);
            approx   = static_cast<fp>(fraction.first) / static_cast<fp>(fraction.second);
            error    = std::abs(abssqrt - approx);

            if (error < tol) { // suitable fraction a/(b * sqrt(2)) found
                if (phase) {
                    if (fraction.second == 1U) {
                        os << " " << fraction.first << "/√2)";
                    } else if (fraction.first == 1U) {
                        os << "/(" << fraction.second << "√2))";
                    } else {
                        os << " " << fraction.first << "/(" << fraction.second << "√2))";
                    }
                } else {
                    if (fraction.second == 1U) {
                        os << fraction.first << "/√2";
                    } else {
                        os << fraction.first << "/(" << fraction.second << "√2)";
                    }
                }
                return;
            }

            auto abspi = absr / PI;

            if (std::abs(abspi - 1.) < tol) { // +- π
                if (phase) {
                    os << " π)";
                } else {
                    os << (std::signbit(r) ? "-" : "") << "π";
                }
                return;
            }

            fraction = getLowestFraction(abspi);
            approx   = static_cast<fp>(fraction.first) / static_cast<fp>(fraction.second);
            error    = std::abs(abspi - approx);

            if (error < tol) { // suitable fraction a/b π found
                if (phase) {
                    if (fraction.second == 1U) {
                        os << " " << fraction.first << "π)";
                    } else if (fraction.first == 1U) {
                        os << " π/" << fraction.second << ")";
                    } else {
                        os << " " << fraction.first << "π/" << fraction.second << ")";
                    }
                } else {
                    if (fraction.second == 1U) {
                        os << fraction.first << "π";
                    } else if (fraction.first == 1U) {
                        os << "π/" << fraction.second;
                    } else {
                        os << fraction.first << "π/" << fraction.second;
                    }
                }
                return;
            }

            if (phase) { // default
                os << " " << absr << ")";
            } else
                os << r;
        }

        static std::string toString(const fp& mag, const fp& phase, bool formatted = true, int precision = -1) {
            std::ostringstream ss{};

            if (precision >= 0) ss << std::setprecision(precision);

            if (std::abs(mag) < MagnitudeTable<>::tolerance()) {
                return "0";
            }

            if (formatted) {
                if (std::abs(std::abs(phase) - 1.0) < PhaseTable<>::tolerance()) {
                    ss << "-";
                    ComplexValue::printFormatted(ss, mag);
                } else {
                    if (std::abs(mag - 1.0) > MagnitudeTable<>::tolerance()) {
                        ComplexValue::printFormatted(ss, mag);
                        if (std::abs(phase) > PhaseTable<>::tolerance()) {
                            ss << " ";
                        }
                    } else if (std::abs(phase) < PhaseTable<>::tolerance()) {
                        ss << "1";
                        return ss.str();
                    }
                    ComplexValue::printFormatted(ss, phase, true);
                }
            } else {
                ss << mag;
                if (std::abs(phase) > PhaseTable<>::tolerance()) {
                    ss << " " << phase;
                }
            }
            return ss.str();
        }

        explicit operator auto() const { return std::polar(mag, phase * PI); }

        ComplexValue& operator+=(const ComplexValue& rhs) {
            auto a = std::polar(mag, phase * PI);
            auto b = std::polar(rhs.mag, rhs.phase * PI);
            a += b;

            mag   = std::abs(a);
            phase = std::arg(a) / PI;
            return *this;
        }

        friend ComplexValue operator+(ComplexValue lhs, const ComplexValue& rhs) {
            lhs += rhs;
            return lhs;
        }
    };

    inline std::ostream& operator<<(std::ostream& os, const ComplexValue& c) {
        os << ComplexValue::toString(c.mag, std::remainder(c.phase, 2.0));
        return os;
    }
} // namespace dd

namespace std {
    template<>
    struct hash<dd::ComplexValue> {
        std::size_t operator()(dd::ComplexValue const& c) const noexcept {
            auto h1 = dd::murmur64(static_cast<std::size_t>(std::round(c.mag / dd::MagnitudeTable<>::tolerance())));
            auto h2 = dd::murmur64(static_cast<std::size_t>(std::round(c.phase / dd::PhaseTable<>::tolerance())));
            return dd::combineHash(h1, h2);
        }
    };
} // namespace std
#endif //DD_PACKAGE_COMPLEXVALUE_HPP
