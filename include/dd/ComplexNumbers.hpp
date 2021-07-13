/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DDcomplex_H
#define DDcomplex_H

#include "Complex.hpp"
#include "ComplexCache.hpp"
#include "ComplexValue.hpp"
#include "Definitions.hpp"
#include "dd/tables/MagnitudeTable.hpp"
#include "dd/tables/PhaseTable.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>

namespace dd {
    struct ComplexNumbers {
        ComplexCache<>   complexCache{};
        MagnitudeTable<> magnitudeTable{};
        PhaseTable<>     phaseTable{};

        ComplexNumbers()  = default;
        ~ComplexNumbers() = default;

        void clear() {
            magnitudeTable.clear();
            phaseTable.clear();
            complexCache.clear();
        }

        static void setTolerance(fp tol) {
            MagnitudeTable<>::setTolerance(tol);
            PhaseTable<>::setTolerance(tol);
        }

        // operations on complex numbers
        // meanings are self-evident from the names
        static void add(Complex& r, const Complex& a, const Complex& b) {
            assert(r != Complex::zero);
            assert(r != Complex::one);
            if (a.approximatelyZero()) {
                r.setVal(b);
            } else if (b.approximatelyZero()) {
                r.setVal(a);
            } else {
                const auto& aMag   = a.mag->value;
                const auto& aPhase = PhaseEntry::val(a.phase);
                const auto& bMag   = b.mag->value;
                const auto& bPhase = PhaseEntry::val(b.phase);

                //                std::cout << ComplexValue::toString(aMag, aPhase) << " + " << ComplexValue::toString(bMag, bPhase) << " = ";

                const auto minmax    = std::minmax(aPhase, bPhase);
                const auto diff      = minmax.second - minmax.first;
                const auto nearby    = std::nearbyint(diff);
                const auto closeness = std::abs(diff - nearby);

                if (closeness < decltype(magnitudeTable)::tolerance()) {
                    if ((static_cast<std::uint64_t>(nearby) % 2) == 0) {
                        // same phase
                        r.mag->value   = aMag + bMag;
                        r.phase->value = aPhase;
                    } else {
                        // opposing phase
                        if (std::abs(aMag - bMag) < decltype(magnitudeTable)::tolerance()) {
                            r.mag->value   = 0.;
                            r.phase->value = 0.;
                        } else if (aMag > bMag) {
                            r.mag->value   = aMag - bMag;
                            r.phase->value = aPhase;
                        } else {
                            r.mag->value   = bMag - aMag;
                            r.phase->value = bPhase;
                        }
                    }
                } else {
                    // general case
                    auto c         = std::polar(aMag, aPhase * PI) + std::polar(bMag, bPhase * PI);
                    r.mag->value   = std::abs(c);
                    r.phase->value = std::arg(c) * PIinv;
                }
                //                std::cout << r << std::endl;
            }
        }
        static void sub(Complex& r, const Complex& a, const Complex& b) {
            assert(r != Complex::zero);
            assert(r != Complex::one);
            if (b.approximatelyZero()) {
                r.setVal(a);
            } else if (a.approximatelyZero()) {
                r.setVal(b);
                r.phase->value += 1.0;
            } else {
                const auto& aMag   = a.mag->value;
                const auto& aPhase = PhaseEntry::val(a.phase);
                const auto& bMag   = b.mag->value;
                const auto& bPhase = PhaseEntry::val(b.phase);

                if (std::abs(aPhase - bPhase) < decltype(magnitudeTable)::tolerance()) {
                    r.mag->value   = aMag - bMag;
                    r.phase->value = aPhase;
                } else {
                    auto c         = std::polar(aMag, aPhase * PI) - std::polar(bMag, bPhase * PI);
                    r.mag->value   = std::abs(c);
                    r.phase->value = std::arg(c) * PIinv;
                }
            }
        }
        static void mul(Complex& r, const Complex& a, const Complex& b) {
            assert(r != Complex::zero);
            assert(r != Complex::one);
            if (a.approximatelyOne()) {
                r.setVal(b);
            } else if (b.approximatelyOne()) {
                r.setVal(a);
            } else if (a.approximatelyZero() || b.approximatelyZero()) {
                r.mag->value   = 0.;
                r.phase->value = 0.;
            } else {
                r.mag->value   = a.mag->value * b.mag->value;
                r.phase->value = PhaseEntry::val(a.phase) + PhaseEntry::val(b.phase);
            }
        }
        static void div(Complex& r, const Complex& a, const Complex& b) {
            assert(r != Complex::zero);
            assert(r != Complex::one);
            if (a.approximatelyEquals(b)) {
                r.mag->value   = 1.;
                r.phase->value = 0.;
            } else if (a.approximatelyZero()) {
                r.mag->value   = 0.;
                r.phase->value = 0.;
            } else if (b.approximatelyOne()) {
                r.setVal(a);
            } else {
                r.mag->value   = a.mag->value / b.mag->value;
                r.phase->value = PhaseEntry::val(a.phase) - PhaseEntry::val(b.phase);
            }
        }
        static inline fp mag2(const Complex& a) {
            auto mag = a.mag->value;
            return mag * mag;
        }
        static inline fp mag(const Complex& a) {
            return a.mag->value;
        }
        static inline fp arg(const Complex& a) {
            return PhaseEntry::val(a.phase) * PI;
        }
        static void conj(Complex& r, const Complex& a) {
            if (r.phase == PhaseTable<>::zeroPtr || r.phase == PhaseTable<>::onePtr) {
                r.mag   = a.mag;
                r.phase = a.phase;
                return;
            }
            r.mag               = a.mag;
            const auto newPhase = static_cast<fp>(2) - PhaseEntry::val(a.phase);
            r.phase             = PhaseEntry::getAlignedPointer(r.phase);
            r.phase->value      = newPhase;
        }
        static void neg(Complex& r, const Complex& a) {
            assert(r.phase != PhaseTable<>::zeroPtr);
            assert(r.phase != PhaseTable<>::onePtr);
            r.mag            = a.mag;
            auto newQuadrant = (PhaseEntry::getQuadrant(a.phase) + 2) % 4;
            r.phase          = PhaseEntry::setQuadrant(a.phase, newQuadrant);
        }

        inline Complex addCached(const Complex& a, const Complex& b) {
            auto c = getCached();
            add(c, a, b);
            return c;
        }

        inline Complex subCached(const Complex& a, const Complex& b) {
            auto c = getCached();
            sub(c, a, b);
            return c;
        }

        inline Complex mulCached(const Complex& a, const Complex& b) {
            auto c = getCached();
            mul(c, a, b);
            return c;
        }

        inline Complex divCached(const Complex& a, const Complex& b) {
            auto c = getCached();
            div(c, a, b);
            return c;
        }

        // lookup a complex value in the complex table; if not found add it
        Complex lookup(const Complex& c) {
            if (c == Complex::zero) {
                return Complex::zero;
            }
            if (c == Complex::one) {
                return Complex::one;
            }

            auto mag   = c.mag->value;
            auto phase = PhaseEntry::val(c.phase);
            return lookup(mag, phase);
        }
        Complex lookup(const fp& mag, const fp& phase) {
            Complex ret{};
            assert(mag >= 0.);
            ret.mag   = magnitudeTable.lookup(mag);
            ret.phase = phaseTable.lookup(phase);
            return ret;
        }
        inline Complex lookup(const ComplexValue& c) { return lookup(c.mag, c.phase); }

        // reference counting and garbage collection
        static void incRef(const Complex& c) {
            // `zero` and `one` are static and never altered
            if (c != Complex::zero && c != Complex::one) {
                MagnitudeTable<>::incRef(c.mag);
                PhaseTable<>::incRef(c.phase);
            }
        }
        static void decRef(const Complex& c) {
            // `zero` and `one` are static and never altered
            if (c != Complex::zero && c != Complex::one) {
                MagnitudeTable<>::decRef(c.mag);
                PhaseTable<>::decRef(c.phase);
            }
        }
        std::size_t garbageCollect(bool force = false) {
            return magnitudeTable.garbageCollect(force) + phaseTable.garbageCollect(force);
        }

        // provide (temporary) cached complex number
        inline Complex getTemporary() {
            return complexCache.getTemporaryComplex();
        }

        inline Complex getTemporary(const fp& mag, const fp& phase) {
            auto c         = complexCache.getTemporaryComplex();
            c.mag->value   = mag;
            c.phase->value = phase;
            return c;
        }

        inline Complex getTemporary(const ComplexValue& c) {
            return getTemporary(c.mag, c.phase);
        }

        inline Complex getCached() {
            return complexCache.getCachedComplex();
        }

        inline Complex getCached(const fp& mag, const fp& phase) {
            auto c         = complexCache.getCachedComplex();
            c.mag->value   = mag;
            c.phase->value = phase;
            return c;
        }

        inline Complex getCached(const ComplexValue& c) {
            return getCached(c.mag, c.phase);
        }

        void returnToCache(Complex& c) {
            complexCache.returnToCache(c);
        }

        [[nodiscard]] std::size_t cacheCount() const {
            return complexCache.getCount();
        }
    };
} // namespace dd
#endif
