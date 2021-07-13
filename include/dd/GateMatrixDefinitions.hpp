/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_GATEMATRIXDEFINITIONS_H
#define DD_PACKAGE_GATEMATRIXDEFINITIONS_H

#include "ComplexValue.hpp"
#include "Definitions.hpp"

#include <array>
#include <cmath>

namespace dd {
    // Complex constants
    constexpr ComplexValue complex_one       = {1., 0.};
    constexpr ComplexValue complex_mone      = {1., 1.};
    constexpr ComplexValue complex_zero      = {0., 0.};
    constexpr ComplexValue complex_i         = {1., 0.5};
    constexpr ComplexValue complex_mi        = {1., -0.5};
    constexpr ComplexValue complex_SQRT2_2   = {SQRT2_2, 0.};
    constexpr ComplexValue complex_mSQRT2_2  = {SQRT2_2, 1.};
    constexpr ComplexValue complex_iSQRT2_2  = {SQRT2_2, 0.5};
    constexpr ComplexValue complex_miSQRT2_2 = {SQRT2_2, -0.5};
    constexpr ComplexValue complex_1plusi    = {1.0, 0.25};
    constexpr ComplexValue complex_1minusi   = {1.0, -0.25};
    constexpr ComplexValue complex_1plusi_2  = {SQRT2_2, 0.25};
    constexpr ComplexValue complex_1minusi_2 = {SQRT2_2, -0.25};

    // Gate matrices
    using GateMatrix = std::array<ComplexValue, NEDGE>;
    constexpr GateMatrix Imat{complex_one, complex_zero, complex_zero, complex_one};
    constexpr GateMatrix Hmat{complex_SQRT2_2, complex_SQRT2_2, complex_SQRT2_2, complex_mSQRT2_2};
    constexpr GateMatrix Xmat{complex_zero, complex_one, complex_one, complex_zero};
    constexpr GateMatrix Ymat{complex_zero, complex_mi, complex_i, complex_zero};
    constexpr GateMatrix Zmat{complex_one, complex_zero, complex_zero, complex_mone};
    constexpr GateMatrix Smat{complex_one, complex_zero, complex_zero, complex_i};
    constexpr GateMatrix Sdagmat{complex_one, complex_zero, complex_zero, complex_mi};
    constexpr GateMatrix Tmat{complex_one, complex_zero, complex_zero, complex_1plusi};
    constexpr GateMatrix Tdagmat{complex_one, complex_zero, complex_zero, complex_1minusi};
    constexpr GateMatrix SXmat{complex_1plusi_2, complex_1minusi_2, complex_1minusi_2, complex_1plusi_2};
    constexpr GateMatrix SXdagmat{complex_1minusi_2, complex_1plusi_2, complex_1plusi_2, complex_1minusi_2};
    constexpr GateMatrix Vmat{complex_SQRT2_2, complex_miSQRT2_2, complex_miSQRT2_2, complex_SQRT2_2};
    constexpr GateMatrix Vdagmat{complex_SQRT2_2, complex_iSQRT2_2, complex_iSQRT2_2, complex_SQRT2_2};

    inline GateMatrix U3mat(fp lambda, fp phi, fp theta) {
        auto cos_2 = std::cos(theta / 2.);
        auto sin_2 = std::sin(theta / 2.);

        auto cosIncrement = 0.0;
        if (cos_2 < 0) {
            cosIncrement = 1.0;
            cos_2 *= -1;
        }

        auto sinIncrement = 0.0;
        if (sin_2 < 0) {
            sinIncrement = 1.0;
            sin_2 *= -1;
        }

        if (lambda >= 0.0) {
            return GateMatrix{{{cos_2, cosIncrement},
                               {sin_2, lambda / PI - 1.0 + sinIncrement},
                               {sin_2, phi / PI + sinIncrement},
                               {cos_2, (lambda + phi) / PI + cosIncrement}}};
        } else {
            return GateMatrix{{{cos_2, cosIncrement},
                               {sin_2, lambda / PI + 1.0 + sinIncrement},
                               {sin_2, phi / PI + sinIncrement},
                               {cos_2, (lambda + phi) / PI + cosIncrement}}};
        }
    }

    inline GateMatrix U2mat(fp lambda, fp phi) {
        if (lambda >= 0.0) {
            return GateMatrix{complex_SQRT2_2,
                              {SQRT2_2, lambda / PI - 1.0},
                              {SQRT2_2, phi / PI},
                              {SQRT2_2, (lambda + phi) / PI}};
        } else {
            return GateMatrix{complex_SQRT2_2,
                              {SQRT2_2, lambda / PI + 1.0},
                              {SQRT2_2, phi / PI},
                              {SQRT2_2, (lambda + phi) / PI}};
        }
    }

    inline GateMatrix Phasemat(fp lambda) {
        return GateMatrix{complex_one, complex_zero, complex_zero, {1.0, lambda / PI}};
    }

    inline GateMatrix RXmat(fp lambda) {
        auto cos_2 = std::cos(lambda / 2.);
        auto sin_2 = std::sin(lambda / 2.);

        auto cosIncrement = 0.0;
        if (cos_2 < 0) {
            cosIncrement = 1.0;
            cos_2 *= -1;
        }

        auto sinIncrement = 0.0;
        if (sin_2 < 0) {
            sinIncrement = 1.0;
            sin_2 *= -1;
        }
        return GateMatrix{{{cos_2, cosIncrement},
                           {sin_2, -0.5 + sinIncrement},
                           {sin_2, -0.5 + sinIncrement},
                           {cos_2, cosIncrement}}};
    }

    inline GateMatrix RYmat(fp lambda) {
        auto cos_2 = std::cos(lambda / 2.);
        auto sin_2 = std::sin(lambda / 2.);

        auto cosIncrement = 0.0;
        if (cos_2 < 0) {
            cosIncrement = 1.0;
            cos_2 *= -1;
        }

        auto sinIncrement = 0.0;
        if (sin_2 < 0) {
            sinIncrement = 1.0;
            sin_2 *= -1;
        }
        return GateMatrix{{{cos_2, cosIncrement},
                           {sin_2, 1. + sinIncrement},
                           {sin_2, sinIncrement},
                           {cos_2, cosIncrement}}};
    }

    inline GateMatrix RZmat(fp lambda) {
        return GateMatrix{{{1.0, -(lambda / 2.) / PI},
                           complex_zero,
                           complex_zero,
                           {1.0, (lambda / 2.) / PI}}};
    }
} // namespace dd
#endif //DD_PACKAGE_GATEMATRIXDEFINITIONS_H
