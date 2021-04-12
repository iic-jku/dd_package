/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DDcomplex_H
#define DDcomplex_H

#include "Definitions.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace dd {
    struct ComplexTableEntry {
        fp                 val;
        ComplexTableEntry* next;
        RefCount           ref;
    };

    struct Complex {
        ComplexTableEntry* r;
        ComplexTableEntry* i;

        bool operator==(const Complex& other) const {
            return r == other.r && i == other.i;
        }

        bool operator!=(const Complex& other) const {
            return !operator==(other);
        }
    };

    struct ComplexValue {
        fp r, i;

        static ComplexValue readBinary(std::istream& is) {
            ComplexValue temp{};
            is.read(reinterpret_cast<char*>(&temp.r), sizeof(fp));
            is.read(reinterpret_cast<char*>(&temp.i), sizeof(fp));
            return temp;
        }

        static ComplexValue from_string(const std::string& real_str, std::string imag_str) {
            fp real = real_str.empty() ? 0. : std::stod(real_str);

            imag_str.erase(remove(imag_str.begin(), imag_str.end(), ' '), imag_str.end());
            imag_str.erase(remove(imag_str.begin(), imag_str.end(), 'i'), imag_str.end());
            if (imag_str == "+" || imag_str == "-") imag_str = imag_str + "1";
            fp imag = imag_str.empty() ? 0. : std::stod(imag_str);
            return ComplexValue{real, imag};
        }

        bool operator==(const ComplexValue& other) const {
            return r == other.r && i == other.i;
        }

        bool operator!=(const ComplexValue& other) const {
            return !operator==(other);
        }
    };

    class ComplexNumbers {
        struct ComplexChunk {
            ComplexTableEntry* entry;
            ComplexChunk*      next;
        };

        static inline ComplexTableEntry zeroEntry{0., nullptr, 1};
        static inline ComplexTableEntry oneEntry{1., nullptr, 1};
        static ComplexTableEntry*       moneEntryPointer;

        ComplexTableEntry* lookupVal(const fp& val);
        ComplexTableEntry* getComplexTableEntry();

    public:
        constexpr static Complex ZERO{(&zeroEntry), (&zeroEntry)};
        constexpr static Complex ONE{(&oneEntry), (&zeroEntry)};
        static constexpr fp      SQRT_2 = 0.707106781186547524400844362104849039284835937688474036588L;
        static constexpr fp      PI     = 3.141592653589793238462643383279502884197169399375105820974L;

        static constexpr unsigned short CACHE_SIZE  = 1800;
        static constexpr unsigned short CHUNK_SIZE  = 2000;
        static constexpr unsigned short NBUCKET     = 32768;
        static constexpr std::size_t    GCLIMIT     = 100000;
        static constexpr std::size_t    GCINCREMENT = 0;
        static inline fp                TOLERANCE   = 1e-13;

        std::array<ComplexTableEntry*, NBUCKET> ComplexTable{};
        ComplexTableEntry*                      Avail                       = nullptr;
        ComplexTableEntry*                      Cache_Avail                 = nullptr;
        ComplexTableEntry*                      Cache_Avail_Initial_Pointer = nullptr;
        ComplexChunk*                           chunks                      = nullptr;

        std::size_t count      = 0;
        long        cacheCount = CACHE_SIZE;
        std::size_t ctCalls    = 0;
        std::size_t ctMiss     = 0;
        std::size_t gcCalls    = 0;
        std::size_t gcRuns     = 0;
        std::size_t gcLimit    = GCLIMIT;

        ComplexNumbers();
        ~ComplexNumbers();

        // TODO: this does not properly clear the whole complex table
        void clear() {
            count      = 0;
            cacheCount = CACHE_SIZE;

            ctCalls = 0;
            ctMiss  = 0;

            gcCalls = 0;
            gcRuns  = 0;
            gcLimit = GCLIMIT;
        }

        static void setTolerance(fp tol) {
            TOLERANCE = tol;
        }

        // operations on complex numbers
        // meanings are self-evident from the names
        static inline fp val(const ComplexTableEntry* x) {
            if (isNegativePointer(x)) {
                return -getAlignedPointer(x)->val;
            }
            return x->val;
        }
        /**
         * The pointer to ComplexTableEntry may encode the sign in the least significant bit, causing mis-aligned access
         * if not handled properly.
         * @param entry pointer to ComplexTableEntry possibly unsafe to deference
         * @return safe pointer for deferencing
         */
        static inline ComplexTableEntry* getAlignedPointer(const ComplexTableEntry* entry) {
            return reinterpret_cast<ComplexTableEntry*>(reinterpret_cast<std::uintptr_t>(entry) & (~1ULL));
        }
        static inline ComplexTableEntry* getNegativePointer(const ComplexTableEntry* entry) {
            return reinterpret_cast<ComplexTableEntry*>(reinterpret_cast<std::uintptr_t>(entry) | 1ULL);
        }
        static inline ComplexTableEntry* flipPointerSign(const ComplexTableEntry* entry) {
            return reinterpret_cast<ComplexTableEntry*>(reinterpret_cast<std::uintptr_t>(entry) ^ 1ULL);
        }
        static inline void setNegativePointer(ComplexTableEntry*& entry) {
            entry = reinterpret_cast<ComplexTableEntry*>(reinterpret_cast<std::uintptr_t>(entry) | 1ULL);
        }
        static inline bool isNegativePointer(const ComplexTableEntry* entry) {
            return reinterpret_cast<std::uintptr_t>(entry) & 1ULL;
        }

        static inline bool equals(const Complex& a, const Complex& b) {
            return std::abs(val(a.r) - val(b.r)) < TOLERANCE && std::abs(val(a.i) - val(b.i)) < TOLERANCE;
        };

        static inline bool equals(const ComplexValue& a, const ComplexValue& b) {
            return std::abs(a.r - b.r) < TOLERANCE && std::abs(a.i - b.i) < TOLERANCE;
        }

        static inline bool equalsZero(const Complex& c) {
            return c == ZERO || (std::abs(val(c.r)) < TOLERANCE && std::abs(val(c.i)) < TOLERANCE);
        }
        static inline bool equalsZero(const ComplexValue& c) {
            return std::abs(c.r) < TOLERANCE && std::abs(c.i) < TOLERANCE;
        }

        static inline bool equalsOne(const Complex& c) {
            return c == ONE || (std::abs(val(c.r) - 1) < TOLERANCE && std::abs(val(c.i)) < TOLERANCE);
        }
        static inline bool equalsOne(const ComplexValue& c) {
            return std::abs(c.r - 1) < TOLERANCE && std::abs(c.i) < TOLERANCE;
        }

        static void add(Complex& r, const Complex& a, const Complex& b);
        static void sub(Complex& r, const Complex& a, const Complex& b);
        static void mul(Complex& r, const Complex& a, const Complex& b);
        static void div(Complex& r, const Complex& a, const Complex& b);

        static inline fp mag2(const Complex& a) {
            auto ar = val(a.r);
            auto ai = val(a.i);

            return ar * ar + ai * ai;
        }
        static inline fp mag(const Complex& a) {
            return std::sqrt(mag2(a));
        }
        static inline fp arg(const Complex& a) {
            auto ar = val(a.r);
            auto ai = val(a.i);
            return std::atan2(ai, ar);
        }
        static Complex conj(const Complex& a);
        static Complex neg(const Complex& a);

        inline Complex addCached(const Complex& a, const Complex& b) {
            auto c = getCachedComplex();
            add(c, a, b);
            return c;
        }

        inline Complex subCached(const Complex& a, const Complex& b) {
            auto c = getCachedComplex();
            sub(c, a, b);
            return c;
        }

        inline Complex mulCached(const Complex& a, const Complex& b) {
            auto c = getCachedComplex();
            mul(c, a, b);
            return c;
        }

        inline Complex divCached(const Complex& a, const Complex& b) {
            auto c = getCachedComplex();
            div(c, a, b);
            return c;
        }

        inline void releaseCached(const Complex& c) {
            assert(c != ZERO);
            assert(c != ONE);
            assert(c.r->ref == 0);
            assert(c.i->ref == 0);
            c.i->next   = Cache_Avail;
            Cache_Avail = c.r;
            cacheCount += 2;
            assert(cacheCount <= CACHE_SIZE);
        }

        static inline unsigned long getKey(const fp& val) {
            assert(val >= 0);
            //const fp hash = 4*val*val*val - 6*val*val + 3*val; // cubic
            const fp hash = val; // linear
            auto     key  = static_cast<unsigned long>(hash * (NBUCKET - 1));
            if (key > NBUCKET - 1) {
                key = NBUCKET - 1;
            }
            return key;
        };

        // lookup a complex value in the complex value table; if not found add it
        Complex lookup(const Complex& c);

        Complex lookup(const fp& r, const fp& i);

        inline Complex lookup(const ComplexValue& c) { return lookup(c.r, c.i); }

        // reference counting and garbage collection
        static void incRef(const Complex& c);
        static void decRef(const Complex& c);
        std::size_t garbageCollect(bool force = false);

        // provide (temporary) cached complex number
        inline Complex getTempCachedComplex() {
            assert(cacheCount >= 2);
            return {Cache_Avail, Cache_Avail->next};
        }

        inline Complex getTempCachedComplex(const fp& r, const fp& i) {
            assert(cacheCount >= 2);
            Cache_Avail->val       = r;
            Cache_Avail->next->val = i;
            return {Cache_Avail, Cache_Avail->next};
        }

        inline Complex getTempCachedComplex(const ComplexValue& c) {
            return getTempCachedComplex(c.r, c.i);
        }

        inline Complex getCachedComplex() {
            assert(cacheCount >= 2);
            cacheCount -= 2;
            Complex c{Cache_Avail, Cache_Avail->next};
            Cache_Avail = Cache_Avail->next->next;
            return c;
        }

        inline Complex getCachedComplex(const fp& r, const fp& i) {
            assert(cacheCount >= 2);
            cacheCount -= 2;
            Complex c{Cache_Avail, Cache_Avail->next};
            c.r->val    = r;
            c.i->val    = i;
            Cache_Avail = Cache_Avail->next->next;
            return c;
        }

        inline Complex getCachedComplex(const ComplexValue& c) {
            return getCachedComplex(c.r, c.i);
        }

        // printing
        static void printFormattedReal(std::ostream& os, fp r, bool imaginary = false);
        void        printComplexTable() const;
        void        statistics() const;

        [[nodiscard]] int cacheSize() const;

        static std::string toString(const Complex& c, bool formatted = true, int precision = -1);
    };

    inline ComplexTableEntry* ComplexNumbers::moneEntryPointer{getNegativePointer(&oneEntry)};

    std::ostream& operator<<(std::ostream& os, const Complex& c);

    std::ostream& operator<<(std::ostream& os, const ComplexValue& c);
} // namespace dd

namespace std {
    template<>
    struct hash<dd::Complex> {
        std::size_t operator()(dd::Complex const& c) const noexcept {
            auto h1 = std::hash<dd::ComplexTableEntry*>{}(c.r);
            auto h2 = std::hash<dd::ComplexTableEntry*>{}(c.i);
            return h1 ^ (h2 << 1);
        }
    };

    template<>
    struct hash<dd::ComplexValue> {
        std::size_t operator()(dd::ComplexValue const& c) const noexcept {
            auto h1 = std::hash<dd::fp>{}(c.r);
            auto h2 = std::hash<dd::fp>{}(c.i);
            return h1 ^ (h2 << 1);
        }
    };
} // namespace std
#endif