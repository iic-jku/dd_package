/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#include "DDpackage.hpp"

#include <cstdlib>
#include <iomanip>

namespace dd {

    Package::vNode Package::vTerminal{
            .next = nullptr,
            .e    = {{{nullptr, CN::ZERO}, {nullptr, CN::ZERO}}},
            .ref  = 0,
            .v    = -1};
    Package::mNode Package::mTerminal{
            .next  = nullptr,
            .e     = {{{nullptr, CN::ZERO}, {nullptr, CN::ZERO}, {nullptr, CN::ZERO}, {nullptr, CN::ZERO}}},
            .ref   = 0,
            .v     = -1,
            .symm  = true,
            .ident = true};

    void Package::printInformation() {
        std::cout << "\n  compiled: " << __DATE__ << " " << __TIME__
                  << "\n  Complex size: " << sizeof(Complex) << " bytes (aligned " << alignof(Complex) << " bytes)"
                  << "\n  ComplexValue size: " << sizeof(ComplexValue) << " bytes (aligned " << alignof(ComplexValue) << " bytes)"
                  << "\n  ComplexNumbers size: " << sizeof(ComplexNumbers) << " bytes (aligned " << alignof(ComplexNumbers) << " bytes)"
                  << "\n  vEdge size: " << sizeof(vEdge) << " bytes (aligned " << alignof(vEdge) << " bytes)"
                  << "\n  vNode size: " << sizeof(vNode) << " bytes (aligned " << alignof(vNode) << " bytes)"
                  << "\n  mEdge size: " << sizeof(mEdge) << " bytes (aligned " << alignof(mEdge) << " bytes)"
                  << "\n  mNode size: " << sizeof(mNode) << " bytes (aligned " << alignof(mNode) << " bytes)"
                  << "\n  CT Vector Add size: " << sizeof(decltype(vectorAdd)::Entry) << " bytes (aligned " << alignof(decltype(vectorAdd)::Entry) << " bytes)"
                  << "\n  CT Matrix Add size: " << sizeof(decltype(matrixAdd)::Entry) << " bytes (aligned " << alignof(decltype(matrixAdd)::Entry) << " bytes)"
                  << "\n  CT Matrix Transpose size: " << sizeof(decltype(matrixTranspose)::Entry) << " bytes (aligned " << alignof(decltype(matrixTranspose)::Entry) << " bytes)"
                  << "\n  CT Conjugate Matrix Transpose size: " << sizeof(decltype(conjugateMatrixTranspose)::Entry) << " bytes (aligned " << alignof(decltype(conjugateMatrixTranspose)::Entry) << " bytes)"
                  << "\n  CT Matrix Multiplication size: " << sizeof(decltype(matrixMultiplication)::Entry) << " bytes (aligned " << alignof(decltype(matrixMultiplication)::Entry) << " bytes)"
                  << "\n  CT Matrix Vector Multiplication size: " << sizeof(decltype(matrixVectorMultiplication)::Entry) << " bytes (aligned " << alignof(decltype(matrixVectorMultiplication)::Entry) << " bytes)"
                  << "\n  CT Vector Inner Product size: " << sizeof(decltype(vectorInnerProduct)::Entry) << " bytes (aligned " << alignof(decltype(vectorInnerProduct)::Entry) << " bytes)"
                  << "\n  CT Vector Kronecker size: " << sizeof(decltype(vectorKronecker)::Entry) << " bytes (aligned " << alignof(decltype(vectorKronecker)::Entry) << " bytes)"
                  << "\n  CT Matrix Kronecker size: " << sizeof(decltype(matrixKronecker)::Entry) << " bytes (aligned " << alignof(decltype(matrixKronecker)::Entry) << " bytes)"
                  << "\n  ToffolieTable::Entry size: " << sizeof(ToffoliTable<mEdge>::Entry) << " bytes (aligned " << alignof(ToffoliTable<mEdge>::Entry) << " bytes)"
                  << "\n  OperationTable::Entry size: " << sizeof(OperationTable<mEdge>::Entry) << " bytes (aligned " << alignof(OperationTable<mEdge>::Entry) << " bytes)"
                  << "\n  Package size: " << sizeof(Package) << " bytes (aligned " << alignof(Package) << " bytes)"
                  << "\n  max variables: " << MAXN
                  << "\n  UniqueTable buckets: " << NBUCKET
                  << "\n  ComputeTable slots: " << CTSLOTS
                  << "\n  ToffoliTable slots: " << TTSLOTS
                  << "\n  OperationTable slots: " << OperationSLOTS
                  << "\n  garbage collection limit: " << GCLIMIT
                  << "\n  garbage collection increment: " << GCINCREMENT
                  << "\n"
                  << std::flush;
    }

    void Package::statistics() {
        std::cout << "DD statistics:" << std::endl
                  << "[vUniqueTable] ";
        vUniqueTable.printStatistics();
        std::cout << "[mUniqueTable] ";
        mUniqueTable.printStatistics();
        std::cout << "[CT Vector Add] ";
        vectorAdd.printStatistics();
        std::cout << "[CT Matrix Add] ";
        matrixAdd.printStatistics();
        std::cout << "[CT Matrix Transpose] ";
        matrixTranspose.printStatistics();
        std::cout << "[CT Conjugate Matrix Transpose] ";
        conjugateMatrixTranspose.printStatistics();
        std::cout << "[CT Matrix Multiplication] ";
        matrixMultiplication.printStatistics();
        std::cout << "[CT Matrix Vector Multiplication] ";
        matrixVectorMultiplication.printStatistics();
        std::cout << "[CT Inner Product] ";
        vectorInnerProduct.printStatistics();
        std::cout << "[CT Vector Kronecker] ";
        vectorKronecker.printStatistics();
        std::cout << "[CT Matrix Kronecker] ";
        matrixKronecker.printStatistics();
        std::cout << "[Toffoli Table] ";
        toffoliTable.printStatistics();
        std::cout << "[Operation Table] ";
        operations.printStatistics();
    }

    Package::vEdge Package::makeVectorNode(const Qubit var, const std::array<vEdge, RADIX>& edge, bool cached) {
        vEdge e{vUniqueTable.getNode(), CN::ONE};
        e.p->v = var;
        e.p->e = edge;
        assert(e.p->ref == 0);
        assert(var - 1 == edge[0].p->v || isTerminal(edge[0]));
        assert(var - 1 == edge[1].p->v || isTerminal(edge[1]));

        // normalize it
        e = normalize(e, cached);
        assert(e.p->v == var || isTerminal(e));

        // look it up in the unique tables
        e = vUniqueTable.lookup(e, false);
        assert(e.p->v == var || isTerminal(e));

        return e;
    }

    Package::vEdge Package::normalize(const Package::vEdge& e, bool cached) {
        auto argmax = -1;

        auto zero = std::array{CN::equalsZero(e.p->e[0].w), CN::equalsZero(e.p->e[1].w)};

        for (auto i = 0U; i < RADIX; i++) {
            if (zero[i] && e.p->e[i].w != CN::ZERO) {
                cn.releaseCached(e.p->e[i].w);
                e.p->e[i] = vZero;
            }
        }

        fp sum = 0.;
        fp div = 0.;
        for (auto i = 0U; i < RADIX; ++i) {
            if (e.p->e[i].p == nullptr || zero[i]) {
                continue;
            }

            if (argmax == -1) {
                argmax = static_cast<decltype(argmax)>(i);
                div    = ComplexNumbers::mag2(e.p->e[i].w);
                sum    = div;
            } else {
                sum += ComplexNumbers::mag2(e.p->e[i].w);
            }
        }

        if (argmax == -1) {
            if (cached) {
                for (auto& i: e.p->e) {
                    if (i.p == nullptr && i.w != CN::ZERO) {
                        cn.releaseCached(i.w);
                    }
                }
            } else if (e.p != vTerminalNode) {
                // If it is not a cached variable, I have to put it pack into the chain
                vUniqueTable.returnNode(e.p);
            }
            return vZero;
        }

        sum = std::sqrt(sum / div);

        auto  r   = e;
        auto& max = r.p->e[argmax];
        if (cached && max.w != CN::ONE) {
            r.w = max.w;
            r.w.r->val *= sum;
            r.w.i->val *= sum;
        } else {
            r.w = cn.lookup(ComplexNumbers::val(max.w.r) * sum, ComplexNumbers::val(max.w.i) * sum);
            if (CN::equalsZero(r.w)) {
                return vZero;
            }
        }
        max.w = cn.lookup(static_cast<fp>(1.0) / sum, 0.);
        if (max.w == CN::ZERO)
            max = vZero;

        auto  argmin = (argmax + 1) % 2;
        auto& min    = r.p->e[argmin];
        if (!zero[argmin]) {
            if (cached) {
                cn.releaseCached(min.w);
                CN::div(min.w, min.w, r.w);
                min.w = cn.lookup(min.w);
                if (min.w == CN::ZERO) {
                    min = vZero;
                }
            } else {
                auto c = cn.getTempCachedComplex();
                CN::div(c, min.w, r.w);
                min.w = cn.lookup(c);
                if (min.w == CN::ZERO) {
                    min = vZero;
                }
            }
        }

        return r;
    }

    Package::vEdge Package::makeZeroState(Qubit mostSignificantQubit) {
        auto f = vOne;
        for (std::size_t p = 0; p <= std::make_unsigned_t<Qubit>(mostSignificantQubit); p++) {
            f = makeVectorNode(static_cast<Qubit>(p), {f, vZero});
        }
        return f;
    }

    Package::vEdge Package::makeBasisState(Qubit mostSignificantQubit, const std::bitset<MAXN>& state) {
        auto f = vOne;
        for (std::size_t p = 0; p <= std::make_unsigned_t<Qubit>(mostSignificantQubit); ++p) {
            if (state[p] == 0) {
                f = makeVectorNode(static_cast<Qubit>(p), {f, vZero});
            } else {
                f = makeVectorNode(static_cast<Qubit>(p), {vZero, f});
            }
        }
        return f;
    }

    Package::vEdge Package::makeBasisState(Qubit mostSignificantQubit, const std::vector<BasisStates>& state) {
        if (state.size() < static_cast<std::size_t>(mostSignificantQubit + 1)) {
            throw std::invalid_argument("Insufficient qubit states provided. Requested " + std::to_string(mostSignificantQubit + 1) + ", but received " + std::to_string(state.size()));
        }

        auto f = vOne;
        for (std::size_t p = 0; p <= std::make_unsigned_t<Qubit>(mostSignificantQubit); ++p) {
            switch (state[p]) {
                case BasisStates::zero:
                    f = makeVectorNode(static_cast<Qubit>(p), {f, vZero});
                    break;
                case BasisStates::one:
                    f = makeVectorNode(static_cast<Qubit>(p), {vZero, f});
                    break;
                case BasisStates::plus:
                    f = makeVectorNode(static_cast<Qubit>(p), {{{f.p, cn.lookup(CN::SQRT_2, 0)}, {f.p, cn.lookup(CN::SQRT_2, 0)}}});
                    break;
                case BasisStates::minus:
                    f = makeVectorNode(static_cast<Qubit>(p), {{{f.p, cn.lookup(CN::SQRT_2, 0)}, {f.p, cn.lookup(-CN::SQRT_2, 0)}}});
                    break;
                case BasisStates::right:
                    f = makeVectorNode(static_cast<Qubit>(p), {{{f.p, cn.lookup(CN::SQRT_2, 0)}, {f.p, cn.lookup(0, CN::SQRT_2)}}});
                    break;
                case BasisStates::left:
                    f = makeVectorNode(static_cast<Qubit>(p), {{{f.p, cn.lookup(CN::SQRT_2, 0)}, {f.p, cn.lookup(0, -CN::SQRT_2)}}});
                    break;
            }
        }
        return f;
    }

    Package::mEdge Package::makeMatrixNode(const Qubit var, const std::array<mEdge, NEDGE>& edge, bool cached) {
        mEdge e{mUniqueTable.getNode(), CN::ONE};
        e.p->v = var;
        e.p->e = edge;
        assert(e.p->ref == 0);
        assert(var - 1 == edge[0].p->v || isTerminal(edge[0]));
        assert(var - 1 == edge[1].p->v || isTerminal(edge[1]));
        assert(var - 1 == edge[2].p->v || isTerminal(edge[2]));
        assert(var - 1 == edge[3].p->v || isTerminal(edge[3]));

        // normalize it
        e = normalize(e, cached);
        assert(e.p->v == var || isTerminal(e));

        // look it up in the unique tables
        auto l = mUniqueTable.lookup(e, false);
        assert(l.p->v == var || isTerminal(l));

        // set specific node properties
        if (l.p == e.p)
            checkSpecialMatrices(l.p);

        return l;
    }

    Package::mEdge Package::normalize(const Package::mEdge& e, bool cached) {
        auto argmax = -1;

        auto zero = std::array{CN::equalsZero(e.p->e[0].w),
                               CN::equalsZero(e.p->e[1].w),
                               CN::equalsZero(e.p->e[2].w),
                               CN::equalsZero(e.p->e[3].w)};

        for (auto i = 0U; i < NEDGE; i++) {
            if (zero[i] && e.p->e[i].w != CN::ZERO) {
                cn.releaseCached(e.p->e[i].w);
                e.p->e[i] = mZero;
            }
        }

        fp   max  = 0;
        auto maxc = CN::ONE;
        // determine max amplitude
        for (auto i = 0U; i < NEDGE; ++i) {
            if (zero[i]) continue;
            if (argmax == -1) {
                argmax = static_cast<decltype(argmax)>(i);
                max    = ComplexNumbers::mag2(e.p->e[i].w);
                maxc   = e.p->e[i].w;
            } else {
                auto mag = ComplexNumbers::mag2(e.p->e[i].w);
                if (mag - max > CN::TOLERANCE) {
                    argmax = static_cast<decltype(argmax)>(i);
                    max    = mag;
                    maxc   = e.p->e[i].w;
                }
            }
        }

        // all equal to zero - make sure to release cached numbers approximately zero, but not exactly zero
        if (argmax == -1) {
            if (cached) {
                for (auto const& i: e.p->e) {
                    if (i.w != CN::ZERO) {
                        cn.releaseCached(i.w);
                    }
                }
            } else if (e.p != mTerminalNode) {
                // If it is not a cached variable, I have to put it pack into the chain
                mUniqueTable.returnNode(e.p);
            }
            return mZero;
        }

        auto r = e;
        // divide each entry by max
        for (auto i = 0U; i < NEDGE; ++i) {
            if (static_cast<decltype(argmax)>(i) == argmax) {
                if (cached) {
                    if (r.w == CN::ONE)
                        r.w = maxc;
                    else
                        CN::mul(r.w, r.w, maxc);
                } else {
                    if (r.w == CN::ONE) {
                        r.w = maxc;
                    } else {
                        auto c = cn.getTempCachedComplex();
                        CN::mul(c, r.w, maxc);
                        r.w = cn.lookup(c);
                    }
                }
                r.p->e[i].w = CN::ONE;
            } else {
                if (zero[i]) {
                    if (cached && r.p->e[i].w != ComplexNumbers::ZERO)
                        cn.releaseCached(r.p->e[i].w);
                    r.p->e[i] = mZero;
                    continue;
                }
                if (cached && !zero[i] && r.p->e[i].w != ComplexNumbers::ONE) {
                    cn.releaseCached(r.p->e[i].w);
                }
                if (CN::equalsOne(r.p->e[i].w))
                    r.p->e[i].w = ComplexNumbers::ONE;
                auto c = cn.getTempCachedComplex();
                CN::div(c, r.p->e[i].w, maxc);
                r.p->e[i].w = cn.lookup(c);
            }
        }
        return r;
    }

    Package::mEdge Package::makeGateDD(const std::array<ComplexValue, NEDGE>& mat, Qubit mostSignificantQubit, const std::array<Qubit, MAXN>& line) {
        std::array<mEdge, NEDGE> em{};
        for (auto i = 0U; i < NEDGE; ++i) {
            if (mat[i].r == 0 && mat[i].i == 0) {
                em[i] = mZero;
            } else {
                em[i] = makeMatrixTerminal(cn.lookup(mat[i]));
            }
        }

        //process lines below target
        std::size_t z = 0;
        for (; line[z] < RADIX; z++) {
            for (auto i1 = 0U; i1 < RADIX; i1++) {
                for (auto i2 = 0U; i2 < RADIX; i2++) {
                    auto i = i1 * RADIX + i2;
                    if (line[z] == 0) { // neg. control
                        em[i] = makeMatrixNode(static_cast<Qubit>(z), {em[i], mZero, mZero,
                                                                       (i1 == i2) ? makeIdent(static_cast<Qubit>(z - 1)) : mZero});
                    } else if (line[z] == 1) { // pos. control
                        em[i] = makeMatrixNode(static_cast<Qubit>(z),
                                               {(i1 == i2) ? makeIdent(static_cast<Qubit>(z - 1)) : mZero, mZero, mZero, em[i]});
                    } else { // not connected
                        em[i] = makeMatrixNode(static_cast<Qubit>(z), {em[i], mZero, mZero, em[i]});
                    }
                }
            }
        }

        // target line
        auto e = makeMatrixNode(static_cast<Qubit>(z), em);

        //process lines above target
        for (z++; z <= std::make_unsigned_t<Qubit>(mostSignificantQubit); z++) {
            if (line[z] == 0) { //  neg. control
                e = makeMatrixNode(static_cast<Qubit>(z), {e, mZero, mZero, makeIdent(static_cast<Qubit>(z - 1))});
            } else if (line[z] == 1) { // pos. control
                e = makeMatrixNode(static_cast<Qubit>(z), {makeIdent(static_cast<Qubit>(z - 1)), mZero, mZero, e});
            } else { // not connected
                e = makeMatrixNode(static_cast<Qubit>(z), {e, mZero, mZero, e});
            }
        }
        return e;
    }

    void Package::checkSpecialMatrices(Package::mNode* p) {
        if (p->v == -1)
            return;

        p->ident = false; // assume not identity
        p->symm  = false; // assume symmetric

        // check if matrix is symmetric
        if (!p->e[0].p->symm || !p->e[3].p->symm) return;
        if (transpose(p->e[1]) != p->e[2]) return;
        p->symm = true;

        // check if matrix resembles identity
        if (!(p->e[0].p->ident) || (p->e[1].w) != CN::ZERO || (p->e[2].w) != CN::ZERO || (p->e[0].w) != CN::ONE ||
            (p->e[3].w) != CN::ONE || !(p->e[3].p->ident))
            return;
        p->ident = true;
    }

    std::size_t Package::hash(const Package::vNode* p) {
        std::uintptr_t key = 0;
        key += (reinterpret_cast<std::uintptr_t>(p->e[0].p) +
                reinterpret_cast<std::uintptr_t>(p->e[0].w.r) +
                (reinterpret_cast<std::uintptr_t>(p->e[0].w.i) >> 1)) &
               HASHMASK;
        key &= HASHMASK;
        key += ((reinterpret_cast<std::uintptr_t>(p->e[1].p) >> 2) +
                (reinterpret_cast<std::uintptr_t>(p->e[1].w.r) >> 2) +
                (reinterpret_cast<std::uintptr_t>(p->e[1].w.i) >> 3)) &
               HASHMASK;
        key &= HASHMASK;
        return key;
    }

    std::size_t Package::hash(const Package::mNode* p) {
        std::uintptr_t key = 0;
        key += ((reinterpret_cast<std::uintptr_t>(p->e[0].p)) +
                (reinterpret_cast<std::uintptr_t>(p->e[0].w.r)) +
                (reinterpret_cast<std::uintptr_t>(p->e[0].w.i) >> 1)) &
               HASHMASK;
        key &= HASHMASK;
        key += ((reinterpret_cast<std::uintptr_t>(p->e[1].p) >> 1) +
                (reinterpret_cast<std::uintptr_t>(p->e[1].w.r) >> 1) +
                (reinterpret_cast<std::uintptr_t>(p->e[1].w.i) >> 2)) &
               HASHMASK;
        key &= HASHMASK;
        key += ((reinterpret_cast<std::uintptr_t>(p->e[2].p) >> 2) +
                (reinterpret_cast<std::uintptr_t>(p->e[2].w.r) >> 2) +
                (reinterpret_cast<std::uintptr_t>(p->e[3].w.i) >> 3)) &
               HASHMASK;
        key &= HASHMASK;
        key += ((reinterpret_cast<std::uintptr_t>(p->e[3].p) >> 3) +
                (reinterpret_cast<std::uintptr_t>(p->e[3].w.r) >> 3) +
                (reinterpret_cast<std::uintptr_t>(p->e[3].w.i) >> 4)) &
               HASHMASK;
        key &= HASHMASK;
        return key;
    }

    void Package::garbageCollect(bool force) {
        [[maybe_unused]] auto vCollect = vUniqueTable.garbageCollect(force);
        [[maybe_unused]] auto mCollect = mUniqueTable.garbageCollect(force);
        [[maybe_unused]] auto cCollect = cn.garbageCollect(force);

        // IMPORTANT clear all compute tables
        clearComputeTables();
    }

    void Package::clearComputeTables() {
        vectorAdd.clear();
        matrixAdd.clear();
        matrixTranspose.clear();
        conjugateMatrixTranspose.clear();
        matrixMultiplication.clear();
        matrixVectorMultiplication.clear();
        vectorInnerProduct.clear();
        vectorKronecker.clear();
        matrixKronecker.clear();

        toffoliTable.clear();

        clearIdentityTable();

        operations.clear();
    }

    Package::vEdge Package::add2(const Package::vEdge& x, const Package::vEdge& y) {
        if (x.p == nullptr) return y;
        if (y.p == nullptr) return x;

        if (x.w == CN::ZERO) {
            if (y.w == CN::ZERO) return y;
            auto r = y;
            r.w    = cn.getCachedComplex(CN::val(y.w.r), CN::val(y.w.i));
            return r;
        }
        if (y.w == CN::ZERO) {
            auto r = x;
            r.w    = cn.getCachedComplex(CN::val(x.w.r), CN::val(x.w.i));
            return r;
        }
        if (x.p == y.p) {
            auto r = y;
            r.w    = cn.addCached(x.w, y.w);
            if (CN::equalsZero(r.w)) {
                cn.releaseCached(r.w);
                return vZero;
            }
            return r;
        }

        auto r = vectorAdd.lookup({x.p, x.w}, {y.p, y.w});
        if (r.p != nullptr) {
            if (CN::equalsZero(r.w)) {
                return vZero;
            } else {
                return {r.p, cn.getCachedComplex(r.w)};
            }
        }

        Qubit w;
        if (isTerminal(x)) {
            w = y.p->v;
        } else {
            w = x.p->v;
            if (!isTerminal(y) && y.p->v > w) {
                w = y.p->v;
            }
        }

        std::array<vEdge, RADIX> edge{};
        for (auto i = 0U; i < RADIX; i++) {
            vEdge e1{};
            if (!isTerminal(x) && x.p->v == w) {
                e1 = x.p->e[i];

                if (e1.w != CN::ZERO) {
                    e1.w = cn.mulCached(e1.w, x.w);
                }
            } else {
                e1 = x;
                if (y.p->e[i].p == nullptr) {
                    e1 = {nullptr, CN::ZERO};
                }
            }
            vEdge e2{};
            if (!isTerminal(y) && y.p->v == w) {
                e2 = y.p->e[i];

                if (e2.w != CN::ZERO) {
                    e2.w = cn.mulCached(e2.w, y.w);
                }
            } else {
                e2 = y;
                if (x.p->e[i].p == nullptr) {
                    e2 = {nullptr, CN::ZERO};
                }
            }

            edge[i] = add2(e1, e2);

            if (!isTerminal(x) && x.p->v == w && e1.w != CN::ZERO) {
                cn.releaseCached(e1.w);
            }

            if (!isTerminal(y) && y.p->v == w && e2.w != CN::ZERO) {
                cn.releaseCached(e2.w);
            }
        }

        auto e = makeVectorNode(w, edge, true);
        vectorAdd.insert({x.p, x.w}, {y.p, y.w}, {e.p, e.w});
        return e;
    }

    Package::mEdge Package::add2(const Package::mEdge& x, const Package::mEdge& y) {
        if (x.p == nullptr) return y;
        if (y.p == nullptr) return x;

        if (x.w == CN::ZERO) {
            if (y.w == CN::ZERO) return y;
            auto r = y;
            r.w    = cn.getCachedComplex(CN::val(y.w.r), CN::val(y.w.i));
            return r;
        }
        if (y.w == CN::ZERO) {
            auto r = x;
            r.w    = cn.getCachedComplex(CN::val(x.w.r), CN::val(x.w.i));
            return r;
        }
        if (x.p == y.p) {
            auto r = y;
            r.w    = cn.addCached(x.w, y.w);
            if (CN::equalsZero(r.w)) {
                cn.releaseCached(r.w);
                return mZero;
            }
            return r;
        }

        auto r = matrixAdd.lookup({x.p, x.w}, {y.p, y.w});
        if (r.p != nullptr) {
            if (CN::equalsZero(r.w)) {
                return mZero;
            } else {
                return {r.p, cn.getCachedComplex(r.w)};
            }
        }

        Qubit w;
        if (isTerminal(x)) {
            w = y.p->v;
        } else {
            w = x.p->v;
            if (!isTerminal(y) && y.p->v > w) {
                w = y.p->v;
            }
        }

        std::array<mEdge, NEDGE> edge{};
        for (auto i = 0U; i < NEDGE; i++) {
            mEdge e1{};
            if (!isTerminal(x) && x.p->v == w) {
                e1 = x.p->e[i];

                if (e1.w != CN::ZERO) {
                    e1.w = cn.mulCached(e1.w, x.w);
                }
            } else {
                e1 = x;
                if (y.p->e[i].p == nullptr) {
                    e1 = {nullptr, CN::ZERO};
                }
            }
            mEdge e2{};
            if (!isTerminal(y) && y.p->v == w) {
                e2 = y.p->e[i];

                if (e2.w != CN::ZERO) {
                    e2.w = cn.mulCached(e2.w, y.w);
                }
            } else {
                e2 = y;
                if (x.p->e[i].p == nullptr) {
                    e2 = {nullptr, CN::ZERO};
                }
            }

            edge[i] = add2(e1, e2);

            if (!isTerminal(x) && x.p->v == w && e1.w != CN::ZERO) {
                cn.releaseCached(e1.w);
            }

            if (!isTerminal(y) && y.p->v == w && e2.w != CN::ZERO) {
                cn.releaseCached(e2.w);
            }
        }

        auto e = makeMatrixNode(w, edge, true);
        matrixAdd.insert({x.p, x.w}, {y.p, y.w}, {e.p, e.w});
        return e;
    }

    Package::mEdge Package::transpose(const Package::mEdge& a) {
        if (a.p == nullptr || isTerminal(a) || a.p->symm) {
            return a;
        }

        // check in compute table
        auto r = matrixTranspose.lookup(a, a);
        if (r.p != nullptr) {
            return r;
        }

        std::array<mEdge, NEDGE> e{};
        // transpose sub-matrices and rearrange as required
        for (auto i = 0U; i < RADIX; ++i) {
            for (auto j = 0U; j < RADIX; ++j) {
                e[RADIX * i + j] = transpose(a.p->e[RADIX * j + i]);
            }
        }
        // create new top node
        r = makeMatrixNode(a.p->v, e);
        // adjust top weight
        auto c = cn.getTempCachedComplex();
        CN::mul(c, r.w, a.w);
        r.w = cn.lookup(c);

        // put in compute table
        matrixTranspose.insert(a, a, r);
        return r;
    }

    Package::mEdge Package::conjugateTranspose(const Package::mEdge& a) {
        if (a.p == nullptr)
            return a;
        if (isTerminal(a)) { // terminal case
            auto r = a;
            r.w    = CN::conj(a.w);
            return r;
        }

        // check if in compute table
        auto r = conjugateMatrixTranspose.lookup(a, a);
        if (r.p != nullptr) {
            return r;
        }

        std::array<mEdge, NEDGE> e{};
        // conjugate transpose submatrices and rearrange as required
        for (auto i = 0U; i < RADIX; ++i) {
            for (auto j = 0U; j < RADIX; ++j) {
                e[RADIX * i + j] = conjugateTranspose(a.p->e[RADIX * j + i]);
            }
        }
        // create new top node
        r = makeMatrixNode(a.p->v, e);

        auto c = cn.getTempCachedComplex();
        // adjust top weight including conjugate
        CN::mul(c, r.w, CN::conj(a.w));
        r.w = cn.lookup(c);

        // put it in the compute table
        conjugateMatrixTranspose.insert(a, a, r);
        return r;
    }

    Package::vEdge Package::multiply2(const Package::mEdge& x, const Package::vEdge& y, Qubit var) {
        if (x.p == nullptr) return {nullptr, CN::ZERO};
        if (y.p == nullptr) return y;

        if (x.w == CN::ZERO || y.w == CN::ZERO) {
            return vZero;
        }

        if (var == -1) {
            return makeVectorTerminal(cn.mulCached(x.w, y.w));
        }

        auto xCopy = x;
        xCopy.w    = CN::ONE;
        auto yCopy = y;
        yCopy.w    = CN::ONE;

        auto r = matrixVectorMultiplication.lookup(xCopy, yCopy);
        if (r.p != nullptr) {
            if (CN::equalsZero(r.w)) {
                return vZero;
            } else {
                auto e = vEdge{r.p, cn.getCachedComplex(r.w)};
                CN::mul(e.w, e.w, x.w);
                CN::mul(e.w, e.w, y.w);
                if (CN::equalsZero(e.w)) {
                    cn.releaseCached(e.w);
                    return vZero;
                }
                return e;
            }
        }

        vEdge e{};
        if (x.p->v == var && x.p->v == y.p->v) {
            if (x.p->ident) {
                e = yCopy;
                matrixVectorMultiplication.insert(xCopy, yCopy, {e.p, e.w});
                e.w = cn.mulCached(x.w, y.w);
                if (CN::equalsZero(e.w)) {
                    cn.releaseCached(e.w);
                    return vZero;
                }
                return e;
            }
        }

        std::array<vEdge, RADIX> edge{};
        for (auto i = 0U; i < RADIX; i++) {
            edge[i] = vZero;
            for (auto k = 0U; k < RADIX; k++) {
                mEdge e1{};
                if (!isTerminal(x) && x.p->v == var) {
                    e1 = x.p->e[RADIX * i + k];
                } else {
                    e1 = xCopy;
                }
                vEdge e2{};
                if (!isTerminal(y) && y.p->v == var) {
                    e2 = y.p->e[k];
                } else {
                    e2 = yCopy;
                }

                auto m = multiply2(e1, e2, static_cast<Qubit>(var - 1));

                if (k == 0 || edge[i].w == CN::ZERO) {
                    edge[i] = m;
                } else if (m.w != CN::ZERO) {
                    auto old_e = edge[i];
                    edge[i]    = add2(edge[i], m);
                    cn.releaseCached(old_e.w);
                    cn.releaseCached(m.w);
                }
            }
        }
        e = makeVectorNode(var, edge, true);

        matrixVectorMultiplication.insert(xCopy, yCopy, {e.p, e.w});

        if (e.w != CN::ZERO && (x.w != CN::ONE || y.w != CN::ONE)) {
            if (e.w == CN::ONE) {
                e.w = cn.mulCached(x.w, y.w);
            } else {
                CN::mul(e.w, e.w, x.w);
                CN::mul(e.w, e.w, y.w);
            }
            if (CN::equalsZero(e.w)) {
                cn.releaseCached(e.w);
                return vZero;
            }
        }
        return e;
    }

    Package::mEdge Package::multiply2(const Package::mEdge& x, const Package::mEdge& y, Qubit var) {
        if (x.p == nullptr) return x;
        if (y.p == nullptr) return y;

        if (x.w == CN::ZERO || y.w == CN::ZERO) {
            return mZero;
        }

        if (var == -1) {
            return makeMatrixTerminal(cn.mulCached(x.w, y.w));
        }

        auto xCopy = x;
        xCopy.w    = CN::ONE;
        auto yCopy = y;
        yCopy.w    = CN::ONE;

        auto r = matrixMultiplication.lookup(xCopy, yCopy);
        if (r.p != nullptr) {
            if (CN::equalsZero(r.w)) {
                return mZero;
            } else {
                auto e = mEdge{r.p, cn.getCachedComplex(r.w)};
                CN::mul(e.w, e.w, x.w);
                CN::mul(e.w, e.w, y.w);
                if (CN::equalsZero(e.w)) {
                    cn.releaseCached(e.w);
                    return mZero;
                }
                return e;
            }
        }

        mEdge e{};
        if (x.p->v == var && x.p->v == y.p->v) {
            if (x.p->ident) {
                if (y.p->ident) {
                    e = makeIdent(var);
                } else {
                    e = yCopy;
                }
                matrixMultiplication.insert(xCopy, yCopy, {e.p, e.w});
                e.w = cn.mulCached(x.w, y.w);
                if (CN::equalsZero(e.w)) {
                    cn.releaseCached(e.w);
                    return mZero;
                }
                return e;
            }
            if (y.p->ident) {
                e = xCopy;
                matrixMultiplication.insert(xCopy, yCopy, {e.p, e.w});
                e.w = cn.mulCached(x.w, y.w);

                if (CN::equalsZero(e.w)) {
                    cn.releaseCached(e.w);
                    return mZero;
                }
                return e;
            }
        }

        std::array<mEdge, NEDGE> edge{};
        mEdge                    e1{}, e2{};
        for (auto i = 0U; i < NEDGE; i += RADIX) {
            for (auto j = 0U; j < RADIX; j++) {
                edge[i + j] = mZero;
                for (auto k = 0U; k < RADIX; k++) {
                    if (!isTerminal(x) && x.p->v == var) {
                        e1 = x.p->e[i + k];
                    } else {
                        e1 = xCopy;
                    }
                    if (!isTerminal(y) && y.p->v == var) {
                        e2 = y.p->e[j + RADIX * k];
                    } else {
                        e2 = yCopy;
                    }

                    auto m = multiply2(e1, e2, static_cast<Qubit>(var - 1));

                    if (k == 0 || edge[i + j].w == CN::ZERO) {
                        edge[i + j] = m;
                    } else if (m.w != CN::ZERO) {
                        auto old_e  = edge[i + j];
                        edge[i + j] = add2(edge[i + j], m);
                        cn.releaseCached(old_e.w);
                        cn.releaseCached(m.w);
                    }
                }
            }
        }
        e = makeMatrixNode(var, edge, true);

        matrixMultiplication.insert(xCopy, yCopy, {e.p, e.w});

        if (e.w != CN::ZERO && (x.w != CN::ONE || y.w != CN::ONE)) {
            if (e.w == CN::ONE) {
                e.w = cn.mulCached(x.w, y.w);
            } else {
                CN::mul(e.w, e.w, x.w);
                CN::mul(e.w, e.w, y.w);
            }
            if (CN::equalsZero(e.w)) {
                cn.releaseCached(e.w);
                return mZero;
            }
        }
        return e;
    }

    ComplexValue Package::innerProduct(const Package::vEdge& x, const Package::vEdge& y) {
        if (x.p == nullptr || y.p == nullptr || CN::equalsZero(x.w) || CN::equalsZero(y.w)) { // the 0 case
            return {0, 0};
        }

        [[maybe_unused]] const auto before = cn.cacheCount;

        auto w = x.p->v;
        if (y.p->v > w) {
            w = y.p->v;
        }
        const auto ip = innerProduct(x, y, static_cast<Qubit>(w + 1));

        [[maybe_unused]] const auto after = cn.cacheCount;
        assert(after == before);

        return ip;
    }

    fp Package::fidelity(const Package::vEdge& x, const Package::vEdge& y) {
        const auto fid = innerProduct(x, y);
        return fid.r * fid.r + fid.i * fid.i;
    }

    ComplexValue Package::innerProduct(const Package::vEdge& x, const Package::vEdge& y, Qubit var) {
        if (x.p == nullptr || y.p == nullptr || CN::equalsZero(x.w) || CN::equalsZero(y.w)) { // the 0 case
            return {0.0, 0.0};
        }

        if (var == 0) {
            auto c = cn.getTempCachedComplex();
            CN::mul(c, x.w, y.w);
            return {c.r->val, c.i->val};
        }

        auto xCopy = x;
        xCopy.w    = CN::ONE;
        auto yCopy = y;
        yCopy.w    = CN::ONE;

        auto r = vectorInnerProduct.lookup(xCopy, yCopy);
        if (r.p != nullptr) {
            auto c = cn.getTempCachedComplex(r.w);
            CN::mul(c, c, x.w);
            CN::mul(c, c, y.w);
            return {CN::val(c.r), CN::val(c.i)};
        }

        auto w = static_cast<Qubit>(var - 1);

        ComplexValue sum{0.0, 0.0};
        for (auto i = 0U; i < RADIX; i++) {
            vEdge e1{};
            if (!isTerminal(x) && x.p->v == w) {
                e1 = x.p->e[i];
            } else {
                e1 = xCopy;
            }
            vEdge e2{};
            if (!isTerminal(y) && y.p->v == w) {
                e2   = y.p->e[i];
                e2.w = CN::conj(e2.w);
            } else {
                e2 = yCopy;
            }
            auto cv = innerProduct(e1, e2, w);
            sum.r += cv.r;
            sum.i += cv.i;
        }

        r.p = vTerminalNode;
        r.w = sum;

        vectorInnerProduct.insert(xCopy, yCopy, r);
        auto c = cn.getTempCachedComplex(sum);
        CN::mul(c, c, x.w);
        CN::mul(c, c, y.w);
        return {CN::val(c.r), CN::val(c.i)};
    }

    Package::mEdge Package::extend(const Package::mEdge& e, Qubit h, Qubit l) {
        auto f = (l > 0) ? kronecker(e, makeIdent(static_cast<Qubit>(l - 1))) : e;
        auto g = (h > 0) ? kronecker(makeIdent(static_cast<Qubit>(h - 1)), f) : f;
        return g;
    }

    Package::vEdge Package::kronecker2(const Package::vEdge& x, const Package::vEdge& y) {
        if (CN::equalsZero(x.w))
            return vZero;

        if (isTerminal(x)) {
            auto r = y;
            r.w    = cn.mulCached(x.w, y.w);
            return r;
        }

        auto r = vectorKronecker.lookup(x, y);
        if (r.p != nullptr) {
            if (CN::equalsZero(r.w)) {
                return vZero;
            } else {
                return {r.p, cn.getCachedComplex(r.w)};
            }
        }

        auto e0 = kronecker2(x.p->e[0], y);
        auto e1 = kronecker2(x.p->e[1], y);

        auto e = makeVectorNode(static_cast<Qubit>(y.p->v + x.p->v + 1), {e0, e1}, true);
        CN::mul(e.w, e.w, x.w);
        vectorKronecker.insert(x, y, {e.p, e.w});
        return e;
    }

    Package::mEdge Package::kronecker2(const Package::mEdge& x, const Package::mEdge& y) {
        if (CN::equalsZero(x.w))
            return mZero;

        if (isTerminal(x)) {
            auto r = y;
            r.w    = cn.mulCached(x.w, y.w);
            return r;
        }

        auto r = matrixKronecker.lookup(x, y);
        if (r.p != nullptr) {
            if (CN::equalsZero(r.w)) {
                return mZero;
            } else {
                return {r.p, cn.getCachedComplex(r.w)};
            }
        }

        if (x.p->ident) {
            auto e = makeMatrixNode(static_cast<Qubit>(y.p->v + 1), {y, mZero, mZero, y});
            for (auto i = 0; i < x.p->v; ++i) {
                e = makeMatrixNode(static_cast<Qubit>(e.p->v + 1), {e, mZero, mZero, e});
            }

            e.w = cn.getCachedComplex(CN::val(y.w.r), CN::val(y.w.i));
            matrixKronecker.insert(x, y, {e.p, e.w});
            return e;
        }

        auto e0 = kronecker2(x.p->e[0], y);
        auto e1 = kronecker2(x.p->e[1], y);
        auto e2 = kronecker2(x.p->e[2], y);
        auto e3 = kronecker2(x.p->e[3], y);

        auto e = makeMatrixNode(static_cast<Qubit>(y.p->v + x.p->v + 1), {e0, e1, e2, e3}, true);
        CN::mul(e.w, e.w, x.w);
        matrixKronecker.insert(x, y, {e.p, e.w});
        return e;
    }

    Package::mEdge Package::partialTrace(const Package::mEdge& a, const std::bitset<MAXN>& eliminate) {
        [[maybe_unused]] const auto before = cn.cacheCount;
        const auto                  result = trace(a, eliminate);
        [[maybe_unused]] const auto after  = cn.cacheCount;
        assert(before == after);
        return result;
    }

    ComplexValue Package::trace(const Package::mEdge& a) {
        auto                        eliminate = std::bitset<MAXN>{}.set();
        [[maybe_unused]] const auto before    = cn.cacheCount;
        const auto                  res       = partialTrace(a, eliminate);
        [[maybe_unused]] const auto after     = cn.cacheCount;
        assert(before == after);
        return {ComplexNumbers::val(res.w.r), ComplexNumbers::val(res.w.i)};
    }

    Package::mEdge Package::trace(const Package::mEdge& a, const std::bitset<MAXN>& eliminate, std::size_t alreadyEliminated) {
        auto v = a.p->v;

        if (CN::equalsZero(a.w)) return mZero;

        if (eliminate.none()) return a;

        // Base case
        if (v == -1) {
            if (isTerminal(a)) return a;
            throw std::runtime_error("Expected terminal node in trace.");
        }

        if (eliminate[v]) {
            auto elims = alreadyEliminated + 1;
            auto r     = mZero;

            auto t0 = trace(a.p->e[0], eliminate, elims);
            r       = add2(r, t0);
            auto r1 = r;

            auto t1 = trace(a.p->e[3], eliminate, elims);
            r       = add2(r, t1);
            auto r2 = r;

            if (r.w == CN::ONE) {
                r.w = a.w;
            } else {
                auto c = cn.getTempCachedComplex();
                CN::mul(c, r.w, a.w);
                r.w = cn.lookup(c); // better safe than sorry. this may result in complex values with magnitude > 1 in the complex table
            }

            if (r1.w != CN::ZERO) {
                cn.releaseCached(r1.w);
            }

            if (r2.w != CN::ZERO) {
                cn.releaseCached(r2.w);
            }

            return r;
        } else {
            auto                     adjustedV = static_cast<Qubit>(a.p->v - (eliminate.count() - alreadyEliminated));
            std::array<mEdge, NEDGE> edge{};
            std::transform(a.p->e.cbegin(),
                           a.p->e.cend(),
                           edge.begin(),
                           [&](const mEdge& e) -> mEdge { return trace(e, eliminate, alreadyEliminated); });
            auto r = makeMatrixNode(adjustedV, edge);

            if (r.w == CN::ONE) {
                r.w = a.w;
            } else {
                auto c = cn.getTempCachedComplex();
                CN::mul(c, r.w, a.w);
                r.w = cn.lookup(c);
            }
            return r;
        }
    }

    Package::mEdge Package::makeIdent(Qubit leastSignificantQubit, Qubit mostSignificantQubit) {
        if (mostSignificantQubit < 0)
            return mOne;

        if (leastSignificantQubit == 0 && IdTable[mostSignificantQubit].p != nullptr) {
            return IdTable[mostSignificantQubit];
        }
        if (mostSignificantQubit >= 1 && (IdTable[mostSignificantQubit - 1]).p != nullptr) {
            IdTable[mostSignificantQubit] = makeMatrixNode(mostSignificantQubit,
                                                           {IdTable[mostSignificantQubit - 1],
                                                            mZero,
                                                            mZero,
                                                            IdTable[mostSignificantQubit - 1]});
            return IdTable[mostSignificantQubit];
        }

        auto e = makeMatrixNode(leastSignificantQubit, {mOne, mZero, mZero, mOne});
        for (auto k = std::make_unsigned_t<Qubit>(leastSignificantQubit + 1); k <= mostSignificantQubit; k++) {
            e = makeMatrixNode(k, {e, mZero, mZero, e});
        }
        if (leastSignificantQubit == 0)
            IdTable[mostSignificantQubit] = e;
        return e;
    }

    Package::mEdge Package::reduceAncillae(Package::mEdge& e, const std::bitset<dd::MAXN>& ancillary, bool regular) {
        // return if no more garbage left
        if (!ancillary.any() || e.p == nullptr) return e;
        Qubit lowerbound = 0;
        for (auto i = 0U; i < ancillary.size(); ++i) {
            if (ancillary.test(i)) {
                lowerbound = static_cast<Qubit>(i);
                break;
            }
        }
        if (e.p->v < lowerbound) return e;
        return reduceAncillaeRecursion(e, ancillary, lowerbound, regular);
    }

    Package::vEdge Package::reduceGarbage(Package::vEdge& e, const std::bitset<dd::MAXN>& garbage) {
        // return if no more garbage left
        if (!garbage.any() || e.p == nullptr) return e;
        Qubit lowerbound = 0;
        for (auto i = 0U; i < garbage.size(); ++i) {
            if (garbage.test(i)) {
                lowerbound = static_cast<Qubit>(i);
                break;
            }
        }
        if (e.p->v < lowerbound) return e;
        return reduceGarbageRecursion(e, garbage, lowerbound);
    }

    Package::mEdge Package::reduceGarbage(Package::mEdge& e, const std::bitset<dd::MAXN>& garbage, bool regular) {
        // return if no more garbage left
        if (!garbage.any() || e.p == nullptr) return e;
        Qubit lowerbound = 0;
        for (auto i = 0U; i < garbage.size(); ++i) {
            if (garbage.test(i)) {
                lowerbound = static_cast<Qubit>(i);
                break;
            }
        }
        if (e.p->v < lowerbound) return e;
        return reduceGarbageRecursion(e, garbage, lowerbound, regular);
    }

    Package::mEdge Package::reduceAncillaeRecursion(Package::mEdge& e, const std::bitset<dd::MAXN>& ancillary, Qubit lowerbound, bool regular) {
        if (e.p->v < lowerbound) return e;

        auto f = e;

        std::array<mEdge, NEDGE> edges{};
        std::bitset<NEDGE>       handled{};
        for (auto i = 0U; i < NEDGE; ++i) {
            if (!handled.test(i)) {
                if (isTerminal(e.p->e[i])) {
                    edges[i] = e.p->e[i];
                } else {
                    edges[i] = reduceAncillaeRecursion(f.p->e[i], ancillary, lowerbound, regular);
                    for (auto j = i + 1; j < NEDGE; ++j) {
                        if (e.p->e[i].p == e.p->e[j].p) {
                            edges[j] = edges[i];
                            handled.set(j);
                        }
                    }
                }
                handled.set(i);
            }
        }
        f = makeMatrixNode(f.p->v, edges);

        // something to reduce for this qubit
        if (f.p->v >= 0 && ancillary.test(f.p->v)) {
            if (regular) {
                if (f.p->e[1].w != CN::ZERO || f.p->e[3].w != CN::ZERO) {
                    f = makeMatrixNode(f.p->v, {f.p->e[0], mZero, f.p->e[2], mZero});
                }
            } else {
                if (f.p->e[2].w != CN::ZERO || f.p->e[3].w != CN::ZERO) {
                    f = makeMatrixNode(f.p->v, {f.p->e[0], f.p->e[1], mZero, mZero});
                }
            }
        }

        auto c = cn.mulCached(f.w, e.w);
        f.w    = cn.lookup(c);
        cn.releaseCached(c);

        // increasing the ref count for safety. TODO: find out if this is necessary
        incRef(f);
        return f;
    }

    Package::vEdge Package::reduceGarbageRecursion(Package::vEdge& e, const std::bitset<dd::MAXN>& garbage, Qubit lowerbound) {
        if (e.p->v < lowerbound) return e;

        auto f = e;

        std::array<vEdge, RADIX> edges{};
        std::bitset<RADIX>       handled{};
        for (auto i = 0U; i < RADIX; ++i) {
            if (!handled.test(i)) {
                if (isTerminal(e.p->e[i])) {
                    edges[i] = e.p->e[i];
                } else {
                    edges[i] = reduceGarbageRecursion(f.p->e[i], garbage, lowerbound);
                    for (auto j = i + 1; j < RADIX; ++j) {
                        if (e.p->e[i].p == e.p->e[j].p) {
                            edges[j] = edges[i];
                            handled.set(j);
                        }
                    }
                }
                handled.set(i);
            }
        }
        f = makeVectorNode(f.p->v, edges);

        // something to reduce for this qubit
        if (f.p->v >= 0 && garbage.test(f.p->v)) {
            if (f.p->e[2].w != CN::ZERO) {
                vEdge g{};
                if (f.p->e[0].w == CN::ZERO && f.p->e[2].w != CN::ZERO) {
                    g = f.p->e[2];
                } else if (f.p->e[2].w != CN::ZERO) {
                    g = add(f.p->e[0], f.p->e[2]);
                } else {
                    g = f.p->e[0];
                }
                f = makeVectorNode(e.p->v, {g, vZero});
            }
        }

        auto c = cn.mulCached(f.w, e.w);
        f.w    = cn.lookup(c);
        cn.releaseCached(c);

        // Quick-fix for normalization bug
        if (CN::mag2(f.w) > 1.0)
            f.w = CN::ONE;

        // increasing the ref count for safety. TODO: find out if this is necessary
        incRef(f);
        return f;
    }

    Package::mEdge Package::reduceGarbageRecursion(Package::mEdge& e, const std::bitset<dd::MAXN>& garbage, Qubit lowerbound, bool regular) {
        if (e.p->v < lowerbound) return e;

        auto f = e;

        std::array<mEdge, NEDGE> edges{};
        std::bitset<NEDGE>       handled{};
        for (auto i = 0U; i < NEDGE; ++i) {
            if (!handled.test(i)) {
                if (isTerminal(e.p->e[i])) {
                    edges[i] = e.p->e[i];
                } else {
                    edges[i] = reduceGarbageRecursion(f.p->e[i], garbage, lowerbound, regular);
                    for (auto j = i + 1; j < NEDGE; ++j) {
                        if (e.p->e[i].p == e.p->e[j].p) {
                            edges[j] = edges[i];
                            handled.set(j);
                        }
                    }
                }
                handled.set(i);
            }
        }
        f = makeMatrixNode(f.p->v, edges);

        // something to reduce for this qubit
        if (f.p->v >= 0 && garbage.test(f.p->v)) {
            if (regular) {
                if (f.p->e[2].w != CN::ZERO || f.p->e[3].w != CN::ZERO) {
                    mEdge g{};
                    if (f.p->e[0].w == CN::ZERO && f.p->e[2].w != CN::ZERO) {
                        g = f.p->e[2];
                    } else if (f.p->e[2].w != CN::ZERO) {
                        g = add(f.p->e[0], f.p->e[2]);
                    } else {
                        g = f.p->e[0];
                    }
                    mEdge h{};
                    if (f.p->e[1].w == CN::ZERO && f.p->e[3].w != CN::ZERO) {
                        h = f.p->e[3];
                    } else if (f.p->e[3].w != CN::ZERO) {
                        h = add(f.p->e[1], f.p->e[3]);
                    } else {
                        h = f.p->e[1];
                    }
                    f = makeMatrixNode(e.p->v, {g, h, mZero, mZero});
                }
            } else {
                if (f.p->e[1].w != CN::ZERO || f.p->e[3].w != CN::ZERO) {
                    mEdge g{};
                    if (f.p->e[0].w == CN::ZERO && f.p->e[1].w != CN::ZERO) {
                        g = f.p->e[1];
                    } else if (f.p->e[1].w != CN::ZERO) {
                        g = add(f.p->e[0], f.p->e[1]);
                    } else {
                        g = f.p->e[0];
                    }
                    mEdge h{};
                    if (f.p->e[2].w == CN::ZERO && f.p->e[3].w != CN::ZERO) {
                        h = f.p->e[3];
                    } else if (f.p->e[3].w != CN::ZERO) {
                        h = add(f.p->e[2], f.p->e[3]);
                    } else {
                        h = f.p->e[2];
                    }
                    f = makeMatrixNode(e.p->v, {g, mZero, h, mZero});
                }
            }
        }

        auto c = cn.mulCached(f.w, e.w);
        f.w    = cn.lookup(c);
        cn.releaseCached(c);

        // Quick-fix for normalization bug
        if (CN::mag2(f.w) > 1.0)
            f.w = CN::ONE;

        // increasing the ref count for safety. TODO: find out if this is necessary
        incRef(f);
        return f;
    }

    ComplexValue Package::getValueByPath(const Package::vEdge& e, size_t i) {
        if (dd::Package::isTerminal(e)) {
            return {CN::val(e.w.r), CN::val(e.w.i)};
        }
        return getValueByPath(e, CN::ONE, i);
    }

    ComplexValue Package::getValueByPath(const Package::vEdge& e, const Complex& amp, size_t i) {
        auto c = cn.mulCached(e.w, amp);

        if (isTerminal(e)) {
            cn.releaseCached(c);
            return {CN::val(c.r), CN::val(c.i)};
        }

        bool one = i & (1 << e.p->v);

        ComplexValue r{};
        if (!one && !CN::equalsZero(e.p->e[0].w)) {
            r = getValueByPath(e.p->e[0], c, i);
        } else if (one && !CN::equalsZero(e.p->e[1].w)) {
            r = getValueByPath(e.p->e[1], c, i);
        }
        cn.releaseCached(c);
        return r;
    }
    ComplexValue Package::getValueByPath(const Package::mEdge& e, size_t i, size_t j) {
        if (isTerminal(e)) {
            return {CN::val(e.w.r), CN::val(e.w.i)};
        }
        return getValueByPath(e, CN::ONE, i, j);
    }
    ComplexValue Package::getValueByPath(const Package::mEdge& e, const Complex& amp, size_t i, size_t j) {
        auto c = cn.mulCached(e.w, amp);

        if (isTerminal(e)) {
            cn.releaseCached(c);
            return {CN::val(c.r), CN::val(c.i)};
        }

        bool row = i & (1 << e.p->v);
        bool col = j & (1 << e.p->v);

        ComplexValue r{};
        if (!row && !col && !CN::equalsZero(e.p->e[0].w)) {
            r = getValueByPath(e.p->e[0], c, i, j);
        } else if (!row && col && !CN::equalsZero(e.p->e[1].w)) {
            r = getValueByPath(e.p->e[1], c, i, j);
        } else if (row && !col && !CN::equalsZero(e.p->e[2].w)) {
            r = getValueByPath(e.p->e[2], c, i, j);
        } else if (row && col && !CN::equalsZero(e.p->e[3].w)) {
            r = getValueByPath(e.p->e[3], c, i, j);
        }
        cn.releaseCached(c);
        return r;
    }

    CVec Package::getVector(const Package::vEdge& e) {
        size_t dim = 1 << (e.p->v + 1);
        // allocate resulting vector
        auto vec = CVec(dim, {0.0, 0.0});
        getVector(e, CN::ONE, 0, vec);
        return vec;
    }

    void Package::getVector(const Package::vEdge& e, const Complex& amp, size_t i, CVec& vec) {
        // calculate new accumulated amplitude
        auto c = cn.mulCached(e.w, amp);

        // base case
        if (isTerminal(e)) {
            vec.at(i) = {CN::val(c.r), CN::val(c.i)};
            cn.releaseCached(c);
            return;
        }

        size_t x = i | (1 << e.p->v);

        // recursive case
        if (!CN::equalsZero(e.p->e[0].w))
            getVector(e.p->e[0], c, i, vec);
        if (!CN::equalsZero(e.p->e[1].w))
            getVector(e.p->e[1], c, x, vec);
        cn.releaseCached(c);
    }

    void Package::printVector(const Package::vEdge& e) {
        unsigned long long element = 2u << e.p->v;
        for (unsigned long long i = 0; i < element; i++) {
            auto amplitude = getValueByPath(e, i);
            for (Qubit j = e.p->v; j >= 0; j--) {
                std::cout << ((i >> j) & 1u);
            }
            std::cout << ": " << amplitude << "\n";
        }
        std::cout << std::flush;
    }

    CMat Package::getMatrix(const Package::mEdge& e) {
        size_t dim = 1 << (e.p->v + 1);
        // allocate resulting matrix
        auto mat = CMat(dim, CVec(dim, {0.0, 0.0}));
        getMatrix(e, CN::ONE, 0, 0, mat);
        return mat;
    }

    void Package::getMatrix(const Package::mEdge& e, const Complex& amp, size_t i, size_t j, CMat& mat) {
        // calculate new accumulated amplitude
        auto c = cn.mulCached(e.w, amp);

        // base case
        if (isTerminal(e)) {
            mat.at(i).at(j) = {CN::val(c.r), CN::val(c.i)};
            cn.releaseCached(c);
            return;
        }

        size_t x = i | (1 << e.p->v);
        size_t y = j | (1 << e.p->v);

        // recursive case
        if (!CN::equalsZero(e.p->e[0].w))
            getMatrix(e.p->e[0], c, i, j, mat);
        if (!CN::equalsZero(e.p->e[1].w))
            getMatrix(e.p->e[1], c, i, y, mat);
        if (!CN::equalsZero(e.p->e[2].w))
            getMatrix(e.p->e[2], c, x, j, mat);
        if (!CN::equalsZero(e.p->e[3].w))
            getMatrix(e.p->e[3], c, x, y, mat);
        cn.releaseCached(c);
    }

} // namespace dd
