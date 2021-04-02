/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#include "DDpackage.hpp"
#include "GateMatrixDefinitions.hpp"

#include <benchmark/benchmark.h>
#include <memory>

using namespace dd;

static void QubitRange(benchmark::internal::Benchmark* b) {
    b->Unit(benchmark::kMicrosecond)->RangeMultiplier(2)->Range(2, 128);
}

///
/// Test class creation
///
/// At the moment packages are allocated with a fixed maximum number of qubits (=128)
/// In the future, the maximum number of qubits can be set at construction time
/// until then, each creation should take equally long
///

static void BM_DDVectorNodeCreation(benchmark::State& state) {
    for (auto _: state) {
        auto node = Package::vNode{};
        benchmark::DoNotOptimize(node);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_DDVectorNodeCreation)->Unit(benchmark::kNanosecond);

static void BM_DDMatrixNodeCreation(benchmark::State& state) {
    for (auto _: state) {
        auto node = Package::mNode{};
        benchmark::DoNotOptimize(node);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_DDMatrixNodeCreation)->Unit(benchmark::kNanosecond);

static void BM_ComplexNumbersCreation(benchmark::State& state) {
    for (auto _: state) {
        auto cn = ComplexNumbers();
        benchmark::DoNotOptimize(cn);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_ComplexNumbersCreation)->Unit(benchmark::kMicrosecond);

static void BM_PackageCreation(benchmark::State& state) {
    [[maybe_unused]] unsigned short nqubits = state.range(0);
    for (auto _: state) {
        auto dd = std::make_unique<Package>();
        benchmark::DoNotOptimize(dd);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_PackageCreation)->Unit(benchmark::kMillisecond)->RangeMultiplier(2)->Range(2, 128);

///
/// Test creation of identity matrix
///

static void BM_MakeIdentCached(benchmark::State& state) {
    auto msq = static_cast<Qubit>(state.range(0) - 1);
    auto dd  = std::make_unique<Package>();
    for (auto _: state) {
        benchmark::DoNotOptimize(dd->makeIdent(msq));
    }
}
BENCHMARK(BM_MakeIdentCached)->Unit(benchmark::kNanosecond)->RangeMultiplier(2)->Range(2, 128);

static void BM_MakeIdent(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(-1);
        line[0] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Imat, msq, line));
    }
}
BENCHMARK(BM_MakeIdent)->Apply(QubitRange);

///
/// Test makeGateDD
///

static void BM_MakeSingleQubitGateDD_TargetTop(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(-1);
        line[msq] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeSingleQubitGateDD_TargetTop)->Apply(QubitRange);

static void BM_MakeSingleQubitGateDD_TargetMiddle(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(-1);
        line[msq / 2] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeSingleQubitGateDD_TargetMiddle)->Apply(QubitRange);

static void BM_MakeSingleQubitGateDD_TargetBottom(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(-1);
        line[0] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeSingleQubitGateDD_TargetBottom)->Apply(QubitRange);

static void BM_MakeControlledQubitGateDD_ControlBottom_TargetTop(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(-1);
        line[0]   = 1;
        line[msq] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeControlledQubitGateDD_ControlBottom_TargetTop)->Apply(QubitRange);

static void BM_MakeControlledQubitGateDD_ControlBottom_TargetMiddle(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(-1);
        line[0]       = 1;
        line[msq / 2] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeControlledQubitGateDD_ControlBottom_TargetMiddle)->Apply(QubitRange);

static void BM_MakeControlledQubitGateDD_ControlTop_TargetMiddle(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(-1);
        line[msq]     = 1;
        line[msq / 2] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeControlledQubitGateDD_ControlTop_TargetMiddle)->Apply(QubitRange);

static void BM_MakeControlledQubitGateDD_ControlTop_TargetBottom(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(-1);
        line[msq] = 1;
        line[0]   = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeControlledQubitGateDD_ControlTop_TargetBottom)->Apply(QubitRange);

static void BM_MakeFullControlledToffoliDD_TargetTop(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(1);
        line[msq] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeFullControlledToffoliDD_TargetTop)->Apply(QubitRange);

static void BM_MakeFullControlledToffoliDD_TargetMiddle(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(1);
        line[msq / 2] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeFullControlledToffoliDD_TargetMiddle)->Apply(QubitRange);

static void BM_MakeFullControlledToffoliDD_TargetBottom(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    for (auto _: state) {
        line.fill(1);
        line[0] = 2;
        benchmark::DoNotOptimize(dd->makeGateDD(Xmat, msq, line));
    }
}
BENCHMARK(BM_MakeFullControlledToffoliDD_TargetBottom)->Apply(QubitRange);

static void BM_MakeSWAPDD(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};
    line.fill(-1);

    for (auto _: state) {
        line[msq] = 1;
        line[0]   = 2;
        auto sv   = dd->makeGateDD(Xmat, msq, line);
        line[msq] = 2;
        line[0]   = 1;
        sv        = dd->multiply(sv, dd->multiply(dd->makeGateDD(Xmat, msq, line), sv));
        benchmark::DoNotOptimize(sv);
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MakeSWAPDD)->Apply(QubitRange);

///
/// Test multiplication
///

static void BM_MxV_X(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};

    for (auto _: state) {
        auto zero = dd->makeZeroState(msq);
        line.fill(-1);
        line[0]  = 2;
        auto x   = dd->makeGateDD(Xmat, msq, line);
        auto sim = dd->multiply(x, zero);
        benchmark::DoNotOptimize(sim);
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxV_X)->Apply(QubitRange);

static void BM_MxV_H(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};

    for (auto _: state) {
        auto zero = dd->makeZeroState(msq);
        line.fill(-1);
        line[0]  = 2;
        auto h   = dd->makeGateDD(Hmat, msq, line);
        auto sim = dd->multiply(h, zero);
        benchmark::DoNotOptimize(sim);
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxV_H)->Apply(QubitRange);

static void BM_MxV_T(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};

    for (auto _: state) {
        auto zero = dd->makeZeroState(msq);
        line.fill(-1);
        line[0]  = 2;
        auto t   = dd->makeGateDD(Tmat, msq, line);
        auto sim = dd->multiply(t, zero);
        benchmark::DoNotOptimize(sim);
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxV_T)->Apply(QubitRange);

static void BM_MxV_CX_ControlTop_TargetBottom(benchmark::State& state) {
    auto msq         = static_cast<Qubit>(state.range(0) - 1);
    auto dd          = std::make_unique<Package>();
    auto line        = std::array<Qubit, MAXN>{};
    auto basisStates = std::vector<BasisStates>{static_cast<std::size_t>(msq + 1), BasisStates::zero};
    basisStates[msq] = BasisStates::plus;

    for (auto _: state) {
        auto plus = dd->makeBasisState(msq, basisStates);
        line.fill(-1);
        line[0]   = 2;
        line[msq] = 1;
        auto cx   = dd->makeGateDD(Xmat, msq, line);
        auto sim  = dd->multiply(cx, plus);
        benchmark::DoNotOptimize(sim);
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxV_CX_ControlTop_TargetBottom)->Apply(QubitRange);

static void BM_MxV_CX_ControlBottom_TargetTop(benchmark::State& state) {
    auto msq         = static_cast<Qubit>(state.range(0) - 1);
    auto dd          = std::make_unique<Package>();
    auto line        = std::array<Qubit, MAXN>{};
    auto basisStates = std::vector<BasisStates>{static_cast<std::size_t>(msq + 1), BasisStates::zero};
    basisStates[0]   = BasisStates::plus;

    for (auto _: state) {
        auto plus = dd->makeBasisState(msq, basisStates);
        line.fill(-1);
        line[msq] = 2;
        line[0]   = 1;
        auto cx   = dd->makeGateDD(Xmat, msq, line);
        auto sim  = dd->multiply(cx, plus);
        benchmark::DoNotOptimize(sim);
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxV_CX_ControlBottom_TargetTop)->Apply(QubitRange);

static void BM_MxV_HadamardLayer(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};

    for (auto _: state) {
        auto sv = dd->makeZeroState(msq);
        line.fill(-1);
        for (int i = 0; i <= msq; ++i) {
            line[std::max(0, i - 1)] = -1;
            line[i]                  = 2;
            auto h                   = dd->makeGateDD(Hmat, msq, line);
            sv                       = dd->multiply(h, sv);
        }
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxV_HadamardLayer)->Apply(QubitRange);

static void BM_MxV_GHZ(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};

    for (auto _: state) {
        auto sv = dd->makeZeroState(msq);
        line.fill(-1);
        line[msq] = 2;
        auto h    = dd->makeGateDD(Hmat, msq, line);
        sv        = dd->multiply(h, sv);
        line[msq] = 1;
        for (int i = msq - 1; i >= 0; --i) {
            line[std::min(msq - 1, i + 1)] = -1;
            line[i]                        = 2;
            auto cx                        = dd->makeGateDD(Xmat, msq, line);
            sv                             = dd->multiply(cx, sv);
        }
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxV_GHZ)->Apply(QubitRange);

static void BM_MxM_Bell(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};

    for (auto _: state) {
        line.fill(-1);
        line[msq] = 2;
        auto h    = dd->makeGateDD(Hmat, msq, line);
        line[msq] = 1;
        line[0]   = 2;
        auto cx   = dd->makeGateDD(Xmat, msq, line);
        auto bell = dd->multiply(cx, h);
        benchmark::DoNotOptimize(bell);
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxM_Bell)->Apply(QubitRange);

static void BM_MxM_GHZ(benchmark::State& state) {
    auto msq  = static_cast<Qubit>(state.range(0) - 1);
    auto dd   = std::make_unique<Package>();
    auto line = std::array<Qubit, MAXN>{};

    for (auto _: state) {
        line.fill(-1);
        line[msq] = 2;
        auto func = dd->makeGateDD(Hmat, msq, line);
        line[msq] = 1;
        for (int i = msq - 1; i >= 0; --i) {
            line[std::min(msq - 1, i + 1)] = -1;
            line[i]                        = 2;
            auto cx                        = dd->makeGateDD(Xmat, msq, line);
            func                           = dd->multiply(cx, func);
        }
        // clear compute table so the next iteration does not find the result cached
        dd->clearComputeTables();
    }
}
BENCHMARK(BM_MxM_GHZ)->Apply(QubitRange);
