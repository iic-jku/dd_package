/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DDpackage_H
#define DDpackage_H

#include "ComputeTable.hpp"
#include "Control.hpp"
#include "DDcomplex.hpp"
#include "Definitions.hpp"
#include "Edge.hpp"
#include "OperationTable.hpp"
#include "ToffoliTable.hpp"
#include "UnaryComputeTable.hpp"
#include "UniqueTable.hpp"

#include <algorithm>
#include <array>
#include <bitset>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dd {
    class Package {
        ///
        /// Complex number handling
        ///
    public:
        ComplexNumbers cn{};

        ///
        /// Construction, destruction, information and reset
        ///
    public:
        static constexpr auto maxPossibleQubits = std::numeric_limits<Qubit>::max() + 1;
        static constexpr auto defaultQubits     = 128;
        explicit Package(std::size_t nq = defaultQubits):
            cn(ComplexNumbers()), nqubits(nq) {
            resize(nq);
        };
        ~Package()                      = default;
        Package(const Package& package) = delete;
        Package& operator=(const Package& package) = delete;

        // resize the package instance
        void resize(std::size_t nq) {
            if (nq > maxPossibleQubits) {
                throw std::invalid_argument("Requested too many qubits from package. Qubit datatype only allows up to " +
                                            std::to_string(maxPossibleQubits) + " qubits, while " +
                                            std::to_string(nq) + " were requested. Please recompile the package with a wider Qubit type!");
            }
            nqubits = nq;
            vUniqueTable.resize(nqubits);
            mUniqueTable.resize(nqubits);
            operationTable.resize(nqubits);
            IdTable.resize(nqubits);
        }

        // print information on package and its members
        static void printInformation();

        // print unique and compute table statistics
        void statistics();

        // reset package state
        void reset() {
            clearUniqueTables();
            clearComputeTables();
        }

        // getter for qubits
        [[nodiscard]] auto qubits() const { return nqubits; }

    private:
        std::size_t nqubits;

        ///
        /// Vector nodes, edges and quantum states
        ///
    public:
        struct vNode {
            std::array<Edge<vNode>, RADIX> e{};   // edges out of this node
            RefCount                       ref{}; // reference count
            Qubit                          v{};   // variable index (nonterminal) value (-1 for terminal)

            static vNode            terminalNode;
            constexpr static vNode* terminal{&terminalNode};

            static constexpr bool isTerminal(const vNode* p) { return p == terminal; }
        };
        using vEdge       = Edge<vNode>;
        using vCachedEdge = CachedEdge<vNode>;

        vEdge normalize(const vEdge& e, bool cached);

        // generate |0...0> with n qubits
        vEdge makeZeroState(QubitCount n);
        // generate computational basis state |i> with n qubits
        vEdge makeBasisState(QubitCount n, const std::vector<bool>& state);
        // generate general basis state with n qubits
        vEdge makeBasisState(QubitCount n, const std::vector<BasisStates>& state);

        ///
        /// Matrix nodes, edges and quantum gates
        ///
    public:
        struct mNode {
            std::array<Edge<mNode>, NEDGE> e{};           // edges out of this node
            RefCount                       ref{};         // reference count
            Qubit                          v{};           // variable index (nonterminal) value (-1 for terminal)
            bool                           symm  = false; // node is symmetric
            bool                           ident = false; // node resembles identity

            static mNode            terminalNode;
            constexpr static mNode* terminal{&terminalNode};

            static constexpr bool isTerminal(const mNode* p) { return p == terminal; }
        };
        using mEdge       = Edge<mNode>;
        using mCachedEdge = CachedEdge<mNode>;

        mEdge normalize(const mEdge& e, bool cached);

        // build matrix representation for a single gate on an n-qubit circuit
        mEdge makeGateDD(const std::array<ComplexValue, NEDGE>& mat, QubitCount n, Qubit target) {
            return makeGateDD(mat, n, std::set<Control>{}, target);
        }
        mEdge makeGateDD(const std::array<ComplexValue, NEDGE>& mat, QubitCount n, const Control& control, Qubit target) {
            return makeGateDD(mat, n, std::set{control}, target);
        }
        mEdge makeGateDD(const std::array<ComplexValue, NEDGE>& mat, QubitCount n, Qubit control, Qubit target) {
            return makeGateDD(mat, n, std::set<Control>{{control}}, target);
        }
        mEdge makeGateDD(const std::array<ComplexValue, NEDGE>& mat, QubitCount n, const std::set<Control>& controls, Qubit target);

    private:
        // check whether node represents a symmetric matrix or the identity
        void checkSpecialMatrices(mNode* p);

        ///
        /// Unique tables, Reference counting and garbage collection
        ///
    public:
        // unique tables
        template<class Node>
        [[nodiscard]] UniqueTable<Node>& getUniqueTable();

        template<class Node>
        void incRef(const Edge<Node>& e) {
            getUniqueTable<Node>().incRef(e);
        }
        template<class Node>
        void decRef(const Edge<Node>& e) {
            getUniqueTable<Node>().decRef(e);
        }

        UniqueTable<vNode> vUniqueTable{nqubits};
        UniqueTable<mNode> mUniqueTable{nqubits};

        void garbageCollect(bool force = false);

        void clearUniqueTables() {
            vUniqueTable.clear();
            mUniqueTable.clear();
        }

        // create a normalized DD node and return an edge pointing to it. The node is not recreated if it already exists.
        template<class Node>
        Edge<Node> makeDDNode(Qubit var, const std::array<Edge<Node>, std::tuple_size_v<decltype(Node::e)>>& edges, bool cached = false) {
            auto&      uniqueTable = getUniqueTable<Node>();
            Edge<Node> e{uniqueTable.getNode(), CN::ONE};
            e.p->v = var;
            e.p->e = edges;

            assert(e.p->ref == 0);
            for ([[maybe_unused]] const auto& edge: edges)
                assert(edge.p->v == var - 1 || edge.isTerminal());

            // normalize it
            e = normalize(e, cached);
            assert(e.p->v == var || e.isTerminal());

            // look it up in the unique tables
            auto l = uniqueTable.lookup(e, false);
            assert(l.p->v == var || l.isTerminal());

            // set specific node properties for matrices
            if constexpr (std::tuple_size_v<decltype(Node::e)> == NEDGE) {
                if (l.p == e.p)
                    checkSpecialMatrices(l.p);
            }

            return l;
        }

        ///
        /// Compute table definitions
        ///
    public:
        void clearComputeTables();

        ///
        /// Addition
        ///
    public:
        ComputeTable<vCachedEdge, vCachedEdge, vCachedEdge> vectorAdd{};
        ComputeTable<mCachedEdge, mCachedEdge, mCachedEdge> matrixAdd{};

        template<class Node>
        [[nodiscard]] ComputeTable<CachedEdge<Node>, CachedEdge<Node>, CachedEdge<Node>>& getAddComputeTable();

        template<class Edge>
        Edge add(const Edge& x, const Edge& y) {
            [[maybe_unused]] const auto before = cn.cacheCount;

            auto result = add2(x, y);

            if (result.w != CN::ZERO) {
                cn.releaseCached(result.w);
                result.w = cn.lookup(result.w);
            }

            [[maybe_unused]] const auto after = cn.cacheCount;
            assert(after == before);

            return result;
        }

    private:
        template<class Node>
        Edge<Node> add2(const Edge<Node>& x, const Edge<Node>& y) {
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
                    return Edge<Node>::zero;
                }
                return r;
            }

            auto& computeTable = getAddComputeTable<Node>();
            auto  r            = computeTable.lookup({x.p, x.w}, {y.p, y.w});
            if (r.p != nullptr) {
                if (CN::equalsZero(r.w)) {
                    return Edge<Node>::zero;
                } else {
                    return {r.p, cn.getCachedComplex(r.w)};
                }
            }

            Qubit w;
            if (x.isTerminal()) {
                w = y.p->v;
            } else {
                w = x.p->v;
                if (!y.isTerminal() && y.p->v > w) {
                    w = y.p->v;
                }
            }

            constexpr std::size_t     N = std::tuple_size_v<decltype(x.p->e)>;
            std::array<Edge<Node>, N> edge{};
            for (auto i = 0U; i < N; i++) {
                Edge<Node> e1{};
                if (!x.isTerminal() && x.p->v == w) {
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
                Edge<Node> e2{};
                if (!y.isTerminal() && y.p->v == w) {
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

                if (!x.isTerminal() && x.p->v == w && e1.w != CN::ZERO) {
                    cn.releaseCached(e1.w);
                }

                if (!y.isTerminal() && y.p->v == w && e2.w != CN::ZERO) {
                    cn.releaseCached(e2.w);
                }
            }

            auto e = makeDDNode(w, edge, true);
            computeTable.insert({x.p, x.w}, {y.p, y.w}, {e.p, e.w});
            return e;
        }

        ///
        /// Matrix (conjugate) transpose
        ///
    public:
        UnaryComputeTable<mEdge, mEdge> matrixTranspose{};
        UnaryComputeTable<mEdge, mEdge> conjugateMatrixTranspose{};

        mEdge transpose(const mEdge& a);
        mEdge conjugateTranspose(const mEdge& a);

        ///
        /// Multiplication
        ///
    public:
        ComputeTable<mEdge, vEdge, vCachedEdge> matrixVectorMultiplication{};
        ComputeTable<mEdge, mEdge, mCachedEdge> matrixMultiplication{};
        template<class LeftOperand, class RightOperand>
        RightOperand multiply(const LeftOperand& x, const RightOperand& y) {
            [[maybe_unused]] const auto before = cn.cacheCount;

            Qubit var = -1;
            if (!x.isTerminal()) {
                var = x.p->v;
            }
            if (!y.isTerminal() && (y.p->v) > var) {
                var = y.p->v;
            }

            auto e = multiply2(x, y, var);

            if (e.w != CN::ZERO && e.w != CN::ONE) {
                cn.releaseCached(e.w);
                e.w = cn.lookup(e.w);
            }

            [[maybe_unused]] const auto after = cn.cacheCount;
            assert(before == after);

            return e;
        }

    private:
        vEdge multiply2(const mEdge& x, const vEdge& y, Qubit var);
        mEdge multiply2(const mEdge& x, const mEdge& y, Qubit var);

        ///
        /// Inner product and fidelity
        ///
    public:
        ComputeTable<vEdge, vEdge, vCachedEdge> vectorInnerProduct{};

        ComplexValue innerProduct(const vEdge& x, const vEdge& y);
        fp           fidelity(const vEdge& x, const vEdge& y);

    private:
        ComplexValue innerProduct(const vEdge& x, const vEdge& y, Qubit var);

        ///
        /// Kronecker/tensor product
        ///
    public:
        ComputeTable<vEdge, vEdge, vCachedEdge> vectorKronecker{};
        ComputeTable<mEdge, mEdge, mCachedEdge> matrixKronecker{};

        template<class Node>
        [[nodiscard]] ComputeTable<Edge<Node>, Edge<Node>, CachedEdge<Node>>& getKroneckerComputeTable();

        template<class Edge>
        Edge kronecker(const Edge& x, const Edge& y) {
            auto e = kronecker2(x, y);

            if (e.w != CN::ZERO && e.w != CN::ONE) {
                cn.releaseCached(e.w);
                e.w = cn.lookup(e.w);
            }

            return e;
        }

        // extent the DD pointed to by `e` with `h` identities on top and `l` identities at the bottom
        mEdge extend(const mEdge& e, Qubit h = 0, Qubit l = 0);

    private:
        template<class Node>
        Edge<Node> kronecker2(const Edge<Node>& x, const Edge<Node>& y) {
            if (CN::equalsZero(x.w))
                return Edge<Node>::zero;

            if (x.isTerminal()) {
                auto r = y;
                r.w    = cn.mulCached(x.w, y.w);
                return r;
            }

            auto& computeTable = getKroneckerComputeTable<Node>();
            auto  r            = computeTable.lookup(x, y);
            if (r.p != nullptr) {
                if (CN::equalsZero(r.w)) {
                    return Edge<Node>::zero;
                } else {
                    return {r.p, cn.getCachedComplex(r.w)};
                }
            }

            constexpr std::size_t N = std::tuple_size_v<decltype(x.p->e)>;
            // special case handling for matrices
            if constexpr (N == NEDGE) {
                if (x.p->ident) {
                    auto e = makeDDNode(static_cast<Qubit>(y.p->v + 1), std::array{y, Edge<Node>::zero, Edge<Node>::zero, y});
                    for (auto i = 0; i < x.p->v; ++i) {
                        e = makeDDNode(static_cast<Qubit>(e.p->v + 1), std::array{e, Edge<Node>::zero, Edge<Node>::zero, e});
                    }

                    e.w = cn.getCachedComplex(CN::val(y.w.r), CN::val(y.w.i));
                    computeTable.insert(x, y, {e.p, e.w});
                    return e;
                }
            }

            std::array<Edge<Node>, N> edge{};
            for (auto i = 0U; i < N; ++i) {
                edge[i] = kronecker2(x.p->e[i], y);
            }

            auto e = makeDDNode(static_cast<Qubit>(y.p->v + x.p->v + 1), edge, true);
            CN::mul(e.w, e.w, x.w);
            computeTable.insert(x, y, {e.p, e.w});
            return e;
        }

        ///
        /// (Partial) trace
        ///
    public:
        mEdge        partialTrace(const mEdge& a, const std::vector<bool>& eliminate);
        ComplexValue trace(const mEdge& a);

    private:
        /// TODO: introduce a compute table for the trace?
        mEdge trace(const mEdge& a, const std::vector<bool>& eliminate, std::size_t alreadyEliminated = 0);

        ///
        /// Toffoli gates
        ///
    public:
        ToffoliTable<mEdge> toffoliTable{};

        ///
        /// Identity matrices
        ///
    public:
        // create n-qubit identity DD. makeIdent(n) === makeIdent(0, n-1)
        mEdge makeIdent(QubitCount n) { return makeIdent(0, static_cast<Qubit>(n - 1)); }
        mEdge makeIdent(Qubit leastSignificantQubit, Qubit mostSignificantQubit);

        // identity table access and reset
        [[nodiscard]] const auto& getIdentityTable() const { return IdTable; }

        void clearIdentityTable() {
            for (auto& entry: IdTable) entry.p = nullptr;
        }

    private:
        std::vector<mEdge> IdTable{};

        ///
        /// Operations (related to noise)
        ///
    public:
        OperationTable<mEdge> operationTable{nqubits};

    private:
        ///
        /// Decision diagram size
        ///
    public:
        template<class Edge>
        unsigned int size(const Edge& e) {
            static constexpr unsigned int            NODECOUNT_BUCKETS = 200000;
            static std::unordered_set<decltype(e.p)> visited{NODECOUNT_BUCKETS}; // 2e6
            visited.max_load_factor(10);
            visited.clear();
            return nodeCount(e, visited);
        }

    private:
        template<class Edge>
        unsigned int nodeCount(const Edge& e, std::unordered_set<decltype(e.p)>& v) const {
            v.insert(e.p);
            unsigned int sum = 1;
            if (!e.isTerminal()) {
                for (const auto& edge: e.p->e) {
                    if (edge.p != nullptr && !v.count(edge.p)) {
                        sum += nodeCount(edge, v);
                    }
                }
            }
            return sum;
        }

        ///
        /// Ancillary and garbage reduction
        ///
    public:
        mEdge reduceAncillae(mEdge& e, const std::vector<bool>& ancillary, bool regular = true);

        // Garbage reduction works for reversible circuits --- to be thoroughly tested for quantum circuits
        vEdge reduceGarbage(vEdge& e, const std::vector<bool>& garbage);
        mEdge reduceGarbage(mEdge& e, const std::vector<bool>& garbage, bool regular = true);

    private:
        mEdge reduceAncillaeRecursion(mEdge& e, const std::vector<bool>& ancillary, Qubit lowerbound, bool regular = true);

        vEdge reduceGarbageRecursion(vEdge& e, const std::vector<bool>& garbage, Qubit lowerbound);
        mEdge reduceGarbageRecursion(mEdge& e, const std::vector<bool>& garbage, Qubit lowerbound, bool regular = true);

        ///
        /// Vector and matrix extraction from DDs
        ///
    public:
        /// Get a single element of the vector or matrix represented by the dd with root edge e
        /// \tparam Edge type of edge to use (vector or matrix)
        /// \param e edge to traverse
        /// \param elements string {0, 1, 2, 3}^n describing which outgoing edge should be followed
        ///        (for vectors entries are limitted to 0 and 1)
        ///        If string is longer than required, the additional characters are ignored.
        /// \return the complex amplitude of the specified element
        template<class Edge>
        ComplexValue getValueByPath(const Edge& e, const std::string& elements) {
            if (e.isTerminal()) {
                return {CN::val(e.w.r), CN::val(e.w.i)};
            }

            auto c = cn.getTempCachedComplex(1, 0);
            auto r = e;
            do {
                CN::mul(c, c, r.w);
                std::size_t tmp = elements.at(r.p->v) - '0';
                assert(tmp <= r.p->e.size());
                r = r.p->e.at(tmp);
            } while (!r.isTerminal());
            CN::mul(c, c, r.w);

            return {CN::val(c.r), CN::val(c.i)};
        }
        ComplexValue getValueByPath(const vEdge& e, size_t i);
        ComplexValue getValueByPath(const vEdge& e, const Complex& amp, size_t i);
        ComplexValue getValueByPath(const mEdge& e, size_t i, size_t j);
        ComplexValue getValueByPath(const mEdge& e, const Complex& amp, size_t i, size_t j);

        CVec getVector(const vEdge& e);
        void getVector(const vEdge& e, const Complex& amp, size_t i, CVec& vec);
        void printVector(const vEdge& e);

        CMat getMatrix(const mEdge& e);
        void getMatrix(const mEdge& e, const Complex& amp, size_t i, size_t j, CMat& mat);

        ///
        /// Deserialization
        ///
    public:
        template<class Node, class Edge = Edge<Node>, std::size_t N = std::tuple_size_v<decltype(Node::e)>>
        Edge deserialize(std::istream& is, bool readBinary = false) {
            auto         result = Edge::zero;
            ComplexValue rootweight{};

            std::unordered_map<std::int_fast64_t, Node*> nodes{};
            std::int_fast64_t                            node_index;
            Qubit                                        v;
            std::array<ComplexValue, N>                  edge_weights{};
            std::array<std::int_fast64_t, N>             edge_indices{};
            edge_indices.fill(-2);

            if (readBinary) {
                double version;
                is.read(reinterpret_cast<char*>(&version), sizeof(double));
                if (version != SERIALIZATION_VERSION) {
                    throw std::runtime_error("Wrong Version of serialization file version. version of file: " + std::to_string(version) + "; current version: " + std::to_string(SERIALIZATION_VERSION));
                }

                if (!is.eof()) {
                    rootweight = ComplexValue::readBinary(is);
                }

                while (is.read(reinterpret_cast<char*>(&node_index), sizeof(std::int_fast64_t))) {
                    is.read(reinterpret_cast<char*>(&v), sizeof(Qubit));
                    for (auto i = 0U; i < N; i++) {
                        is.read(reinterpret_cast<char*>(&edge_indices[i]), sizeof(std::int_fast64_t));
                        edge_weights[i] = ComplexValue::readBinary(is);
                    }
                    result = deserializeNode(node_index, v, edge_indices, edge_weights, nodes);
                }
            } else {
                std::string version;
                std::getline(is, version);
                if (std::stod(version) != SERIALIZATION_VERSION) {
                    throw std::runtime_error("Wrong Version of serialization file version. version of file: " + version + "; current version: " + std::to_string(SERIALIZATION_VERSION));
                }

                std::string line;
                std::string complex_real_regex = R"(([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?![ \d\.]*(?:[eE][+-])?\d*[iI]))?)";
                std::string complex_imag_regex = R"(( ?[+-]? ?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?[iI])?)";
                std::string edge_regex         = " \\(((-?\\d+) (" + complex_real_regex + complex_imag_regex + "))?\\)";
                std::regex  complex_weight_regex(complex_real_regex + complex_imag_regex);
                std::string line_construct = "(\\d+) (\\d+)";
                for (auto i = 0U; i < N; ++i) {
                    line_construct += "(?:" + edge_regex + ")";
                }
                line_construct += " *(?:#.*)?";
                std::regex  line_regex(line_construct);
                std::smatch m;

                if (std::getline(is, line)) {
                    if (!std::regex_match(line, m, complex_weight_regex)) {
                        throw std::runtime_error("Regex did not match second line: " + line);
                    }
                    rootweight = ComplexValue::from_string(m.str(1), m.str(2));
                }

                while (std::getline(is, line)) {
                    if (line.empty() || line.size() == 1) continue;

                    if (!std::regex_match(line, m, line_regex)) {
                        throw std::runtime_error("Regex did not match line: " + line);
                    }

                    // match 1: node_idx
                    // match 2: qubit_idx

                    // repeats for every edge
                    // match 3: edge content
                    // match 4: edge_target_idx
                    // match 5: real + imag (without i)
                    // match 6: real
                    // match 7: imag (without i)
                    node_index = std::stoi(m.str(1));
                    v          = static_cast<Qubit>(std::stoi(m.str(2)));

                    for (auto edge_idx = 3U, i = 0U; i < N; i++, edge_idx += 5) {
                        if (m.str(edge_idx).empty()) continue;

                        edge_indices[i] = std::stoi(m.str(edge_idx + 1));
                        edge_weights[i] = ComplexValue::from_string(m.str(edge_idx + 3), m.str(edge_idx + 4));
                    }

                    result = deserializeNode(node_index, v, edge_indices, edge_weights, nodes);
                }
            }

            auto w = cn.getCachedComplex(rootweight.r, rootweight.i);
            CN::mul(w, result.w, w);
            result.w = cn.lookup(w);
            cn.releaseCached(w);

            return result;
        }

        template<class Node, class Edge = Edge<Node>>
        Edge deserialize(const std::string& inputFilename, bool readBinary) {
            auto ifs = std::ifstream(inputFilename);

            if (!ifs.good()) {
                std::cerr << "Cannot open serialized file: " << inputFilename << std::endl;
                return Edge::zero;
            }

            return deserialize<Node>(ifs, readBinary);
        }

    private:
        template<class Node, class Edge = Edge<Node>, std::size_t N = std::tuple_size_v<decltype(Node::e)>>
        Edge deserializeNode(std::int_fast64_t index, Qubit v, std::array<std::int_fast64_t, N>& edge_idx, std::array<ComplexValue, N>& edge_weight, std::unordered_map<std::int_fast64_t, Node*>& nodes) {
            if (index == -1) {
                return Edge::zero;
            }

            std::array<Edge, N> edges{};
            for (auto i = 0U; i < N; ++i) {
                if (edge_idx[i] == -2) {
                    edges[i] = Edge::zero;
                } else {
                    if (edge_idx[i] == -1) {
                        edges[i] = Edge::one;
                    } else {
                        edges[i].p = nodes[edge_idx[i]];
                    }
                    edges[i].w = cn.lookup(edge_weight[i]);
                }
            }

            auto newedge = makeDDNode(v, edges);
            nodes[index] = newedge.p;

            // reset
            edge_idx.fill(-2);

            return newedge;
        }

        ///
        /// Debugging
        ///
    public:
        template<class Node>
        void debugnode(const Node* p) const {
            if (Node::isTerminal(p)) {
                std::clog << "terminal\n";
                return;
            }
            std::clog << "Debug node: " << debugnode_line(p) << "\n";
            for (const auto& edge: p->e) {
                std::clog << "  " << std::hexfloat
                          << std::setw(22) << CN::val(edge.w.r) << " "
                          << std::setw(22) << CN::val(edge.w.i) << std::defaultfloat
                          << "i --> " << debugnode_line(edge.p) << "\n";
            }
            std::clog << std::flush;
        }

        template<class Node>
        std::string debugnode_line(const Node* p) const {
            if (Node::isTerminal(p)) {
                return "terminal";
            }
            std::stringstream sst;
            sst << "0x" << std::hex << reinterpret_cast<std::uintptr_t>(p) << std::dec
                << "[v=" << static_cast<std::int_fast64_t>(p->v)
                << " ref=" << p->ref
                //                << " hash=" << hash(p) TODO: reinclude this
                << "]";
            return sst.str();
        }

        template<class Edge>
        bool isLocallyConsistent(const Edge& e) {
            assert(CN::ONE.r->val == 1 && CN::ONE.i->val == 0);
            assert(CN::ZERO.r->val == 0 && CN::ZERO.i->val == 0);

            const bool result = isLocallyConsistent2(e);

            if (!result) {
                export2Dot(e, "locally_inconsistent.dot", true);
            }

            return result;
        }

        template<class Edge>
        bool isGloballyConsistent(const Edge& e) {
            std::map<ComplexTableEntry*, std::size_t> weight_counter{};
            std::map<decltype(e.p), std::size_t>      node_counter{};
            fillConsistencyCounter(e, weight_counter, node_counter);
            checkConsistencyCounter(e, weight_counter, node_counter);
            return true;
        }

    private:
        template<class Edge>
        bool isLocallyConsistent2(const Edge& e) {
            const auto ptr_r = CN::getAlignedPointer(e.w.r);
            const auto ptr_i = CN::getAlignedPointer(e.w.i);

            if ((ptr_r->ref == 0 || ptr_i->ref == 0) && e.w != CN::ONE && e.w != CN::ZERO) {
                std::clog << "\nLOCAL INCONSISTENCY FOUND\nOffending Number: " << e.w << " (" << ptr_r->ref << ", " << ptr_i->ref << ")\n\n";
                debugnode(e.p);
                return false;
            }

            if (e.isTerminal()) {
                return true;
            }

            if (!e.isTerminal() && e.p->ref == 0) {
                std::clog << "\nLOCAL INCONSISTENCY FOUND: RC==0\n";
                debugnode(e.p);
                return false;
            }

            for (const auto& child: e.p->e) {
                if (child.p->v + 1 != e.p->v && !child.isTerminal()) {
                    std::clog << "\nLOCAL INCONSISTENCY FOUND: Wrong V\n";
                    debugnode(e.p);
                    return false;
                }
                if (!child.isTerminal() && child.p->ref == 0) {
                    std::clog << "\nLOCAL INCONSISTENCY FOUND: RC==0\n";
                    debugnode(e.p);
                    return false;
                }
                if (!isLocallyConsistent2(child)) {
                    return false;
                }
            }
            return true;
        }

        template<class Edge>
        void fillConsistencyCounter(const Edge& edge, std::map<ComplexTableEntry*, std::size_t>& weight_map, std::map<decltype(edge.p), std::size_t>& node_map) {
            weight_map[CN::getAlignedPointer(edge.w.r)]++;
            weight_map[CN::getAlignedPointer(edge.w.i)]++;

            if (edge.isTerminal()) {
                return;
            }
            node_map[edge.p]++;
            for (auto& child: edge.p->e) {
                if (node_map[child.p] == 0) {
                    fillConsistencyCounter(child, weight_map, node_map);
                } else {
                    node_map[child.p]++;
                    weight_map[CN::getAlignedPointer(child.w.r)]++;
                    weight_map[CN::getAlignedPointer(child.w.i)]++;
                }
            }
        }

        template<class Edge>
        void checkConsistencyCounter(const Edge& edge, const std::map<ComplexTableEntry*, std::size_t>& weight_map, const std::map<decltype(edge.p), std::size_t>& node_map) {
            auto* r_ptr = CN::getAlignedPointer(edge.w.r);
            auto* i_ptr = CN::getAlignedPointer(edge.w.i);

            if (weight_map.at(r_ptr) > r_ptr->ref && r_ptr != CN::ONE.r && r_ptr != CN::ZERO.i) {
                std::clog << "\nOffending weight: " << edge.w << "\n";
                std::clog << "Bits: " << std::hexfloat << CN::val(edge.w.r) << " " << CN::val(edge.w.i) << std::defaultfloat << "\n";
                debugnode(edge.p);
                throw std::runtime_error("Ref-Count mismatch for " + std::to_string(r_ptr->val) + "(r): " + std::to_string(weight_map.at(r_ptr)) + " occurences in DD but Ref-Count is only " + std::to_string(r_ptr->ref));
            }

            if (weight_map.at(i_ptr) > i_ptr->ref && i_ptr != CN::ZERO.i && i_ptr != CN::ONE.r) {
                std::clog << edge.w << "\n";
                debugnode(edge.p);
                throw std::runtime_error("Ref-Count mismatch for " + std::to_string(i_ptr->val) + "(i): " + std::to_string(weight_map.at(i_ptr)) + " occurences in DD but Ref-Count is only " + std::to_string(i_ptr->ref));
            }

            if (edge.isTerminal()) {
                return;
            }

            if (node_map.at(edge.p) != edge.p->ref) {
                debugnode(edge.p);
                throw std::runtime_error("Ref-Count mismatch for node: " + std::to_string(node_map.at(edge.p)) + " occurences in DD but Ref-Count is " + std::to_string(edge.p->ref));
            }
            for (auto child: edge.p->e) {
                if (!child.isTerminal() && child.p->v != edge.p->v - 1) {
                    std::clog << "child.p->v == " << child.p->v << "\n";
                    std::clog << " edge.p->v == " << edge.p->v << "\n";
                    debugnode(child.p);
                    debugnode(edge.p);
                    throw std::runtime_error("Variable level ordering seems wrong");
                }
                checkConsistencyCounter(child, weight_map, node_map);
            }
        }
    };

    inline Package::vNode Package::vNode::terminalNode{.e   = {{{nullptr, CN::ZERO}, {nullptr, CN::ZERO}}},
                                                       .ref = 0,
                                                       .v   = -1};

    inline Package::mNode Package::mNode::terminalNode{
            .e     = {{{nullptr, CN::ZERO}, {nullptr, CN::ZERO}, {nullptr, CN::ZERO}, {nullptr, CN::ZERO}}},
            .ref   = 0,
            .v     = -1,
            .symm  = true,
            .ident = true};

    template<>
    [[nodiscard]] inline UniqueTable<Package::vNode>& Package::getUniqueTable() { return vUniqueTable; }

    template<>
    [[nodiscard]] inline UniqueTable<Package::mNode>& Package::getUniqueTable() { return mUniqueTable; }

    template<>
    [[nodiscard]] inline ComputeTable<Package::vCachedEdge, Package::vCachedEdge, Package::vCachedEdge>& Package::getAddComputeTable() { return vectorAdd; }

    template<>
    [[nodiscard]] inline ComputeTable<Package::mCachedEdge, Package::mCachedEdge, Package::mCachedEdge>& Package::getAddComputeTable() { return matrixAdd; }

    template<>
    [[nodiscard]] inline ComputeTable<Package::vEdge, Package::vEdge, Package::vCachedEdge>& Package::getKroneckerComputeTable() { return vectorKronecker; }

    template<>
    [[nodiscard]] inline ComputeTable<Package::mEdge, Package::mEdge, Package::mCachedEdge>& Package::getKroneckerComputeTable() { return matrixKronecker; }
} // namespace dd
#endif