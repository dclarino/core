/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#pragma once

#include "dd/CachedEdge.hpp"
#include "dd/Complex.hpp"
#include "dd/ComplexNumbers.hpp"
#include "dd/ComplexValue.hpp"
#include "dd/ComputeTable.hpp"
#include "dd/DDDefinitions.hpp"
#include "dd/DDpackageConfig.hpp"
#include "dd/DensityNoiseTable.hpp"
#include "dd/Edge.hpp"
#include "dd/MemoryManager.hpp"
#include "dd/Node.hpp"
#include "dd/Package_fwd.hpp" // IWYU pragma: export
#include "dd/RealNumber.hpp"
#include "dd/RealNumberUniqueTable.hpp"
#include "dd/StochasticNoiseOperationTable.hpp"
#include "dd/UnaryComputeTable.hpp"
#include "dd/UniqueTable.hpp"
#include "ir/Definitions.hpp"
#include "ir/Permutation.hpp"
#include "ir/operations/Control.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <regex>
#include <stack>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace dd {

class Package {
public:
  static constexpr std::size_t MAX_POSSIBLE_QUBITS =
      static_cast<std::size_t>(std::numeric_limits<Qubit>::max()) + 1U;
  static constexpr std::size_t DEFAULT_QUBITS = 32U;
  explicit Package(std::size_t nq = DEFAULT_QUBITS,
                   const DDPackageConfig& config = DDPackageConfig{});
  ~Package() = default;
  Package(const Package& package) = delete;

  Package& operator=(const Package& package) = delete;

  void resize(std::size_t nq);
  void reset();
  [[nodiscard]] auto qubits() const { return nqubits; }

private:
  std::size_t nqubits;
  DDPackageConfig config_;

public:
  MemoryManager vMemoryManager{
      MemoryManager::create<vNode>(config_.utVecInitialAllocationSize)};
  MemoryManager mMemoryManager{
      MemoryManager::create<mNode>(config_.utMatInitialAllocationSize)};
  MemoryManager dMemoryManager{
      MemoryManager::create<dNode>(config_.utDmInitialAllocationSize)};
  MemoryManager cMemoryManager{MemoryManager::create<RealNumber>()};

  template <class T> [[nodiscard]] auto& getMemoryManager() {
    if constexpr (std::is_same_v<T, vNode>) {
      return vMemoryManager;
    } else if constexpr (std::is_same_v<T, mNode>) {
      return mMemoryManager;
    } else if constexpr (std::is_same_v<T, dNode>) {
      return dMemoryManager;
    } else if constexpr (std::is_same_v<T, RealNumber>) {
      return cMemoryManager;
    }
  }

  void resetMemoryManagers(bool resizeToTotal = false);

  UniqueTable vUniqueTable{vMemoryManager, {0U, config_.utVecNumBucket}};
  UniqueTable mUniqueTable{mMemoryManager, {0U, config_.utMatNumBucket}};
  UniqueTable dUniqueTable{dMemoryManager, {0U, config_.utDmNumBucket}};
  RealNumberUniqueTable cUniqueTable{cMemoryManager};
  ComplexNumbers cn{cUniqueTable};

  template <class T> [[nodiscard]] auto& getUniqueTable() {
    if constexpr (std::is_same_v<T, vNode>) {
      return vUniqueTable;
    } else if constexpr (std::is_same_v<T, mNode>) {
      return mUniqueTable;
    } else if constexpr (std::is_same_v<T, dNode>) {
      return dUniqueTable;
    } else if constexpr (std::is_same_v<T, RealNumber>) {
      return cUniqueTable;
    }
  }

  void clearUniqueTables();

  template <class Node> void incRef(const Edge<Node>& e) noexcept {
    if (Edge<Node>::trackingRequired(e)) {
      roots.addToRoots(e);
    }
  }

  template <class Node> void decRef(const Edge<Node>& e) {
    if (Edge<Node>::trackingRequired(e)) {
      roots.removeFromRoots(e);
    }
  }

  template <class Node> [[nodiscard]] auto& getRootSet() noexcept {
    return roots.getRoots<Node>();
  }

private:
  struct RootSetManager {
  public:
    template <class Node>
    using RootSet = std::unordered_map<Edge<Node>, std::size_t>;

    template <class Node> void addToRoots(const Edge<Node>& e) noexcept {
      getRoots<Node>()[e]++;
    }
    template <class Node> void removeFromRoots(const Edge<Node>& e) {
      auto& set = getRoots<Node>();
      auto it = set.find(e);
      if (it == set.end()) {
        throw std::invalid_argument("Edge is not part of the root set.");
      }
      if (--it->second == 0U) {
        set.erase(it);
      }
    }
    template <class Result, typename Fn> Result execute(Fn& op) noexcept {
      mark();
      Result res = op();
      unmark();
      return res;
    }
    void reset() {
      vRoots.clear();
      mRoots.clear();
      dRoots.clear();
    }

  private:
    template <class Node> static void mark(const RootSet<Node>& roots) {
      for (auto& [edge, _] : roots) {
        edge.mark();
      }
    }
    template <class Node> static void unmark(const RootSet<Node>& roots) {
      for (auto& [edge, _] : roots) {
        edge.unmark();
      }
    }
    void mark() noexcept {
      RootSetManager::mark(vRoots);
      RootSetManager::mark(mRoots);
      RootSetManager::mark(dRoots);
    }
    void unmark() noexcept {
      RootSetManager::unmark(vRoots);
      RootSetManager::unmark(mRoots);
      RootSetManager::unmark(dRoots);
    }
    template <class Node,
              std::enable_if_t<std::is_same_v<Node, vNode>, bool> = true>
    auto& getRoots() noexcept {
      return vRoots;
    }
    template <class Node,
              std::enable_if_t<std::is_same_v<Node, mNode>, bool> = true>
    auto& getRoots() noexcept {
      return mRoots;
    }
    template <class Node,
              std::enable_if_t<std::is_same_v<Node, dNode>, bool> = true>
    auto& getRoots() noexcept {
      return dRoots;
    }
    RootSet<vNode> vRoots;
    RootSet<mNode> mRoots;
    RootSet<dNode> dRoots;
    template <class Node> friend auto& Package::getRootSet() noexcept;
  };

  RootSetManager roots;

public:
  bool garbageCollect(bool force = false);

  struct ActiveCounts {
    std::size_t vector = 0U;
    std::size_t matrix = 0U;
    std::size_t density = 0U;
    std::size_t reals = 0U;
  };
  [[nodiscard]] ActiveCounts computeActiveCounts();

  dEdge makeZeroDensityOperator(std::size_t n);

  mEdge makeGateDD(const GateMatrix& mat, qc::Qubit target);
  mEdge makeGateDD(const GateMatrix& mat, const qc::Control& control,
                   qc::Qubit target);
  mEdge makeGateDD(const GateMatrix& mat, const qc::Controls& controls,
                   qc::Qubit target);

  mEdge makeTwoQubitGateDD(const TwoQubitGateMatrix& mat, qc::Qubit target0,
                           qc::Qubit target1);
  mEdge makeTwoQubitGateDD(const TwoQubitGateMatrix& mat,
                           const qc::Control& control, qc::Qubit target0,
                           qc::Qubit target1);
  mEdge makeTwoQubitGateDD(const TwoQubitGateMatrix& mat,
                           const qc::Controls& controls, qc::Qubit target0,
                           qc::Qubit target1);

  mEdge makeDDFromMatrix(const CMat& matrix);

private:
  mCachedEdge makeDDFromMatrix(const CMat& matrix, Qubit level,
                               std::size_t rowStart, std::size_t rowEnd,
                               std::size_t colStart, std::size_t colEnd);

public:
  template <class Node, template <class> class EdgeType>
  EdgeType<Node>
  makeDDNode(const Qubit var,
             const std::array<EdgeType<Node>,
                              std::tuple_size_v<decltype(Node::e)>>& edges,
             [[maybe_unused]] const bool generateDensityMatrix = false);

  template <class Node>
  Edge<Node> deleteEdge(const Edge<Node>& e, const Qubit v,
                        const std::size_t edgeIdx);
  template <class Node>
  Edge<Node> deleteEdge(const Edge<Node>& e, const Qubit v,
                        const std::size_t edgeIdx,
                        std::unordered_map<Node*, Edge<Node>>& nodes);

  void clearComputeTables();

  std::string measureAll(vEdge& rootEdge, bool collapse, std::mt19937_64& mt,
                         fp epsilon = 0.001);

private:
  static fp assignProbabilities(const vEdge& edge,
                                std::unordered_map<const vNode*, fp>& probs);

public:
  static std::pair<fp, fp>
  determineMeasurementProbabilities(const vEdge& rootEdge, Qubit index);

  char measureOneCollapsing(vEdge& rootEdge, Qubit index, std::mt19937_64& mt,
                            fp epsilon = 0.001);

  char measureOneCollapsing(dEdge& e, Qubit index, std::mt19937_64& mt);

  void performCollapsingMeasurement(vEdge& rootEdge, Qubit index,
                                    fp probability, bool measureZero);

  ComputeTable<vCachedEdge, vCachedEdge, vCachedEdge> vectorAdd{
      config_.ctVecAddNumBucket};
  ComputeTable<mCachedEdge, mCachedEdge, mCachedEdge> matrixAdd{
      config_.ctMatAddNumBucket};
  ComputeTable<dCachedEdge, dCachedEdge, dCachedEdge> densityAdd{
      config_.ctDmAddNumBucket};

  template <class Node> [[nodiscard]] auto& getAddComputeTable() {
    if constexpr (std::is_same_v<Node, vNode>) {
      return vectorAdd;
    } else if constexpr (std::is_same_v<Node, mNode>) {
      return matrixAdd;
    } else if constexpr (std::is_same_v<Node, dNode>) {
      return densityAdd;
    }
  }

  ComputeTable<vCachedEdge, vCachedEdge, vCachedEdge> vectorAddMagnitudes{
      config_.ctVecAddMagNumBucket};
  ComputeTable<mCachedEdge, mCachedEdge, mCachedEdge> matrixAddMagnitudes{
      config_.ctMatAddMagNumBucket};

  template <class Node> [[nodiscard]] auto& getAddMagnitudesComputeTable() {
    if constexpr (std::is_same_v<Node, vNode>) {
      return vectorAddMagnitudes;
    } else if constexpr (std::is_same_v<Node, mNode>) {
      return matrixAddMagnitudes;
    }
  }

  template <class Node>
  Edge<Node> add(const Edge<Node>& x, const Edge<Node>& y);

  template <class Node>
  CachedEdge<Node> add2(const CachedEdge<Node>& x, const CachedEdge<Node>& y,
                        const Qubit var);

  template <class Node>
  CachedEdge<Node> addMagnitudes(const CachedEdge<Node>& x,
                                 const CachedEdge<Node>& y, const Qubit var);

  UnaryComputeTable<vNode*, vCachedEdge> conjugateVector{
      config_.ctVecConjNumBucket};

  vEdge conjugate(const vEdge& a);
  vCachedEdge conjugateRec(const vEdge& a);

  UnaryComputeTable<mNode*, mCachedEdge> conjugateMatrixTranspose{
      config_.ctMatConjTransNumBucket};

  mEdge conjugateTranspose(const mEdge& a);
  mCachedEdge conjugateTransposeRec(const mEdge& a);

  ComputeTable<mNode*, vNode*, vCachedEdge> matrixVectorMultiplication{
      config_.ctMatVecMultNumBucket};
  ComputeTable<mNode*, mNode*, mCachedEdge> matrixMatrixMultiplication{
      config_.ctMatMatMultNumBucket};
  ComputeTable<dNode*, dNode*, dCachedEdge> densityDensityMultiplication{
      config_.ctDmDmMultNumBucket};

  template <class RightOperandNode>
  [[nodiscard]] auto& getMultiplicationComputeTable() {
    if constexpr (std::is_same_v<RightOperandNode, vNode>) {
      return matrixVectorMultiplication;
    } else if constexpr (std::is_same_v<RightOperandNode, mNode>) {
      return matrixMatrixMultiplication;
    } else if constexpr (std::is_same_v<RightOperandNode, dNode>) {
      return densityDensityMultiplication;
    }
  }

  VectorDD applyOperation(const MatrixDD& operation, const VectorDD& e);
  MatrixDD applyOperation(const MatrixDD& operation, const MatrixDD& e,
                          bool applyFromLeft = true);
  dEdge applyOperationToDensity(dEdge& e, const mEdge& operation);

  template <class LeftOperandNode, class RightOperandNode>
  Edge<RightOperandNode>
  multiply(const Edge<LeftOperandNode>& x, const Edge<RightOperandNode>& y,
           [[maybe_unused]] const bool generateDensityMatrix = false);

private:
  template <class LeftOperandNode, class RightOperandNode>
  CachedEdge<RightOperandNode>
  multiply2(const Edge<LeftOperandNode>& x, const Edge<RightOperandNode>& y,
            const Qubit var,
            [[maybe_unused]] const bool generateDensityMatrix = false);

public:
  ComputeTable<vNode*, vNode*, vCachedEdge> vectorInnerProduct{
      config_.ctVecInnerProdNumBucket};

  ComplexValue innerProduct(const vEdge& x, const vEdge& y);
  fp fidelity(const vEdge& x, const vEdge& y);
  static fp
  fidelityOfMeasurementOutcomes(const vEdge& e, const SparsePVec& probs,
                                const qc::Permutation& permutation = {});

private:
  ComplexValue innerProduct(const vEdge& x, const vEdge& y, Qubit var);

  static fp fidelityOfMeasurementOutcomesRecursive(
      const vEdge& e, const SparsePVec& probs, std::size_t i,
      const qc::Permutation& permutation, std::size_t nQubits);

public:
  fp expectationValue(const mEdge& x, const vEdge& y);

  ComputeTable<vNode*, vNode*, vCachedEdge> vectorKronecker{
      config_.ctVecKronNumBucket};
  ComputeTable<mNode*, mNode*, mCachedEdge> matrixKronecker{
      config_.ctMatKronNumBucket};

  template <class Node> [[nodiscard]] auto& getKroneckerComputeTable() {
    if constexpr (std::is_same_v<Node, vNode>) {
      return vectorKronecker;
    } else {
      return matrixKronecker;
    }
  }

  template <class Node>
  Edge<Node> kronecker(const Edge<Node>& x, const Edge<Node>& y,
                       const std::size_t yNumQubits, const bool incIdx = true);

private:
  template <class Node>
  CachedEdge<Node> kronecker2(const Edge<Node>& x, const Edge<Node>& y,
                              const std::size_t yNumQubits,
                              const bool incIdx = true);

public:
  UnaryComputeTable<dNode*, dCachedEdge> densityTrace{
      config_.ctDmTraceNumBucket};
  UnaryComputeTable<mNode*, mCachedEdge> matrixTrace{
      config_.ctMatTraceNumBucket};

  template <class Node> [[nodiscard]] auto& getTraceComputeTable() {
    if constexpr (std::is_same_v<Node, mNode>) {
      return matrixTrace;
    } else {
      return densityTrace;
    }
  }

  mEdge partialTrace(const mEdge& a, const std::vector<bool>& eliminate);

  template <class Node>
  ComplexValue trace(const Edge<Node>& a, const std::size_t numQubits);

  [[nodiscard]] bool isCloseToIdentity(const mEdge& m, fp tol = 1e-10,
                                       const std::vector<bool>& garbage = {},
                                       bool checkCloseToOne = true) const;

private:
  template <class Node>
  CachedEdge<Node> trace(const Edge<Node>& a,
                         const std::vector<bool>& eliminate, std::size_t level,
                         std::size_t alreadyEliminated = 0);

  static bool isCloseToIdentityRecursive(
      const mEdge& m, std::unordered_set<decltype(m.p)>& visited, fp tol,
      const std::vector<bool>& garbage, bool checkCloseToOne);

public:
  static mEdge makeIdent();
  mEdge createInitialMatrix(const std::vector<bool>& ancillary);

  StochasticNoiseOperationTable<mEdge> stochasticNoiseOperationCache{
      nqubits, config_.stochasticCacheOps};
  DensityNoiseTable<dEdge, dEdge> densityNoise{config_.ctDmNoiseNumBucket};

  mEdge reduceAncillae(mEdge e, const std::vector<bool>& ancillary,
                       bool regular = true);
  vEdge reduceGarbage(vEdge& e, const std::vector<bool>& garbage,
                      bool normalizeWeights = false);
  mEdge reduceGarbage(const mEdge& e, const std::vector<bool>& garbage,
                      bool regular = true, bool normalizeWeights = false);

private:
  mCachedEdge reduceAncillaeRecursion(mNode* p,
                                      const std::vector<bool>& ancillary,
                                      Qubit lowerbound, bool regular = true);

  vCachedEdge reduceGarbageRecursion(vNode* p, const std::vector<bool>& garbage,
                                     Qubit lowerbound,
                                     bool normalizeWeights = false);
  mCachedEdge reduceGarbageRecursion(mNode* p, const std::vector<bool>& garbage,
                                     Qubit lowerbound, bool regular = true,
                                     bool normalizeWeights = false);

public:
  template <class Node> Edge<Node> transfer(Edge<Node>& original);

  template <class Node, class Edge = Edge<Node>,
            std::size_t N = std::tuple_size_v<decltype(Node::e)>>
  Edge deserialize(std::istream& is, const bool readBinary = false);

  template <class Node, class Edge = Edge<Node>>
  Edge deserialize(const std::string& inputFilename, const bool readBinary);

private:
  template <class Node, std::size_t N = std::tuple_size_v<decltype(Node::e)>>
  CachedEdge<Node>
  deserializeNode(const std::int64_t index, const Qubit v,
                  std::array<std::int64_t, N>& edgeIdx,
                  const std::array<ComplexValue, N>& edgeWeight,
                  std::unordered_map<std::int64_t, Node*>& nodes);
};

} // namespace dd
