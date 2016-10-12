// Copyright (C) 2010 Michael E. Houle
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef INDEXING_RCT_PUBLIC_RANK_COVER_TREE_H_
#define INDEXING_RCT_PUBLIC_RANK_COVER_TREE_H_

#include <random>
#include <vector>
  
namespace indexing {
namespace rct {

// TODO: Bring back load/save support.
// TODO: Extract caching logic and make it compatible with any DistanceType.
template <typename DomainType, typename DistanceType,
          DistanceType (*DistanceFn)(const DomainType&, const DomainType&)>
class RankCoverTree {
 public:
  // TODO: Migrate to factory pattern.
  RankCoverTree() {}
  virtual ~RankCoverTree() {}

  // Builds an RCT from a collection of data items.
  //
  // TODO: Extract into factory class.
  // TODO: Specify desired height rather than sample probability.
  bool Build(const std::vector<const DomainType*>& data,
             const double build_coverage, const size_t num_parents,
             const double sample_probability);

  // Accept a new query item.
  //
  // TODO: Fold this in with the actual query operation.
  void SetNewQuery(const DomainType* query);

  // Prints statistics related to the RCT construction.
  void PrintStats() const;


  //! Loads a previously-saved RCT from a file.
  // int build(const char* fileName, DomainType** inputData, const int&
  // numItems);

  //! Retrieve the fraction of edges in the RCT which are well-formed.
  // double getFractionOfWellformedEdges();

  //! Create an RCT with a specific seed.
  // RCT(const unsigned long& seed = 3141569UL);

  //! Set the sample rate.
  // void setSampleRate(const float& sampleRate);

  //! Perform an approximate nearest-neighbor query.
  // int findNear(DomainType* query, const int& howMany = 1,
  //              const float& scaleFactor = 1.0f, const int& sampleLevel = 0);

  //! Perform an exact nearest-neighbor query.
  // int findNearest(DomainType* query, const int& howMany = 1,
  //                 const int& sampleLevel = 0);

  //! Retrieve the average node degree.
  // float getAvgDegree() const;

  //! Retrieve the build scale factor.
  // float getBuildScaleFactor() const;

  //! Retrieve the coverage parameter.
  // float getCoverageParameter() const;

  //! Retrieve the data items.
  // DomainType** getData();

  //! Retrieve a mapping from external to internal item indices.
  // int getExternToInternMapping(int* result, int capacity) const;

  //! Retrieve a mapping from internal to external item indices.
  // int getInternToExternMapping(int* result, int capacity) const;

  //! Retrieve the level set sizes.
  // int getLevelSetSizes(int* result, int capacity) const;

  //! Retrieve the maximum node degree.
  // int getMaxDegree() const;

  //! Retrieve the height of the items in the RCT.
  // int getMaxLevelAssignment(int* result, int capacity) const;

  //! Retrieve the maximum number of parents allowed.
  // int getMaxParents() const;

  //! Retrieve the number of items in the RCT.
  // int getNumItems() const;

  //! Retrieve the number of levels in the RCT.
  // int getNumLevels() const;

  //! Retrieve the number of nodes in the RCT.
  // int getNumNodes() const;

  //! Computes the recall accuracy.
  // float getResultAcc(float* exactDistList, int howMany) const;

  //! Retrieve the query-to-neighbor distances.
  // int getResultDists(float* result, int capacity) const;

  //! Retrieve the number of distance comparisons performed.
  // unsigned long getResultDistComps() const;

  //! Retrieve the indices of the query result items.
  // int getResultIndices(int* result, int capacity) const;

  //! Retrieve the number of results found.
  // int getResultNumFound() const;

  //! Retrieve the sample size used in the query.
  // int getResultSampleSize() const;

  //! Retrieve the random number generator seed.
  // unsigned long getRNGSeed() const;

  //! Save the RCT to a file.
  // int saveToFile(const char* fileName) const;

  //! Set the coverage parameter.
  // bool setCoverageParameter(const float& coverageParameter);

  //! Set the verbosity of the RCT.
  // void setVerbosity(const int& verbosity);

 private:
  // Builds an RCT on a set of data items.
  //
  // TODO: Move this into a factory.
  void BuildIncrementally();

  // Sets up random leveling, reserves RCT storage, and sets up index
  // parameters. As a result of this operation, the RCT size, number of levels,
  // etc, are set.
  void SetupLevels(const size_t num_items, const size_t num_parents,
                   const double sample_probability);

  // Allocate storage for the RCT and its data.
  //
  // TODO: Migrate this to a factory class.
  void AllocateStorage();

  // Computes query-distance of the provided object.
  DistanceType ComputeDistFromQuery(const int index);

  // Performs approximate nearest neighbor search with respect to a given level
  // of the RCT.
  int FindNear(const int how_many, const int sample_level, const double coverage);

  // Partial quicksort.
  //
  // TODO: Deprecate this in favor of std::partial_sort with mutable
  // zip-iterators.
  int PartialQuickSort(int howMany, float* distList, int* indexList,
                       int rangeFirst, int rangeLast);

  // //! Resets the query.
  // void resetQuery();

  // //! Build an RCT on data items.
  // void doBuild();

  // //! Performs an exact range query.
  // int doFindAllInRange(float limit, int sampleLevel);

  // //! Performs an approximate range query.
  // int doFindMostInRange(float limit, int sampleLevel, float scaleFactor);

  // //! Performs an exact nearest-neighbor query.
  // int doFindNearest(int howMany, int sampleLevel);


  // //! Quicksort.
  // void quickSort(float* distList, int* indexList, int rangeFirst,
  //                int rangeLast);

  // Data elements stored in the RCT.
  //
  // The elements are not owned by the index structure and it is the
  // responsibility of the user to ensure the data elements remain alive until
  // the RCT is destroyed.
  std::vector<const DomainType*> data_;

  // Random number generator used by the RCT.
  std::mt19937_64 rand_;

  // Sample rate.
  //
  // The probability of an item in level set i also occurring in level set i + 1
  // is the inverse of the sample rate. Therefore, a sample rate of 2.0 is
  // expected to yield a random leveling where each successive random level
  // contains roughly half as many elements as the level below it. A using the
  // j-th root of the number of items provides a random leveling with expected
  // height j.
  double sample_rate_;

  // Size of the level sets.
  //
  // Stores the number of element in the level sets L_0,L_1,... depending on how
  // the items are distributed across levels by the random leveling process.
  // 
  // TODO: Use managed container.
  int* levelSetSizeList;

  // Maximum number of parents per node.
  //
  // This value fixes the maximum number of parents per node. To obtain a rank
  // cover tree this value must be set to 1. Using larger values than 1 results
  // in a non-tree structure. Larger values can improve search accuracies at the
  // cost of query and construction time complexity.
  size_t max_parents_;

  // Maximum node degree.
  //
  // The maximum number of children per node. These include the virtual
  // super-node constructed on top of the highest level in the tree (in case
  // that
  // level contains more than one item).
  int max_degree_;

  // Average node degree.
  //
  // The average number of children per node. These are determined as a result
  // of the build step.
  float average_degree_;

  // Mapping between internal and external indices.
  //
  // Stores the mapping from internal item indices to external (input) indices.
  // This mapping represents the internal permutation of data items done by the
  // RCT during the random leveling. The original data array is not changed.
  std::vector<size_t> intern_to_extern_mapping_;

  // Height of the RCT.
  //
  // Number of sample levels in the rank cover tree (other than the artificial
  // root's). The bottom rank cover tree level has index 0.
  size_t levels;

  // Total number of nodes.
  //
  // The total number of nodes that is present in the rank cover tree. This
  // value is expected to be in O(n) whenever the sample rate used in the
  // construction is at least 2.0.
  //
  // Note, while in design the RCT uses copies of nodes at different levels,
  // these copies are not made in the actual implementation! Therefore, the
  // number of nodes becomes a symbolic property of the RCT.
  size_t num_nodes_;

  // Lists of tentative parents per node.
  //
  // The parent lists store the tentative parents which are sought during the
  // RCT construction phase. The number of parent sought depend on the setting
  // of maxParents.
  //
  // TODO: Use managed containers.
  // TODO: Move this into an RCT factory.
  int*** parentIndexLLList;

  // Tentative parent list lengths.
  //
  // These arrays store the number of tentative parents found for each node
  // during the construction phase.
  //
  // TODO: Remove this in favor of 'tentative_parent_list_[].size()'
  // TODO: Move this into an RCT factory.
  int** parentLSizeLList;

  // List of children per node.
  //
  // The children lists store the children that any node has at any given level.
  // The lists represent an 'inversion' of the parent lists used during
  // construction and are used to descend the RCT during search operations.
  //
  // TODO: Use managed containers.
  int*** childIndexLLList;

  // TODO: Remove once obsolete, use 'child_index_list_[].size()' instead!!
  int** childLSizeLList;

  // The current query object.
  //
  // The data item that is subject of the currently processed (most recently
  // processed) query. This marker is also used to identify repeated queries
  // where computation time can be saved.
  //
  // TODO: Move per-query state out of the index.
  const DomainType* query_;

  // Distance cache.
  //
  // Stores distances from the current query item to some data items in the RCT
  // that have been performed during the most recent query.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Use managed container.
  float* distFromQueryList;

  // Item cache.
  //
  // Stores indices of items whose distance to the current query object has been
  // determined during the current query.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Use managed container.
  int* storedDistIndexList;

  // Cache length.
  //
  // Stores the number of items in the current query cache.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Remove once obsolete.
  size_t num_stored_dists_;

  // Number of distance comparisons performed.
  //
  // Stores the number of distance comparisons that have been performed during
  // the most recent RCT operation.
  //
  // TODO: Move per-query state out of the index.
  size_t num_dist_comps_;

  // Level quota list.
  //
  // Stores the maximum number of nodes that are retained in the cover set
  // during a search. The quotas are calculated for each search based on the
  // number of neighbors sought.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Use managed container.
  int* levelQuotaList;

  // Coverage parameter.
  //
  // The coverage parameter determines the amount of items that are sought
  // during search and construction operations. Although a higher value is can
  // increase search accuracy and is more likely to produce a well-formed RCT,
  // it comes at the cost of increased search and construction time complexity.
  double coverage_;

  // Query result items.
  //
  // Stores the indices of the items found during the most recent search
  // operation.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Use managed container.
  int* queryResultIndexList;

  // Query result distances.
  //
  // Stores the query-to-result distances of the items found during the most
  // recent search operation.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Use managed container.
  float* queryResultDistList;

  // Query result size.
  //
  // Stores the number of results found in the most recent search operation.
  //
  // TODO: Move per-query state out of the index.
  size_t query_result_size_;

  // Query result sample size.
  //
  // The number of sample items within which the most recent similarity search
  // was performed.
  //
  // TODO: Move per-query state out of the index.
  size_t query_result_sample_size_;

  // Temporary index list.
  //
  // This list stores nodes visited during the search operation.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Use managed container.
  int* visitedNodeIndexList;

  // Temporary distance list.
  //
  // This list stores the distance to the query item of nodes already visited
  // during the current search.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Use managed container.
  float* tempResultDistList;

  // Temporary result list.
  //
  // This list stores result candidates during the search operation.
  //
  // TODO: Move per-query state out of the index.
  // TODO: Use managed container.
  int* tempResultIndexList;

  // Verbosity level.
  //
  // The verbosity level of the RCT determines the amount of feedback provided
  // to the user. Values equal to or less than 0 result in no feedback at all. A
  // value of 1 makes the RCT report only errors, whereas a value of 2 generates
  // messages for errors as well as progress. Finally, a value of 3 reports
  // errors, progress and debug output.
  //
  // TODO: Use logging libraries VLOG levels instead.
  int verbosity_;

  // Build scale factor.
  //
  // The build scale factor determine the scaling that was used during the
  // construction of the RCT.
  double build_coverage_;

  // Random number generator seed.
  //
  // The seed value used to initialize the random number generator.
  //
  // TODO: Use proper type.
  unsigned long seed_;

  RankCoverTree(const RankCoverTree&) = delete;
  RankCoverTree& operator=(const RankCoverTree&) = delete;
};

}  // namespace rct
}  // namespace indexing

#include "indexing/rct/public/rank_cover_tree-inl.h"

#endif  // INDEXING_RCT_PUBLIC_RANK_COVER_TREE_H_
