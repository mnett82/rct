#include "indexing/rct/public/rank_cover_tree.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>

#include <glog/logging.h>

#if defined(_RCT_TEMPLATE_DECL) || defined(_RCT_TEMPLATE_DEFN)
#error "Macros '_RCT_TEMPLATE_DECL' or '_RCT_TEMPLATE_DEFN' already defined!"
#endif
#define _RCT_TEMPLATE_DECL                              \
  template <typename DomainType, typename DistanceType, \
            DistanceType (*DistanceFn)(const DomainType&, const DomainType&)>
#define _RCT_TEMPLATE_DEFN RankCoverTree<DomainType, DistanceType, DistanceFn>

namespace indexing {
namespace rct {
namespace {

constexpr size_t kNoIndex = std::numeric_limits<size_t>::max();

// TODO: Get rid of these constants.
constexpr int kNone = -1;
constexpr float kUnknown = -1.f;
constexpr int kBufSize = 1024;
constexpr char kVersion[] = "1.0";
constexpr int FALSE = 0;
constexpr int TRUE = 1;

}  // namespace

_RCT_TEMPLATE_DECL
bool _RCT_TEMPLATE_DEFN::Build(const std::vector<const DomainType*>& data,
                               const double build_coverage,
                               const size_t num_parents,
                               const double sample_probability) {
  // TODO: The RCT should handle these cases gracefully and not as an error.
  if (data.empty()) {
    LOG(ERROR) << "RCT::Build() called with empty data set.";
    return false;
  }
  if (data.size() == 1) {
    LOG(ERROR) << "RCT::Build() called with singleton data set.";
    return false;
  }

  VLOG(2) << "Building RCT from " << data.size() << " data elements...";
  data_ = data;
  sample_rate_ = 1.0 / sample_probability;
  SetupLevels(data_.size(), num_parents, sample_probability);
  AllocateStorage();

  // Randomly assign data items to RCT levels.
  std::random_shuffle(
      intern_to_extern_mapping_.begin(), intern_to_extern_mapping_.end(),
      [this](const size_t range) {
        return std::uniform_int_distribution<size_t>(0, range - 1)(rand_);
      });

  // Build the index.
  num_dist_comps_ = 0;
  build_coverage_ = build_coverage;
  BuildIncrementally();
  PrintStats();
  return true;  // TODO: should not return anything.
}

_RCT_TEMPLATE_DECL
void _RCT_TEMPLATE_DEFN::SetNewQuery(const DomainType* query) {
  if (query != query_) {
    for (int i = 0; i < num_stored_dists_; ++i) {
      distFromQueryList[storedDistIndexList[i]] = kUnknown;
    }
    num_stored_dists_ = 0;
    query_ = query;
  }
}

_RCT_TEMPLATE_DECL
void _RCT_TEMPLATE_DEFN::PrintStats() const {
  VLOG(2) << "RCT build statistics:";
  VLOG(2) << "  size                   == " << data_.size();
  VLOG(2) << "  levels                 == " << levels;
  VLOG(2) << "  total nodes            == " << num_nodes_;
  VLOG(2) << "  max parents per node   == " << max_parents_;
  VLOG(2) << "  max node degree        == " << max_degree_;
  VLOG(2) << "  avg node degree        == " << average_degree_;
  VLOG(2) << "  distance comparisons   == " << num_dist_comps_;
  VLOG(2) << "  RNG seed               == " << seed_;
}

// TODO: Extract a function like BuildLevel(i) and call that in a loop.
_RCT_TEMPLATE_DECL
void _RCT_TEMPLATE_DEFN::BuildIncrementally() {
  // Build the top level of the RCT as a special case by promoting the first
  // data element to the artificial root node of the RCT.
  parentLSizeLList[levels][0] = 0;
  parentIndexLLList[levels][0] = NULL;
  size_t num_lower_items = levelSetSizeList[levels - 1];
  size_t num_upper_items = 0;
  childLSizeLList[levels][0] = num_lower_items;
  childIndexLLList[levels][0] = new int[num_lower_items];
  max_degree_ = num_lower_items;
  size_t total_degree = num_lower_items;
  for (int i = 0; i < num_lower_items; i++) {
    childIndexLLList[levels][0][i] = i;
  }
  VLOG(2) << "RCT root level constructed.";

  for (int lvl = levels - 2; lvl >= 0; lvl--) {
    num_upper_items = levelSetSizeList[lvl + 1];
    num_lower_items = levelSetSizeList[lvl];

    // Connect the bottom level of the current partial RCT to the items on this
    // level. The elements on the previously constructed level (lvl + 1) are to
    // become parents of the items on the current level.
    //
    // Also, temporarily store (in "childLSizeLList") the number of times each
    // node is requested as a parent.
    for (int child = 0; child < num_lower_items; child++) {
      if (child % 5000 == 4999) {
        VLOG(2) << "Inserting item " << child + 1 << " (out of "
                << num_lower_items << ") at level " << lvl << "...";
      }

      // Find some parents for the current child. If only one parent is
      // requested and the child has a copy at the level above, then just choose
      // it directly. Otherwise, do a search with respect to the partial RCT.
      if ((child < num_upper_items) && (max_parents_ == 1)) {
        query_result_size_ = 1;
        queryResultIndexList[0] = child;
      } else {
        SetNewQuery(data_[intern_to_extern_mapping_[child]]);
        // TODO: doFindNear(max_parents_, lvl + 1, buildScaleFactor);
      }

      // Connect links from child to parents. If a copy of the child also exists
      // at the upper level, then make sure that it is listed as the first
      // parent.
      parentLSizeLList[lvl][child] = query_result_size_;
      parentIndexLLList[lvl][child] = new int[max_parents_];
      int offset = 0;
      if ((child < num_upper_items) && (queryResultIndexList[0] != child)) {
        // The first query result should have been a copy of the child, but
        // wasn't. Repair this situation by explicitly placing a copy of the
        // child at the head of the list, and shifting the remaining query
        // result elements to accommodate the child copy.
        offset = 1;
        parentIndexLLList[lvl][child][0] = child;
        childLSizeLList[lvl + 1][child]++;
      } else {
        // Either the query result contains a copy of the child at its head, or
        // the child isn't supposed to appear in the query result anyway.
        offset = 0;
        parentIndexLLList[lvl][child][0] = queryResultIndexList[0];
        childLSizeLList[lvl + 1][queryResultIndexList[0]]++;
      }
      for (int i = 1; i < query_result_size_; i++) {
        // Apply an offset shift only until a copy of the child is found (one
        // may not necessarily be found). This is to avoid picking up this copy
        // more than once.
        if (queryResultIndexList[i] == child) {
          offset = 0;
        } else {
          parentIndexLLList[lvl][child][i] = queryResultIndexList[i - offset];
          childLSizeLList[lvl + 1][queryResultIndexList[i - offset]]++;
        }
      }
    }

    // For each parent, reserve storage for its child lists.
    for (int parent = 0; parent < num_upper_items; parent++) {
      int childLSize = childLSizeLList[lvl + 1][parent];

      if (childLSize > 0) {
        childIndexLLList[lvl + 1][parent] =
            new int[childLSizeLList[lvl + 1][parent]];
        childLSizeLList[lvl + 1][parent] = 0;

        total_degree += (long int)childLSize;

        if (childLSize > max_degree_) {
          max_degree_ = childLSize;
        }
      }
    }

    // Construct child lists for each of the parents, by reversing the
    // child-to-parent edges. Since the child-to-parent edges are no longer
    // needed, delete them.
    for (int child = 0; child < num_lower_items; child++) {
      for (int i = parentLSizeLList[lvl][child] - 1; i >= 0; i--) {
        int parent = parentIndexLLList[lvl][child][i];
        int j = childLSizeLList[lvl + 1][parent];
        childIndexLLList[lvl + 1][parent][j] = child;
        childLSizeLList[lvl + 1][parent]++;
      }

      if (parentLSizeLList[lvl][child] > 0) {
        delete[] parentIndexLLList[lvl][child];
        parentIndexLLList[lvl][child] = NULL;
        parentLSizeLList[lvl][child] = 0;
      }
    }

    // The RCT has grown by one level.
    VLOG(2) << "RCT level " << lvl << "constructed.";
  }

  average_degree_ =
      (float)(((double)total_degree) / (num_nodes_ - data_.size()));
}

_RCT_TEMPLATE_DECL
void _RCT_TEMPLATE_DEFN::SetupLevels(const size_t num_items,
                                     const size_t num_parents,
                                     const double sample_probability) {
  CHECK_LT(0, num_parents);
  max_parents_ = num_parents;

  // Determine the maximum level that each individual data item occurs at.
  std::vector<size_t> max_level_list(num_items);
  std::geometric_distribution<size_t> dist(1.0 - sample_probability);
  std::generate(max_level_list.begin(), max_level_list.end(),
                [this, &dist]() { return dist(rand_); });
  levels = *std::max_element(max_level_list.begin(), max_level_list.end());

  // Allocate storage for various containers.
  levelQuotaList = new int[levels + 1];
  levelSetSizeList = new int[levels + 1];

  // Generate the sizes of individual levels in the RCT. The top-most level
  // contains only a single artificially promoted node.
  levelQuotaList = new int[levels + 1];
  levelSetSizeList = new int[levels + 1];
  for (int lvl = 0; lvl < levels; lvl++) {
    levelQuotaList[lvl] = 0;
    levelSetSizeList[lvl] = 0;
  }
  levelQuotaList[levels] = 0;
  levelSetSizeList[levels] = 1;
  for (int i = 0; i < data_.size(); i++) {
    for (int lvl = max_level_list[i]; lvl >= 0; lvl--) {
      levelSetSizeList[lvl]++;
    }
  }

  // Calculate the number of nodes in the RCT.
  num_nodes_ = 0;
  for (int lvl = 0; lvl < levels; ++lvl) {
    num_nodes_ += levelSetSizeList[lvl];
  }
}

_RCT_TEMPLATE_DECL
void _RCT_TEMPLATE_DEFN::AllocateStorage() {
  // Reserve storage for the mapping between internal and external data indices.
  intern_to_extern_mapping_.resize(data_.size());
  for (int i = 0; i < data_.size(); i++) {
    intern_to_extern_mapping_[i] = i;
  }

  // Set up storage for child-to-parent edges and parent-to-child edges.
  parentIndexLLList = new int**[levels + 1];
  parentLSizeLList = new int*[levels + 1];
  childIndexLLList = new int**[levels + 1];
  childLSizeLList = new int*[levels + 1];
  for (int lvl = 0; lvl <= levels; lvl++) {
    parentIndexLLList[lvl] = new int*[levelSetSizeList[lvl]];
    parentLSizeLList[lvl] = new int[levelSetSizeList[lvl]];
    childIndexLLList[lvl] = new int*[levelSetSizeList[lvl]];
    childLSizeLList[lvl] = new int[levelSetSizeList[lvl]];
    for (int i = levelSetSizeList[lvl] - 1; i >= 0; i--) {
      parentIndexLLList[lvl][i] = NULL;
      parentLSizeLList[lvl][i] = 0;
      childIndexLLList[lvl][i] = NULL;
      childLSizeLList[lvl][i] = 0;
    }
  }

  // Set up storage for managing distance computations and query results.
  distFromQueryList = new float[data_.size()];
  storedDistIndexList = new int[data_.size()];
  num_stored_dists_ = 0;
  queryResultDistList = new float[data_.size()];
  queryResultIndexList = new int[data_.size()];
  query_result_size_ = 0;
  query_result_sample_size_ = 0;
  visitedNodeIndexList = new int[data_.size()];
  tempResultIndexList = new int[data_.size()];
  tempResultDistList = new float[data_.size()];
  for (int i = 0; i < data_.size(); i++) {
    distFromQueryList[i] = kUnknown;
    storedDistIndexList[i] = kNone;
    queryResultDistList[i] = kUnknown;
    queryResultIndexList[i] = kNone;
    visitedNodeIndexList[i] = FALSE;
    tempResultIndexList[i] = kNone;
    tempResultDistList[i] = kUnknown;
  }
}

_RCT_TEMPLATE_DECL
DistanceType _RCT_TEMPLATE_DEFN::ComputeDistFromQuery(const int index) {
  if (distFromQueryList[index] == kUnknown) {
    distFromQueryList[index] =
        DistanceFn(*query_, *data_[intern_to_extern_mapping_[index]]);
    storedDistIndexList[num_stored_dists_] = index;
    ++num_stored_dists_;
    ++num_dist_comps_;
  }
  return distFromQueryList[index];
}

_RCT_TEMPLATE_DECL
int _RCT_TEMPLATE_DEFN::FindNear(const int how_many, const int sample_level,
                                 const double coverage) {
  // Compute quota of items to be retained at every level.
  // Rank cover tree rules.levelQuotaList
  double var_quota = (double)how_many;
  for (int lvl = sample_level; lvl < levels; lvl++) {
    levelQuotaList[lvl] = (int)((coverage * var_quota) + 0.999999F);
    if (levelQuotaList[lvl] < coverage) {
      levelQuotaList[lvl] = (int)(coverage + 0.999999F);
    }
    var_quota /= sample_rate_;
  }
  if (how_many > levelQuotaList[sample_level]) {
    levelQuotaList[sample_level] = how_many;
  }

  // Load the root as the tentative sole member of the query result list.
  query_result_size_ = 0;
  queryResultDistList[0] = ComputeDistFromQuery(0);
  queryResultIndexList[0] = 0;
  int num_retained = 1;

  // From the root, search out other nodes to place in the query result.
  for (int lvl = levels - 1; lvl >= sample_level; lvl--) {
    // For every node at the active level, load its children
    //   into the scratch list, and compute their distances to the query.
    int num_found = 0;

    for (int i = 0; i < num_retained; i++) {
      int node_index = queryResultIndexList[i];
      int num_children = childLSizeLList[lvl + 1][node_index];
      int* childList = childIndexLLList[lvl + 1][node_index];

      for (int j = 0; j < num_children; j++) {
        int child = childList[j];

        if (visitedNodeIndexList[child] != TRUE) {
          visitedNodeIndexList[child] = TRUE;
          tempResultIndexList[num_found] = child;
          tempResultDistList[num_found] = ComputeDistFromQuery(child);
          num_found++;
        }
      }
    }

    for (int i = 0; i < num_found; i++) {
      visitedNodeIndexList[tempResultIndexList[i]] = FALSE;
    }

    // Extract the closest nodes from the list of accumulated children,
    //   and keep them as the tentative parents of the query.

    if (num_found > levelQuotaList[lvl]) {
      num_retained = levelQuotaList[lvl];
    } else {
      num_retained = num_found;
    }

    num_retained = PartialQuickSort(num_retained, tempResultDistList,
                                    tempResultIndexList, 0, num_found - 1);

    for (int i = 0; i < num_retained; i++) {
      queryResultIndexList[i] = tempResultIndexList[i];
      queryResultDistList[i] = tempResultDistList[i];
    }
  }

  // Select the final number of neighbors needed.
  if (num_retained > how_many) {
    query_result_size_ = how_many;
  } else {
    query_result_size_ = num_retained;
  }
  return query_result_size_;
}

// =================[ DO NOT CLEAN UP! WILL BE REMOVED SOON! ]=================
_RCT_TEMPLATE_DECL
int _RCT_TEMPLATE_DEFN::PartialQuickSort(int how_many, float* distList,
                                         int* indexList, int rangeFirst,
                                         int rangeLast) {
  int i;
  int pivotLoc = 0;
  int pivotIndex = 0;
  float pivotDist = 0.0F;
  int tempIndex = 0;
  float tempDist = 0.0F;
  int low = 0;
  int high = 0;
  int num_found = 0;
  int numDuplicatesToReplace = 0;
  int tieBreakIndex = 0;

  // If the range is empty, or if we've been asked to sort no
  //   items, then return immediately.

  if ((rangeLast < rangeFirst) || (how_many < 1)) {
    return 0;
  }

  // If there is exactly one element, then again there is nothing
  //   that need be done.

  if (rangeLast == rangeFirst) {
    return 1;
  }

  // If the range to be sorted is small, just do an insertion sort.

  if (rangeLast - rangeFirst < 7) {
    std::uniform_int_distribution<size_t> distribution(0,
                                                       rangeLast - rangeFirst);
    high = rangeFirst + 1;
    tieBreakIndex = indexList[distribution(rand_)];

    // The outer while loop considers each item in turn (starting
    //   with the second item in the range), for insertion into
    //   the sorted list of items that precedes it.

    while (high <= rangeLast) {
      // Copy the next item to be inserted, as the "pivot".
      // Start the insertion tests with its immediate predecessor.

      pivotDist = distList[high];
      pivotIndex = indexList[high];
      low = high - 1;

      // Work our way down through previously-sorted items
      //   towards the start of the range.

      while (low >= rangeFirst) {
        // Compare the item to be inserted (the "pivot") with
        //   the current item.

        if (distList[low] < pivotDist) {
          // The current item precedes the pivot in the sorted order.
          // Break out of the loop - we have found the insertion point.

          break;
        } else if (distList[low] > pivotDist) {
          // The current item follows the pivot in the sorted order.
          // Shift the current item one spot upwards, to make room
          //   for inserting the pivot below it.

          distList[low + 1] = distList[low];
          indexList[low + 1] = indexList[low];
          low--;
        } else {
          if (indexList[low] != pivotIndex) {
            // The items have the same sort value but are not identical.
            // Break the tie pseudo-randomly.

            if (((tieBreakIndex >= pivotIndex) &&
                 ((indexList[low] < pivotIndex) ||
                  (tieBreakIndex < indexList[low]))) ||
                ((tieBreakIndex < pivotIndex) &&
                 ((indexList[low] < pivotIndex) &&
                  (tieBreakIndex < indexList[low])))) {
              // The current item precedes the pivot in the sorted order.
              // Break out of the loop - we have found the insertion point.

              break;
            } else {
              // The current item follows the pivot in the sorted order.
              // Shift the current item one spot upwards, to make room
              //   for inserting the pivot below it.

              distList[low + 1] = distList[low];
              indexList[low + 1] = indexList[low];
              low--;
            }
          } else {
            // Oh no!
            // We opened up an empty slot for the pivot,
            //   only to find that it's a duplicate of the current item!
            // Close the slot up again, and eliminate the duplicate.

            for (i = low + 1; i < high; i++) {
              distList[i] = distList[i + 1];
              indexList[i] = indexList[i + 1];
            }

            // To eliminate the duplicate, overwrite its location with the
            //   item from the end of the range, and then shrink the range
            //   by one.

            distList[high] = distList[rangeLast];
            indexList[high] = indexList[rangeLast];
            rangeLast--;

            // The next iteration must not advance "high", since we've
            //   just put a new element into it which needs to be processed.
            // Decrementing it here will cancel out with the incrementation
            //   of the next iteration.

            high--;

            // When we break the loop, the pivot element will be put
            //   in its proper place ("low" + 1)
            // Here, the proper place is where rangeLast used to be.
            // To achieve this, we need to adjust "low" here.

            low = rangeLast;

            break;
          }
        }
      }

      // If we've made it to here, we've found the insertion
      //   spot for the current element.
      // Perform the insertion.

      low++;
      distList[low] = pivotDist;
      indexList[low] = pivotIndex;

      // Move to the next item to be inserted in the growing sorted list.

      high++;
    }

    // Return the number of sorted items found.

    num_found = rangeLast - rangeFirst + 1;

    if (num_found > how_many) {
      num_found = how_many;
    }

    return num_found;
  }

  // The range to be sorted is large, so do a partial quicksort.
  // Select a pivot item, and swap it with the item at the beginning
  //   of the range.

  std::uniform_int_distribution<size_t> distribution(0, rangeLast - rangeFirst);
  pivotLoc = rangeFirst + distribution(rand_);
  tieBreakIndex = indexList[distribution(rand_)];

  pivotDist = distList[pivotLoc];
  distList[pivotLoc] = distList[rangeFirst];
  distList[rangeFirst] = pivotDist;

  pivotIndex = indexList[pivotLoc];
  indexList[pivotLoc] = indexList[rangeFirst];
  indexList[rangeFirst] = pivotIndex;

  // Eliminate all duplicates of the pivot.
  // Any duplicates found are pushed to the end of the range, and
  //   the range shrunk by one (thereby excluding them).

  i = rangeFirst + 1;

  while (i <= rangeLast) {
    if ((pivotIndex == indexList[i]) && (pivotDist == distList[i])) {
      distList[i] = distList[rangeLast];
      indexList[i] = indexList[rangeLast];
      rangeLast--;
    } else {
      i++;
    }
  }

  // Partition the remaining items with respect to the pivot.
  // This efficient method is adapted from the one outlined in
  //   Cormen, Leiserson & Rivest.
  // The range is scanned from both ends.
  // Items with small distances are placed below "low", and those
  //   with large distances are placed above "high".
  // Where "low" and "high" meet, the pivot item is inserted.

  low = rangeFirst;
  high = rangeLast + 1;

  while (TRUE) {
    // Move the "high" endpoint down until it meets either the pivot,
    //   or something that belongs on the "low" side.
    // If the key values are tied, decide pseudo-randomly.

    do {
      high--;
    } while ((distList[high] > pivotDist) ||
             ((distList[high] == pivotDist) && (high > low) &&
              (((tieBreakIndex >= pivotIndex) &&
                ((pivotIndex < indexList[high]) &&
                 (indexList[high] <= tieBreakIndex))) ||
               ((tieBreakIndex < pivotIndex) &&
                ((pivotIndex < indexList[high]) ||
                 (indexList[high] <= tieBreakIndex))))));

    // Move the "low" endpoint up until it meets either the "high" endpoint,
    //   or something that belongs on the "high" side.
    // If the key values are tied, decide pseudo-randomly.

    do {
      low++;
    } while ((low < high) && ((distList[low] < pivotDist) ||
                              ((distList[low] == pivotDist) &&
                               (((tieBreakIndex >= pivotIndex) &&
                                 ((indexList[low] <= pivotIndex) ||
                                  (tieBreakIndex < indexList[low]))) ||
                                ((tieBreakIndex < pivotIndex) &&
                                 ((indexList[low] <= pivotIndex) &&
                                  (tieBreakIndex < indexList[low])))))));

    // Have the "low" and "high" endpoints crossed?
    // If not, we still have more work to do.

    if (low < high) {
      // Swap the misplaced items, and try again.

      tempDist = distList[low];
      distList[low] = distList[high];
      distList[high] = tempDist;

      tempIndex = indexList[low];
      indexList[low] = indexList[high];
      indexList[high] = tempIndex;
    } else {
      // We found the cross-over point.

      break;
    }
  }

  // The pivot value ends up at the location referenced by "high".
  // Swap it with the pivot (which resides at the beginning of the range).

  distList[rangeFirst] = distList[high];
  distList[high] = pivotDist;

  indexList[rangeFirst] = indexList[high];
  indexList[high] = pivotIndex;

  pivotLoc = high;

  // The partition is complete.
  // Recursively sort the items with smaller distance.

  num_found =
      PartialQuickSort(how_many, distList, indexList, rangeFirst, pivotLoc - 1);

  // If we found enough items (including the pivot), then we are done.
  // Make sure the pivot is in its correct position, if it is used.

  if (num_found >= how_many - 1) {
    if (num_found == how_many - 1) {
      distList[rangeFirst + num_found] = pivotDist;
      indexList[rangeFirst + num_found] = pivotIndex;
    }

    return how_many;
  }

  // We didn't find enough items, even taking the pivot into account.
  // Were any duplicates discovered during this call?

  if (num_found < pivotLoc - rangeFirst) {
    // Duplicates were discovered!
    // Figure out the minimum number of duplicates that must be
    //   replaced by items from the end of the range in order to
    //   leave the non-duplicates in contiguous locations.

    numDuplicatesToReplace = pivotLoc - rangeFirst - num_found;
    high = rangeLast;

    if (numDuplicatesToReplace > rangeLast - pivotLoc) {
      numDuplicatesToReplace = rangeLast - pivotLoc;
      rangeLast = rangeFirst + num_found + numDuplicatesToReplace;
    } else {
      rangeLast -= numDuplicatesToReplace;
    }

    // Replace the required number of duplicates by items from
    //   the end of the range.
    // The size of the range will shrink as a result.

    low = rangeFirst + num_found + 1;

    for (i = 0; i < numDuplicatesToReplace; i++) {
      distList[low] = distList[high];
      indexList[low] = indexList[high];
      low++;
      high--;
    }
  }

  // Put the pivot element in its proper place.

  distList[rangeFirst + num_found] = pivotDist;
  indexList[rangeFirst + num_found] = pivotIndex;

  // Finish up by sorting larger-distance items.
  // Note that the number of sorted items needed has dropped.

  return num_found + 1 + PartialQuickSort(how_many - num_found - 1, distList,
                                          indexList, rangeFirst + num_found + 1,
                                          rangeLast);
}

}  // namespace rct
}  // namespace indexing

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// #ifndef kNone
// #define kNone (-1)
// #endif

// #ifndef kUnknown
// #define kUnknown (-1.0F)
// #endif

// #ifndef kVersion
// #define kVersion ("1.0")
// #endif

// /*!
//  * Sets the sample rate.
//  *
//  * @note Must be called before construction.
//  *
//  * @param sample_rate_ The desired sample rate.
//  */
// _RCT_TMPL_DECL
// void RCT::setSampleRate(const float& sample_rate_) {
//   (*this).sample_rate_ = sample_rate_;
// }

// /*!
//  * Constructor using seed for random number generator initialization.
//  */
// _RCT_TMPL_DECL
// RCT::RCT(const unsigned long& seed) : seed(seed) {
//   data = NULL;
//   size = 0;
//   max_parents_ = 1;
//   max_degree_ = 0;
//   avgDegree = 0.0F;
//   intern_to_extern_mapping_ = NULL;
//   levelSetSizeList = NULL;
//   levels = 0;
//   num_nodes_ = 0;
//   parentIndexLLList = NULL;
//   parentLSizeLList = NULL;
//   childIndexLLList = NULL;
//   childLSizeLList = NULL;
//   query = NULL;
//   distFromQueryList = NULL;
//   storedDistIndexList = NULL;
//   num_stored_dists_ = 0;
//   numDistComps = 0UL;
//   levelQuotaList = NULL;
//   coverageParameter = 1.0f;
//   queryResultIndexList = NULL;
//   queryResultDistList = NULL;
//   visitedNodeIndexList = NULL;
//   tempResultIndexList = NULL;
//   tempResultDistList = NULL;
//   query_result_size_ = 0;
//   query_result_sample_size_ = 0;
//   verbosity = 0;
//   coverageParameter = 1.0f;
//   sample_rate_ = 2.0f;
//   buildScaleFactor = kUnknown;
//   rand_.seed(seed);
// }

// /*!
//  * The destructor releases memory allocated by the RCT prior to its
//  * destruction.
//  */
// _RCT_TMPL_DECL
// RCT::~RCT() {
//   int i;
//   int lvl;
//   int tempLength = 0;
//   data = NULL;
//   if (intern_to_extern_mapping_ != NULL) {
//     delete[] intern_to_extern_mapping_;
//     intern_to_extern_mapping_ = NULL;
//   }

//   if (parentIndexLLList != NULL) {
//     for (lvl = 0; lvl <= levels; lvl++) {
//       if (parentIndexLLList[lvl] != NULL) {
//         tempLength = levelSetSizeList[lvl];

//         for (i = 0; i < tempLength; i++) {
//           if (parentIndexLLList[lvl][i] != NULL) {
//             delete[] parentIndexLLList[lvl][i];
//             parentIndexLLList[lvl][i] = NULL;
//           }
//         }

//         delete[] parentIndexLLList[lvl];
//         parentIndexLLList[lvl] = NULL;
//       }
//     }

//     delete[] parentIndexLLList;
//     parentIndexLLList = NULL;
//   }

//   if (parentLSizeLList != NULL) {
//     for (lvl = 0; lvl <= levels; lvl++) {
//       if (parentLSizeLList[lvl] != NULL) {
//         delete[] parentLSizeLList[lvl];
//         parentLSizeLList[lvl] = NULL;
//       }
//     }

//     delete[] parentLSizeLList;
//     parentLSizeLList = NULL;
//   }

//   if (childIndexLLList != NULL) {
//     for (lvl = 0; lvl <= levels; lvl++) {
//       if (childIndexLLList[lvl] != NULL) {
//         tempLength = levelSetSizeList[lvl];

//         for (i = 0; i < tempLength; i++) {
//           if (childIndexLLList[lvl][i] != NULL) {
//             delete[] childIndexLLList[lvl][i];
//             childIndexLLList[lvl][i] = NULL;
//           }
//         }

//         delete[] childIndexLLList[lvl];
//         childIndexLLList[lvl] = NULL;
//       }
//     }

//     delete[] childIndexLLList;
//     childIndexLLList = NULL;
//   }

//   if (childLSizeLList != NULL) {
//     for (lvl = 0; lvl <= levels; lvl++) {
//       if (childLSizeLList[lvl] != NULL) {
//         delete[] childLSizeLList[lvl];
//         childLSizeLList[lvl] = NULL;
//       }
//     }

//     delete[] childLSizeLList;
//     childLSizeLList = NULL;
//   }

//   if (levelSetSizeList != NULL) {
//     delete[] levelSetSizeList;
//     levelSetSizeList = NULL;
//   }

//   query = NULL;

//   if (distFromQueryList != NULL) {
//     delete[] distFromQueryList;
//     distFromQueryList = NULL;
//   }

//   if (storedDistIndexList != NULL) {
//     delete[] storedDistIndexList;
//     storedDistIndexList = NULL;
//   }

//   if (levelQuotaList != NULL) {
//     delete[] levelQuotaList;
//     levelQuotaList = NULL;
//   }

//   if (queryResultIndexList != NULL) {
//     delete[] queryResultIndexList;
//     queryResultIndexList = NULL;
//   }

//   if (queryResultDistList != NULL) {
//     delete[] queryResultDistList;
//     queryResultDistList = NULL;
//   }

//   if (tempResultIndexList != NULL) {
//     delete[] tempResultIndexList;
//     tempResultIndexList = NULL;
//   }

//   if (tempResultDistList != NULL) {
//     delete[] tempResultDistList;
//     tempResultDistList = NULL;
//   }

//   if (visitedNodeIndexList != NULL) {
//     delete[] visitedNodeIndexList;
//     visitedNodeIndexList = NULL;
//   }
// }

// /*!
//  * Loads a previously-computed RCT from the specified file. The original data
//  * set must also be provided (as well as the number of items in the data
//  set).
//  *
//  * @note The extension ".rctf" is automatically appended to the file name.
//  *
//  * @return If successful, the number of RCT items is returned. Otherwise,
//  zero
//  *         is returned.
//  */
// _RCT_TMPL_DECL
// int RCT::build(const char* fileName, DistData** inputData,
//                const int& numItems) {
//   int i = 0;
//   int j = 0;
//   int lvl = 0;
//   int loc = 0;
//   int inLevel = 0;
//   // int temp = 0;
//   int inSize = 0;
//   int inLevels = 0;
//   int inNumNodes = 0;
//   int inMaxParents = 0;
//   int inMaxDegree = 0;
//   float inAvgDegree = 0.0f;
//   float inCoverageParameter = 0.0f;
//   float inBuildScaleFactor = 0.0f;
//   int num_children = 0;
//   int* childList = NULL;
//   ifstream inFile;
//   int levelSetSize = 0;

//   // If the data set is empty, then abort.
//   if ((fileName == NULL) || (numItems <= 0) || (inputData == NULL)) {
//     if (verbosity > 0) {
//       if (numItems == 1) {
//         LOG(ERROR) << "Data set has only 1 item.";
//       } else {
//         LOG(ERROR) << "Empty data set or filename.";
//       }
//     }
//     return 0;
//   }

//   if (verbosity >= 2) {
//     LOG(INFO) << "Loading RCT from file '" << fileName << ".rctf'...";
//   }
//   data = inputData;

//   // Open the file containing the RCT. If we fail to open the file, abort.
//   ostringstream fullFileName;
//   fullFileName << fileName << ".rctf";
//   inFile.open(fullFileName.str().c_str(), ios::in);
//   if (!inFile.is_open()) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "File '" << fullFileName.str() << "' could not be
//       opened.";
//     }
//     return 0;
//   }

//   // Skip two comment lines.
//   string buffer;
//   getline(inFile, buffer);
//   assert((buffer[0] == '%') && (buffer[1] == '%'));
//   getline(inFile, buffer);
//   assert((buffer[0] == '%') && (buffer[1] == '%'));

//   // Read in basic parameters.
//   inFile >> inSize >> inLevels >> inNumNodes >> inMaxParents >> inMaxDegree;
//   inFile >> inAvgDegree >> inCoverageParameter >> inBuildScaleFactor;
//   inFile >> seed;

//   // Are these parameter values what we expected? If not, then abort!
//   if (inSize != numItems) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Unexpected RCT parameters in file '" <<
//       fullFileName.str()
//                  << "'.";
//     }
//     inFile.close();
//     return 0;
//   }
//   getline(inFile, buffer);

//   // Assign properties.
//   size = inSize;
//   levels = inLevels;
//   num_nodes_ = inNumNodes;
//   max_parents_ = inMaxParents;
//   max_degree_ = inMaxDegree;
//   avgDegree = inAvgDegree;
//   coverageParameter = inCoverageParameter;
//   buildScaleFactor = inBuildScaleFactor;

//   // Skip another comment line.
//   getline(inFile, buffer);
//   assert((buffer[0] == '%') && (buffer[1] == '%'));

//   // Fetch the level set sizes.
//   levelSetSizeList = new int[levels + 1];
//   levelQuotaList = new int[levels + 1];
//   for (lvl = 0; lvl <= levels; lvl++) {
//     inFile >> levelSetSizeList[lvl];
//     levelQuotaList[lvl] = 0;
//   }
//   getline(inFile, buffer);

//   // Reserve RCT storage. After this operation, the expected RCT size,
//   // number of levels, etc, are set.
//   reserveStorage();

//   // Skip yet another comment line.
//   getline(inFile, buffer);
//   assert((buffer[0] == '%') && (buffer[1] == '%'));

//   // Read in information for each node: level set, internal index,
//   // external index, number of children, and indices of children.
//   // Build the list of children, if any exist.
//   for (lvl = levels; lvl >= 0; lvl--) {
//     levelSetSize = levelSetSizeList[lvl];
//     for (i = 0; i < levelSetSize; i++) {
//       inFile >> inLevel >> loc;
//       if ((loc != i) || (lvl != inLevel)) {
//         if (verbosity > 0) {
//           LOG(ERROR) << "Invalid entry in file '" << fileName << ".rctf'.";
//         }
//         inFile.close();
//         return 0;
//       }
//       inFile >> intern_to_extern_mapping_[i] >> num_children;
//       if (num_children > 0) {
//         childList = new int[num_children];
//       } else {
//         childList = NULL;
//       }
//       childIndexLLList[lvl][i] = childList;
//       childLSizeLList[lvl][i] = num_children;
//       for (j = 0; j < num_children; j++) {
//         inFile >> childList[j];
//       }
//       getline(inFile, buffer);
//     }
//   }
//   return size;
// }

// /*!
//  * Perform an approximate nearest-neighbor query for the specified item. The
//  * number of desired nearest neighbors <code>how_many</code> (default 1)
//  * can be specified. The search is relative to the random level
//  * sample_level (default is 0, i.e. the complete
//  * data set). The method also makes use of a parameter (<code>coverage
//  * </code>) that influences the trade-off between time and accuracy.
//  *
//  * @param coverage Additional search efforts.
//  * @param query The query location.
//  * @param how_many The desired number of neighbors.
//  * @param sample_level The sample level with respect to which the query is
//  *                    performed.
//  *
//  * @return The number of elements actually found.
//  *
//  * @note The query result can be obtained via calls to the following methods:
//  *       <code>getResultAcc</code>, <code>getResultDists</code>,
//  *       <code>getResultDistComps</code>, <code>getResultIndices</code> and
//  *       <code>getResultNumFound</code>. The result items are sorted in
//  *       increasing order of their distances to the query.
//  */
// _RCT_TMPL_DECL
// int RCT::findNear(DistData* query, const int& how_many, const float&
// coverage,
//                   const int& sample_level) {
//   query_result_size_ = 0;
//   query_result_sample_size_ = 0;
//   numDistComps = 0UL;
//   if ((size <= 0) || (query == NULL) || (how_many <= 0) || (sample_level < 0)
//   ||
//       ((sample_level >= levels) && (size > 1)) || (coverage <= 0.0F)) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Invalid argument(s).";
//     }
//     return 0;
//   }
//   setNewQuery(query);
//   return doFindNear(how_many, sample_level, coverage);
// }

// /*!
//  * Perform an exact nearest-neighbor query for the specified item. The
//  * number of desired nearest neighbors <code>how_many</code> (default 1)
//  * can be specified. The search is relative to the random level
//  * sample_level (default is 0, i.e. the complete
//  * data set).
//  *
//  * @param query The query location.
//  * @param how_many The desired number of neighbors.
//  * @param sample_level The sample level with respect to which the query is
//  *                    performed.
//  *
//  * @return The number of elements actually found.
//  *
//  * @note The query result can be obtained via calls to the following methods:
//  *       <code>getResultAcc</code>, <code>getResultDists</code>,
//  *       <code>getResultDistComps</code>, <code>getResultIndices</code> and
//  *       <code>getResultNumFound</code>. The result items are sorted in
//  *       increasing order of their distances to the query.
//  */
// _RCT_TMPL_DECL
// int RCT::findNearest(DistData* query, const int& how_many,
//                      const int& sample_level) {
//   query_result_size_ = 0;
//   query_result_sample_size_ = 0;
//   numDistComps = 0UL;
//   if ((size <= 0) || (query == NULL) || (how_many <= 0) || (sample_level < 0)
//   ||
//       ((sample_level >= levels) && (size > 1))) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Invalid argument(s).";
//     }
//     return 0;
//   }
//   setNewQuery(query);
//   return doFindNearest(how_many, sample_level);
// }

// /*!
//  * Retrieve the average degree of a node in the RCT. This value is roughly
//  * equal to sample_rate_.
//  *
//  * @return The average degree of a node in the RCT.
//  */
// _RCT_TMPL_DECL
// float RCT::getAvgDegree() const { return avgDegree; }

// /*!
//  * Retrieve the scale factor used during the consturction of the RCT.
//  *
//  * @return The scale factor used during RCT construction.
//  */
// _RCT_TMPL_DECL
// float RCT::getBuildScaleFactor() const { return buildScaleFactor; }

// /*!
//  * Retrieve the coverage parameter set for the RCT.
//  *
//  * @return The coverage parameter of the RCT.
//  */
// _RCT_TMPL_DECL
// float RCT::getCoverageParameter() const { return coverageParameter; }

// /*!
//  * Retrieve a pointer to the data array used to assemble the
//  * RCT.
//  *
//  * @note Since the random leveling process only permutes the
//  *       data items internally, be aware that only the item
//  *       indices retrieved by the <code>getResultIndices</code>
//  *       function refer to external (original) indices. All
//  *       internally used indices reffer to the permuted
//  *       ordering of data items.
//  *
//  * @return A pointer to the data items in the RCT.
//  */
// _RCT_TMPL_DECL
// DistData** RCT::getData() { return data; }

// /*!
//  * Fills the supplied list with the mapping from external item indices to
//  * internal RCT indices.
//  *
//  * @param result An array to store the mapping in.
//  * @param capacity The capacity of the target array.
//  *
//  * @return If successful, the number of RCT items is returned. Otherwise,
//  zero
//  *         is returned.
//  */
// _RCT_TMPL_DECL
// int RCT::getExternToInternMapping(int* result, int capacity) const {
//   int i;
//   if ((result == NULL) || (capacity < size)) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Result list capacity is too small.";
//     }
//     return 0;
//   }
//   for (i = 0; i < size; i++) {
//     result[intern_to_extern_mapping_[i]] = i;
//   }
//   return size;
// }

// /*!
//  * Fills the supplied list with the mapping from internal RCT item indices to
//  * external indices.
//  *
//  * @param result An array to store the mapping in.
//  * @param capacity The capacity of the target array.
//  *
//  * @return If successful, the number of RCT items is returned. Otherwise,
//  zero
//  *         is returned.
//  */
// _RCT_TMPL_DECL
// int RCT::getInternToExternMapping(int* result, int capacity) const {
//   int i;
//   if ((result == NULL) || (capacity < size)) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Result list capacity is too small.";
//     }
//     return 0;
//   }
//   for (i = 0; i < size; i++) {
//     result[i] = intern_to_extern_mapping_[i];
//   }
//   return size;
// }

// /*!
//  * Fills the supplied list with the RCT level set sizes, from smallest to
//  * largest. The result does not include the sample consisting solely of
//  * the virtual RCT root item.
//  *
//  * @param result An array to hold the level set sizes.
//  * @param capacity The capacity of the target array.
//  *
//  * @return If successful, the number of RCT levels is returned (excluding
//  that
//  *         of the root). If unsuccessful, zero is returned.
//  */
// _RCT_TMPL_DECL
// int RCT::getLevelSetSizes(int* result, int capacity) const {
//   int lvl;
//   if ((result == NULL) || (capacity < levels)) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Result list capacity is too small.";
//     }
//     return 0;
//   }
//   for (lvl = 0; lvl < levels; lvl++) {
//     result[lvl] = levelSetSizeList[lvl];
//   }
//   return levels;
// }

// /*!
//  * Retrieve the maximum node degree occurring in the RCT.
//  *
//  * @return The degree of the node with the largest number of children.
//  */
// _RCT_TMPL_DECL
// int RCT::getMaxDegree() const { return max_degree_; }

// /*!
//  * Fills the supplied list with the mapping from external item indices to
//  * internal RCT level sets.
//  *
//  * @param result An array to hold the item heights.
//  * @param capacity The capacity of the target array.
//  *
//  * @return If successful, the number of RCT items is returned. If
//  unsuccessful,
//  *         zero is returned.
//  */
// _RCT_TMPL_DECL
// int RCT::getMaxLevelAssignment(int* result, int capacity) const {
//   int i;
//   int lvl;
//   if ((result == NULL) || (capacity < size)) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Result list capacity is too small.";
//     }
//     return 0;
//   }
//   for (lvl = 0; lvl < levels; lvl++) {
//     for (i = levelSetSizeList[lvl + 1]; i < levelSetSizeList[lvl]; i++) {
//       result[intern_to_extern_mapping_[i]] = lvl;
//     }
//   }
//   result[intern_to_extern_mapping_[0]] = levels;
//   return size;
// }

// /*!
//  * Retrieve the maximum number of parents a node in the RCT is allowed to
//  * have.
//  *
//  * @return The maxmimum number of allowed parents per node.
//  */
// _RCT_TMPL_DECL
// int RCT::getMaxParents() const { return max_parents_; }

// /*!
//  * Retrieve the number of data items stored in the RCT.
//  *
//  * @return The number of data items the RCT is built on.
//  */
// _RCT_TMPL_DECL
// int RCT::getNumItems() const { return size; }

// /*!
//  * Retrieve the number of levels in the RCT.
//  *
//  * @note Since the RCT is built on level sets whose membership in items
//  *       is purely randomized one can observe variance in height when
//  *       rebuilding an RCT on the same data with a different random leveling
//  *       (using different seeds for example).
//  */
// _RCT_TMPL_DECL
// int RCT::getNumLevels() const { return levels; }

// /*!
//  * Retrieve the number of nodes in the RCT.
//  *
//  * @note Although the model of the RCT relies on using copies of nodes
//  *       in different level sets, the actual implementation does not do this.
//  *       The number of nodes should therefore be seen as a symbolic property
//  *       of the random leveling and not as a contributor to the search
//  *       or construction time complexity.
//  *
//  * @return The number of nodes in the RCt.
//  */
// _RCT_TMPL_DECL
// int RCT::getNumNodes() const { return num_nodes_; }

// /*!
//  * Computes the recall accuracy of the most recent query result. A list of
//  the
//  * exact distances must be provided, sorted from smallest to largest. The
//  * number of exact distances provided determines the size of the
//  neighbourhood
//  * within which the accuracy is assessed.
//  *
//  * @note The list must contain at least as many entries as the number of
//  items
//  *       found in the query result.
//  *
//  * @return If unsuccessful, a negative value is returned.
//  */
// _RCT_TMPL_DECL
// float RCT::getResultAcc(float* exactDistList, int how_many) const {
//   int i;
//   int loc = 0;
//   if ((exactDistList == NULL) || (how_many < query_result_size_)) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Exact distance list is too small.";
//     }
//     return kUnknown;
//   }
//   for (i = 0; i < how_many; i++) {
//     if ((loc < query_result_size_) &&
//         (queryResultDistList[loc] <= exactDistList[i])) {
//       loc++;
//     }
//   }
//   return ((float)loc) / how_many;
// }

// /*!
//  * Fills the supplied list with the query-to-neighbour distances found in the
//  * most recent RCT query.
//  *
//  * @param result An array to store the distance values.
//  * @param capacity The capacity of the target array.
//  *
//  * @return If successful, the number of items found is returned. Otherwise,
//  *         zero is returned.
//  */
// _RCT_TMPL_DECL
// int RCT::getResultDists(float* result, int capacity) const {
//   int i;
//   if ((result == NULL) || (capacity < query_result_size_)) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Result list capacity is too small.";
//     }
//     return 0;
//   }
//   for (i = 0; i < query_result_size_; i++) {
//     result[i] = queryResultDistList[i];
//   }
//   return query_result_size_;
// }

// /*!
//  * Retrieve the number of distance comparisons performed during the most
//  * recent operation.
//  *
//  * @return The number of distance comparisons.
//  */
// _RCT_TMPL_DECL
// unsigned long RCT::getResultDistComps() const { return numDistComps; }

// /*!
//  * Fills the supplied list with the (external) indices of the items found in
//  * the most recent RCT query.
//  *
//  * @param result An array to store the item indices.
//  * @param capacity The capacity of the target array.
//  *
//  * @return If successful, the number of items found is returned. Otherwise,
//  *         zero is returned.
//  */
// _RCT_TMPL_DECL
// int RCT::getResultIndices(int* result, int capacity) const {
//   int i;
//   if ((result == NULL) || (capacity < query_result_size_)) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Result list capacity is too small.";
//     }
//     return 0;
//   }
//   for (i = 0; i < query_result_size_; i++) {
//     result[i] = intern_to_extern_mapping_[queryResultIndexList[i]];
//   }
//   return query_result_size_;
// }

// /*!
//  * Returns the number of items found in the most recent query.
//  *
//  * @return Number of results found in the most recent query.
//  */
// _RCT_TMPL_DECL
// int RCT::getResultNumFound() const { return query_result_size_; }

// /*!
//  * Returns the sample size used in the most recent query.
//  *
//  * @return Sample size used in the most recent query.
//  */
// _RCT_TMPL_DECL
// int RCT::getResultSampleSize() const { return query_result_sample_size_; }

// /*!
//  * Retrieve the seed value used to initialize the random number generator.
//  */
// _RCT_TMPL_DECL
// unsigned long RCT::getRNGSeed() const { return seed; }

// /*!
//  * Resets the current query object to NULL. This has the effect of
//  * clearing any saved distances - subsequent <code>findNear</code> and
//  * <code>findNearest</code> operations would be forced to compute
//  * all needed distances from scratch.
//  */
// _RCT_TMPL_DECL
// void RCT::resetQuery() { setNewQuery(NULL); }

// /*!
//  * Save the RCT to the specified file. The extension ".rctf" is automatically
//  * appended to the file name.
//  *
//  * @param fileName The file name to save the RCT under.
//  * @return If successful, the number of RCT items is returned. Otherwise,
//  *         zero is returned.
//  */
// _RCT_TMPL_DECL
// int RCT::saveToFile(const char* fileName) const {
//   int i;
//   int j;
//   int lvl;
//   int num_children = 0;
//   int* childList = NULL;
//   ofstream outFile;
//   int levelSetSize = 0;

//   // If the RCT has not yet been built, abort.
//   if (size <= 0) {
//     return 0;
//   }

//   // Open the file for writing.
//   // If this fails, then abort.
//   if (fileName == NULL) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "Output file name is NULL.";
//     }
//     return 0;
//   }

//   // Attach extension '.rctf'.
//   ostringstream fullFileName;
//   fullFileName << fileName << ".rctf";

//   // Try to open the output file.
//   outFile.open(fullFileName.str().c_str(), ios::out);
//   if (!outFile.is_open()) {
//     if (verbosity > 0) {
//       LOG(ERROR) << "File '" << fullFileName.str() << "' could not be
//       opened.";
//     }
//     return 0;
//   }

//   // Begin writing to the output file. First, write a comment identifying the
//   // RCT version and the output file name.
//   outFile << "%% RCT " << kVersion << ' ' << fileName << endl;

//   // Write the main RCT parameters.
//   outFile << "%% size levels num_nodes_ max_parents_ max_degree_ avgDegree ";
//   outFile << "coverageParameter buildScaleFactor seed" << endl;
//   outFile << size << ' ' << levels << ' ' << num_nodes_ << ' ';
//   outFile << max_parents_ << ' ' << max_degree_ << ' ' << avgDegree;
//   outFile << ' ' << coverageParameter << ' ' << buildScaleFactor;
//   outFile << ' ' << seed << endl;

//   // Write the level sizes.
//   outFile << "%% level set sizes:" << endl;
//   for (lvl = 0; lvl < levels; lvl++) {
//     outFile << levelSetSizeList[lvl] << ' ';
//   }
//   outFile << 1 << endl;

//   // For each item at each level, write out:
//   //   its level,
//   //   its index,
//   //   the index of the item in the original input list,
//   //   the number of children of the item,
//   //   and a list of the indices of the children.
//   outFile << "%% level nodeID origItemID num_children c_0 c_1 ..." << endl;
//   for (lvl = levels; lvl >= 0; lvl--) {
//     levelSetSize = levelSetSizeList[lvl];
//     for (i = 0; i < levelSetSize; i++) {
//       num_children = childLSizeLList[lvl][i];
//       childList = childIndexLLList[lvl][i];
//       outFile << lvl << ' ' << i << ' ' << intern_to_extern_mapping_[i];
//       outFile << ' ' << num_children;
//       for (j = 0; j < num_children; j++) {
//         outFile << ' ' << childList[j];
//       }
//       outFile << endl;
//     }
//   }
//   outFile.close();
//   return size;
// }

// /*!
//  * Sets the coverage parameter to a specific value.
//  *
//  * @param coverageParameter The desired value of the coverage parameter.
//  * @return If the coverage parameter was changes true is returned.
//  *         Otherwise, false is returned.
//  */
// _RCT_TMPL_DECL
// bool RCT::setCoverageParameter(const float& coverageParameter) {
//   if ((*this).coverageParameter == coverageParameter) {
//     return false;
//   }
//   (*this).coverageParameter = coverageParameter;
//   return true;
// }

// /*!
//  * Sets the verbosity level for messages. Verbosity of zero or less:
//  * no messages produced. Verbosity of 1: error messages only.
//  * Verbosity of 2: error and progress messages only. Verbosity of
//  * 3 or more: error, progress, and debug messages reported.
//  */
// _RCT_TMPL_DECL
// void RCT::setVerbosity(const int& verbosity) {
//   if (verbosity <= 0) {
//     (*this).avgDegree = 0;
//   } else if (verbosity >= 3) {
//     (*this).verbosity = 3;
//   } else {
//     (*this).verbosity = verbosity;
//   }
// }

// /*!
//  * Returns the distance from the current query object to the specified data
//  * object. If the distance has already been computed and stored, the stored
//  * distance is returned. Otherwise, the distance is computed and stored
//  * before returning it.
//  *
//  * @param itemIndex The internal index of a data item.
//  * @return The distance from the current query to that item.
//  */
// _RCT_TMPL_DECL
// float RCT::computeDistFromQuery(int itemIndex) {
//   if (distFromQueryList[itemIndex] == kUnknown) {
//     distFromQueryList[itemIndex] =
//         query->distanceTo(data_[intern_to_extern_mapping_[itemIndex]]);
//     storedDistIndexList[num_stored_dists_] = itemIndex;
//     num_stored_dists_++;
//     numDistComps++;
//   }
//   return distFromQueryList[itemIndex];
// }

// /*!
//  * Builds an RCT on items in the first locations of the scrambled data array.
//  */

// /*!
//  * Perform an approximate nearest-neighbor query for the specified item. The
//  * number of desired nearest neighbors <code>how_many</code> must be
//  specified.
//  * The search is relative to the random level sample_level.
//  *
//  * @param coverage Additional search efforts.
//  * @param how_many The desired number of neighbors.
//  * @param sample_level The sample level with respect to which the query is
//  *                    performed.
//  *
//  * @return The number of elements actually found.
//  */
// _RCT_TMPL_DECL
// int RCT::doFindNear(int how_many, int sample_level, float coverage) {
//   int i;
//   int j;
//   int lvl;
//   int child = 0;
//   int node_index = 0;
//   int num_children = 0;
//   // int tempQueryResultSize = 0;
//   double var_quota = 0.0F;
//   int num_found = 0;
//   int num_retained = 0;
//   int* childList = NULL;

//   // Compute quota of items to be retained at every level.
//   // Rank cover tree rules.levelQuotaList
//   var_quota = (double)how_many;
//   for (lvl = sample_level; lvl < levels; lvl++) {
//     levelQuotaList[lvl] =
//         (int)((coverage * var_quota * coverageParameter) + 0.999999F);
//     if (levelQuotaList[lvl] < coverage * coverageParameter) {
//       levelQuotaList[lvl] =
//           (int)((coverage * coverageParameter) + 0.999999F);
//     }
//     var_quota /= sample_rate_;
//   }
//   if (how_many > levelQuotaList[sample_level]) {
//     levelQuotaList[sample_level] = how_many;
//   }

//   // Load the root as the tentative sole member of the query result list.
//   query_result_size_ = 0;
//   queryResultDistList[0] = computeDistFromQuery(0);
//   queryResultIndexList[0] = 0;
//   num_retained = 1;

//   // From the root, search out other nodes to place in the query result.
//   for (lvl = levels - 1; lvl >= sample_level; lvl--) {
//     // For every node at the active level, load its children
//     //   into the scratch list, and compute their distances to the query.
//     num_found = 0;

//     for (i = 0; i < num_retained; i++) {
//       node_index = queryResultIndexList[i];
//       num_children = childLSizeLList[lvl + 1][node_index];
//       childList = childIndexLLList[lvl + 1][node_index];

//       for (j = 0; j < num_children; j++) {
//         child = childList[j];

//         if (visitedNodeIndexList[child] != TRUE) {
//           visitedNodeIndexList[child] = TRUE;
//           tempResultIndexList[num_found] = child;
//           tempResultDistList[num_found] = computeDistFromQuery(child);
//           num_found++;
//         }
//       }
//     }

//     for (i = 0; i < num_found; i++) {
//       visitedNodeIndexList[tempResultIndexList[i]] = FALSE;
//     }

//     // Extract the closest nodes from the list of accumulated children,
//     //   and keep them as the tentative parents of the query.

//     if (num_found > levelQuotaList[lvl]) {
//       num_retained = levelQuotaList[lvl];
//     } else {
//       num_retained = num_found;
//     }

//     num_retained = PartialQuickSort(num_retained, tempResultDistList,
//                                    tempResultIndexList, 0, num_found - 1);

//     for (i = 0; i < num_retained; i++) {
//       queryResultIndexList[i] = tempResultIndexList[i];
//       queryResultDistList[i] = tempResultDistList[i];
//     }
//   }

//   // Select the final number of neighbors needed.
//   if (num_retained > how_many) {
//     query_result_size_ = how_many;
//   } else {
//     query_result_size_ = num_retained;
//   }
//   childList = NULL;
//   return query_result_size_;
// }

// /*!
//  * Perform an exact nearest-neighbor query for the specified item. The
//  * number of desired nearest neighbors <code>how_many</code> must be
//  specified.
//  * The search is relative to the random level sample_level.
//  *
//  * @param how_many The desired number of neighbors.
//  * @param sample_level The sample level with respect to which the query is
//  *                    performed.
//  *
//  * @return The number of elements actually found.
//  */
// _RCT_TMPL_DECL
// int RCT::doFindNearest(int how_many, int sample_level) {
//   int i;

//   // Handle the singleton case separately.
//   if (size == 1) {
//     query_result_size_ = 1;
//     queryResultDistList[0] = computeDistFromQuery(0);
//     queryResultIndexList[0] = 0;

//     return 1;
//   }

//   query_result_size_ = levelSetSizeList[sample_level];

//   // Compute distances from the current query to all items.
//   for (i = 0; i < query_result_size_; i++) {
//     queryResultDistList[i] = computeDistFromQuery(i);
//     queryResultIndexList[i] = i;
//   }

//   query_result_size_ =
//       PartialQuickSort(how_many, queryResultDistList, queryResultIndexList,
//       0,
//                        query_result_size_ - 1);

//   return query_result_size_;
// }

// /*!
//  * Sorts the items in the index and distance list in ascending order w.r.t.
//  * their distances.
//  *
//  * @param distList List of item distances.
//  * @param indexList List of item indices.
//  * @param rangeFirst First element of range.
//  * @param rangeLast Last element of range.
//  */
// _RCT_TMPL_DECL
// void RCT::quickSort(float* distList, int* indexList, int rangeFirst,
//                     int rangeLast) {
//   int pivotLoc = 0;
//   float pivotDist;
//   int pivotIndex;
//   float tempDist;
//   int tempIndex;
//   int low = 0;
//   int high = 0;
//   float tieBreakDist;

//   // If the range is empty, or if it contains only one item,
//   //   then return immediately.

//   if (rangeLast <= rangeFirst) {
//     return;
//   }

//   // If the range to be sorted is small, just do an insertion sort.

//   if (rangeLast - rangeFirst < 7) {
//     std::uniform_int_distribution<size_t> distribution(0,
//                                                        rangeLast -
//                                                        rangeFirst);
//     high = rangeFirst + 1;
//     tieBreakDist = distList[distribution(rand_)];

//     // The outer while loop considers each item in turn (starting
//     //   with the second item in the range), for insertion into
//     //   the sorted list of items that precedes it.

//     while (high <= rangeLast) {
//       // Copy the next item to be inserted, as the "pivot".
//       // Start the insertion tests with its immediate predecessor.

//       pivotDist = distList[high];
//       pivotIndex = indexList[high];
//       low = high - 1;

//       // Work our way down through previously-sorted items
//       //   towards the start of the range.

//       while (low >= rangeFirst) {
//         // Compare the item to be inserted (the "pivot") with
//         //   the current item.

//         if (distList[low] < pivotDist) {
//           // The current item precedes the pivot in the sorted order.
//           // Break out of the loop - we have found the insertion point.

//           break;
//         } else if (distList[low] > pivotDist) {
//           // The current item follows the pivot in the sorted order.
//           // Shift the current item one spot upwards, to make room
//           //   for inserting the pivot below it.

//           distList[low + 1] = distList[low];
//           indexList[low + 1] = indexList[low];
//           low--;
//         } else {
//           // Break the tie pseudo-randomly.

//           if (((tieBreakDist < pivotDist) && (tieBreakDist < distList[low])
//           &&
//                (distList[low] < pivotDist)) ||
//               ((tieBreakDist >= pivotDist) && ((tieBreakDist < distList[low])
//               ||
//                                                (distList[low] < pivotDist))))
//                                                {
//             // The current item precedes the pivot in the sorted order.
//             // Break out of the loop - we have found the insertion point.

//             break;
//           } else {
//             // The current item follows the pivot in the sorted order.
//             // Shift the current item one spot upwards, to make room
//             //   for inserting the pivot below it.

//             distList[low + 1] = distList[low];
//             indexList[low + 1] = indexList[low];
//             low--;
//           }
//         }
//       }

//       // If we've made it to here, we've found the insertion
//       //   spot for the current element.
//       // Perform the insertion.

//       low++;
//       distList[low] = pivotDist;
//       indexList[low] = pivotIndex;

//       // Move to the next item to be inserted in the growing sorted list.

//       high++;
//     }

//     return;
//   }

//   // The range to be sorted is large, so do a quicksort.
//   // Select a pivot item, and swap it with the item at the beginning
//   //   of the range.

//   std::uniform_int_distribution<size_t> distribution(0, rangeLast -
//   rangeFirst);
//   pivotLoc = rangeFirst + distribution(rand_);
//   tieBreakDist = distList[distribution(rand_)];

//   pivotDist = distList[pivotLoc];
//   distList[pivotLoc] = distList[rangeFirst];
//   distList[rangeFirst] = pivotDist;

//   pivotIndex = indexList[pivotLoc];
//   indexList[pivotLoc] = indexList[rangeFirst];
//   indexList[rangeFirst] = pivotIndex;

//   // Partition the items with respect to the pivot.
//   // This efficient method is adapted from the one outlined in
//   //   Cormen, Leiserson & Rivest.
//   // The range is scanned from both ends.
//   // Items with small distances are placed below "low", and those
//   //   with large distances are placed above "high".
//   // Where "low" and "high" meet, the pivot item is inserted.

//   low = rangeFirst;
//   high = rangeLast + 1;

//   while (TRUE) {
//     // Move the "high" endpoint down until it meets either the pivot,
//     //   or something that belongs on the "low" side.
//     // If the key values are tied, decide pseudo-randomly.

//     do {
//       high--;
//     } while (
//         (distList[high] > pivotDist) ||
//         ((distList[high] == pivotDist) &&
//          (((tieBreakDist >= pivotDist) && (pivotDist < distList[high]) &&
//            (distList[high] <= tieBreakDist)) ||
//           ((tieBreakDist < pivotDist) && ((pivotDist < distList[high]) ||
//                                           (distList[high] <=
//                                           tieBreakDist))))));

//     // Move the "low" endpoint up until it meets either the pivot,
//     //   or something that belongs on the "high" side.
//     // If the key values are tied, decide pseudo-randomly.

//     do {
//       low++;
//     } while (
//         (low < high) &&
//         ((distList[low] < pivotDist) ||
//          ((distList[low] == pivotDist) &&
//           (((tieBreakDist < pivotDist) && (tieBreakDist < distList[low]) &&
//             (distList[low] < pivotDist)) ||
//            ((tieBreakDist >= pivotDist) && ((tieBreakDist < distList[low]) ||
//                                             (distList[low] <
//                                             pivotDist)))))));

//     // Have the "low" and "high" endpoints crossed?
//     // If not, we still have more work to do.

//     if (low < high) {
//       // Swap the misplaced items, and try again.

//       tempDist = distList[low];
//       distList[low] = distList[high];
//       distList[high] = tempDist;

//       tempIndex = indexList[low];
//       indexList[low] = indexList[high];
//       indexList[high] = tempIndex;
//     } else {
//       // We found the cross-over point.

//       break;
//     }
//   }

//   // The pivot value ends up at the location referenced by "high".
//   // Swap it with the pivot (which resides at the beginning of the range).

//   distList[rangeFirst] = distList[high];
//   distList[high] = pivotDist;

//   indexList[rangeFirst] = indexList[high];
//   indexList[high] = pivotIndex;

//   // The partition is complete.
//   // Recursively sort the remaining items.

//   quickSort(distList, indexList, rangeFirst, high - 1);
//   quickSort(distList, indexList, high + 1, rangeLast);
// }

// /*!
//  * Accepts a new item as the query object for future distance comparisons.
//  * Any previously-stored distances are cleared by this operation, except in
//  the
//  * case where the previous query object is identical to the current query
//  * object.
//  *
//  * @param query The new query item.
//  */
// _RCT_TMPL_DECL
// void RCT::setNewQuery(DistData* query) {
//   int i;
//   if (query != this->query) {
//     for (i = 0; i < num_stored_dists_; i++) {
//       distFromQueryList[storedDistIndexList[i]] = kUnknown;
//     }
//     this->query = query;
//     num_stored_dists_ = 0;
//   }
// }

// /*!
//  * Get the fraction of edges in the RCT that are well-formed. An edge from a
//  * parent p on level i to a child c is consider
//  * well-formed if and only if there is no item p' on level i,
//  * such that dist(p',c) is strictly less than dist(p,c).
//  *
//  * @note The RCT is well-formed if and only if all of its edges are
//  well-formed.
//  *       Edges between the artificial root node and its children are by
//  *       definition always well-formed and are not checked.
//  *
//  * @note Be aware that checking well-formedness is very expensive due to the
//  *       required exact nearest-neighbor queries!
//  *
//  * @return Fraction of edges that are well-formed.
//  */
// _RCT_TMPL_DECL
// double RCT::getFractionOfWellformedEdges() {
//   // Count the well-formed edges.
//   long long wellFormedEdges = 0LL;
//   long long edgesChecked = 0LL;

//   // Check all levels:
//   for (int L = levels - 1; L > 0; --L) {
//     // Check all nodes on level L.
//     for (int i = 0; i < levelSetSizeList[L]; ++i) {
//       // Select an item on level L as the parent.
//       DistData* parent = data_[intern_to_extern_mapping_[i]];

//       // Investigate each child of that parent.
//       int num_children = childLSizeLList[L][i];
//       int* children = childIndexLLList[L][i];
//       for (int j = 0; j < num_children; ++j) {
//         // Retrieve the child item.
//         DistData* child = data_[intern_to_extern_mapping_[children[j]]];

//         // Is this a copy of the parent?
//         if (parent == child) {
//           ++wellFormedEdges;
//           ++edgesChecked;
//           continue;
//         }

//         // Determine the actual parent-child distance.
//         float actualParentChildDistance = parent->distanceTo(child);

//         // Find the nearest-neighbor's distance of the child at the parent
//         level
//         // L.
//         setNewQuery(child);
//         doFindNearest(1, L);
//         float correctParentChildDistance;
//         getResultDists(&correctParentChildDistance, 1);

//         // Do the distances match?
//         if (actualParentChildDistance == correctParentChildDistance) {
//           ++wellFormedEdges;
//         }
//         ++edgesChecked;
//       }
//     }
//   }

//   // How many edges do we have in the RCT?
//   return static_cast<double>(wellFormedEdges) /
//          static_cast<double>(edgesChecked);
// }

#undef _RCT_TEMPLATE_DECL
#undef _RCT_TEMPLATE_DEFN
