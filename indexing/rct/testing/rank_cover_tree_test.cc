#include "indexing/rct/public/rank_cover_tree.h"

#include <iterator>
#include <algorithm>
#include <functional>
#include <memory>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "indexing/rct/public/dense_feature_vector.h"
#include "indexing/rct/public/metrics.h"

using std::bind;
using std::generate;
using std::mt19937;
using std::normal_distribution;
using std::unique_ptr;
using std::vector;
using testing::UnitTest;
using std::transform;
using std::back_inserter;

namespace indexing {
namespace rct {
namespace testing {
namespace {

constexpr int kNumElements = 1000;
constexpr int kNumFeatures = 1000;

}  // namespace

class RankCoverTreeTest : public ::testing::Test {
 protected:
  RankCoverTreeTest()
      : data_(kNumElements), rand_(UnitTest::GetInstance()->random_seed()) {}
  ~RankCoverTreeTest() override {}

  vector<const DenseFeatureVector<double>*> RandomDataSet() {
    for (size_t i = 0; i < data_.size(); ++i) {
      data_[i].reset(new DenseFeatureVector<double>(RandomRealVector()));
    }
    vector<const DenseFeatureVector<double>*> unowned_ptrs;
    transform(data_.begin(), data_.end(), back_inserter(unowned_ptrs),
              [](const unique_ptr<DenseFeatureVector<double>>& element) {
                return element.get();
              });
    return unowned_ptrs;
  }

 private:
  vector<double> RandomRealVector() {
    vector<double> reals(kNumFeatures);
    generate(reals.begin(), reals.end(),
             bind(normal_distribution<double>(0.0, 10.0), rand_));
    return reals;
  }

  vector<unique_ptr<DenseFeatureVector<double>>> data_;
  mt19937 rand_;
};

TEST_F(RankCoverTreeTest, Build) {
  for (int i = 0; i < 5; ++i) {
    RankCoverTree<DenseFeatureVector<double>, double, EuclideanDistance> rct;

    EXPECT_TRUE(rct.Build(RandomDataSet(), 10.0, 1, 0.5));
  }
}

// TEST_F(MetricsTest, SquaredEuclideanDistanceIdentifiesIndiscernables) {
//   for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
//     const DenseFeatureVector<double> v(RandomRealVector());

//     EXPECT_DOUBLE_EQ(0.0, SquaredEuclideanDistance(v, v));
//   }
// }

// TEST_F(MetricsTest, SquaredEuclideanDistanceIsSymmetric) {
//   for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
//     const DenseFeatureVector<double> v(RandomRealVector());
//     const DenseFeatureVector<double> w(RandomRealVector());

//     EXPECT_DOUBLE_EQ(SquaredEuclideanDistance(v, w),
//                      SquaredEuclideanDistance(w, v));
//   }
// }

// TEST_F(MetricsTest, SquaredEuclideanDistanceIsNonNegative) {
//   for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
//     const DenseFeatureVector<double> v(RandomRealVector());
//     const DenseFeatureVector<double> w(RandomRealVector());

//     EXPECT_LE(0.0, SquaredEuclideanDistance(v, w));
//   }
// }

// TEST_F(MetricsTest, EuclideanDistanceIdentifiesIndiscernables) {
//   for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
//     const DenseFeatureVector<double> v(RandomRealVector());

//     EXPECT_DOUBLE_EQ(0.0, EuclideanDistance(v, v));
//   }
// }

// TEST_F(MetricsTest, EuclideanDistanceIsSymmetric) {
//   for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
//     const DenseFeatureVector<double> v(RandomRealVector());
//     const DenseFeatureVector<double> w(RandomRealVector());

//     EXPECT_DOUBLE_EQ(EuclideanDistance(v, w), EuclideanDistance(w, v));
//   }
// }

// TEST_F(MetricsTest, EuclideanDistanceIsNonNegative) {
//   for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
//     const DenseFeatureVector<double> v(RandomRealVector());
//     const DenseFeatureVector<double> w(RandomRealVector());

//     EXPECT_LE(0.0, EuclideanDistance(v, w));
//   }
// }

// TEST_F(MetricsTest, EuclideanDistanceIsSubadditive) {
//   for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
//     const DenseFeatureVector<double> u(RandomRealVector());
//     const DenseFeatureVector<double> v(RandomRealVector());
//     const DenseFeatureVector<double> w(RandomRealVector());

//     EXPECT_LE(EuclideanDistance(u, w),
//               EuclideanDistance(u, v) + EuclideanDistance(v, w));
//   }
// }

}  // namespace testing
}  // namespace rct
}  // namespace indexing
