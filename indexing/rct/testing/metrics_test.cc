#include <gtest/gtest.h>
#include <algorithm>
#include <functional>
#include <random>
#include <vector>

#include "indexing/rct/public/dense_feature_vector.h"
#include "indexing/rct/public/metrics.h"

namespace indexing {
namespace rct {
namespace testing {
namespace {

constexpr int kMaxLength = 1000;
constexpr int kMaxAttemptsPerTest = 1000;

}  // namespace

class MetricsTest : public ::testing::Test {
 protected:
  MetricsTest() : rand_(::testing::UnitTest::GetInstance()->random_seed()) {}
  ~MetricsTest() override {}

  std::vector<double> RandomRealVector() {
    std::vector<double> reals(
        std::uniform_int_distribution<>(0, kMaxLength)(rand_));
    std::generate(
        std::begin(reals), std::end(reals),
        std::bind(std::normal_distribution<double>(0.0, 10.0), rand_));
    return reals;
  }

 private:
  std::mt19937 rand_;
};

TEST_F(MetricsTest, SquaredEuclideanDistanceIdentifiesIndiscernables) {
  for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
    const DenseFeatureVector<double> v(RandomRealVector());

    EXPECT_DOUBLE_EQ(0.0, SquaredEuclideanDistance(v, v));
  }
}

TEST_F(MetricsTest, SquaredEuclideanDistanceIsSymmetric) {
  for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
    const DenseFeatureVector<double> v(RandomRealVector());
    const DenseFeatureVector<double> w(RandomRealVector());

    EXPECT_DOUBLE_EQ(SquaredEuclideanDistance(v, w),
                     SquaredEuclideanDistance(w, v));
  }
}

TEST_F(MetricsTest, SquaredEuclideanDistanceIsNonNegative) {
  for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
    const DenseFeatureVector<double> v(RandomRealVector());
    const DenseFeatureVector<double> w(RandomRealVector());

    EXPECT_LE(0.0, SquaredEuclideanDistance(v, w));
  }
}

TEST_F(MetricsTest, EuclideanDistanceIdentifiesIndiscernables) {
  for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
    const DenseFeatureVector<double> v(RandomRealVector());

    EXPECT_DOUBLE_EQ(0.0, EuclideanDistance(v, v));
  }
}

TEST_F(MetricsTest, EuclideanDistanceIsSymmetric) {
  for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
    const DenseFeatureVector<double> v(RandomRealVector());
    const DenseFeatureVector<double> w(RandomRealVector());

    EXPECT_DOUBLE_EQ(EuclideanDistance(v, w), EuclideanDistance(w, v));
  }
}

TEST_F(MetricsTest, EuclideanDistanceIsNonNegative) {
  for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
    const DenseFeatureVector<double> v(RandomRealVector());
    const DenseFeatureVector<double> w(RandomRealVector());

    EXPECT_LE(0.0, EuclideanDistance(v, w));
  }
}

TEST_F(MetricsTest, EuclideanDistanceIsSubadditive) {
  for (int i = 0; i < kMaxAttemptsPerTest; ++i) {
    const DenseFeatureVector<double> u(RandomRealVector());
    const DenseFeatureVector<double> v(RandomRealVector());
    const DenseFeatureVector<double> w(RandomRealVector());

    EXPECT_LE(EuclideanDistance(u, w),
              EuclideanDistance(u, v) + EuclideanDistance(v, w));
  }
}

}  // namespace testing
}  // namespace rct
}  // namespace indexing
