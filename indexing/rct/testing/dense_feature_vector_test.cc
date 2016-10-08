#include <gtest/gtest.h>

#include "indexing/rct/public/dense_feature_vector.h"

namespace indexing {
namespace rct {
namespace testing {

// Note: This test is mostly here to make sure the template class is being
// instantiated; keep improving the tests as the feature vector object gains
// more complex logic.
TEST(DenseFeatureVectorTest, SanityCheck) {
  DenseFeatureVector<int> v({1, 2, 3});

  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 2);
  EXPECT_EQ(v[2], 3);
}

}  // namespace testing
}  // namespace rct
}  // namespace indexing
