#ifndef INDEXING_RCT_PUBLIC_DENSE_FEATURE_VECTOR_H_
#define INDEXING_RCT_PUBLIC_DENSE_FEATURE_VECTOR_H_

#include <cstddef>
#include <utility>
#include <vector>

namespace indexing {
namespace rct {

// Implements a dense feature vector.
//
// TODO: Extract attribute storage type as template parameter.
//
// This class is thread-safe.
template <typename T>
class DenseFeatureVector {
 public:
  virtual ~DenseFeatureVector() {}

  explicit DenseFeatureVector(std::vector<T> attributes)
      : attributes_(std::move(attributes)) {}

  // Return the number of attributes in the vector.
  size_t size() const { return attributes_.size(); }

  // Return the i-th attribute of the vector.
  typename std::vector<T>::const_reference operator[](const size_t i) const {
    return attributes_[i];
  }

 private:
  // Attributes of the dense feature vector.
  const std::vector<T> attributes_;

  DenseFeatureVector() = delete;
  DenseFeatureVector(const DenseFeatureVector&) = delete;
  DenseFeatureVector& operator=(const DenseFeatureVector&) = delete;
};

}  // namespace rct
}  // namespace indexing

#endif  // INDEXING_RCT_PUBLIC_DENSE_FEATURE_VECTOR_H_
