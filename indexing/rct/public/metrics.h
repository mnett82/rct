#ifndef INDEXING_RCT_CONTRIB_METRICS_H_
#define INDEXING_RCT_CONTRIB_METRICS_H_

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "indexing/rct/public/dense_feature_vector.h"

namespace indexing {
namespace rct {
namespace detail {

// Type trait expressing whether a type 'T' is closed under addition.
template <typename T>
using is_closed_under_addition =
    ::std::is_same<T, decltype(::std::declval<T>() + ::std::declval<T>())>;

// Type trait expressing whether a type 'T' is closed under subtraction.
template <typename T>
using is_closed_under_subtraction =
    ::std::is_same<T, decltype(::std::declval<T>() - ::std::declval<T>())>;

// Type trait expressing whether a type 'T' is closed under multiplication.
template <typename T>
using is_closed_under_multiplication =
    ::std::is_same<T, decltype(::std::declval<T>() * ::std::declval<T>())>;

// Type trait expressing whether a type 'T' is a complex value.
template <typename T>
struct is_complex : std::false_type {};
template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

}  // namespace detail

// Implements squared Euclidean distance for dense feature vectors.
//
// When the arguments are dense feature vectors of distinct lengths, the
// computation substitutes the missing components of the shorter vector by
// default-constructed values of type 'T'.
template <typename T, typename Enabled = std::enable_if<
                          detail::is_closed_under_addition<T>::value &&
                          detail::is_closed_under_subtraction<T>::value &&
                          detail::is_closed_under_multiplication<T>::value>>
T SquaredEuclideanDistance(const DenseFeatureVector<T>& lhs,
                           const DenseFeatureVector<T>& rhs) {
  const size_t num_shared_attrs = std::min(lhs.size(), rhs.size());
  T accum = T();
  for (size_t i = 0; i < num_shared_attrs; ++i) {
    const T diff = lhs[i] - rhs[i];
    accum += diff * diff;
  }
  for (size_t i = num_shared_attrs; i < lhs.size(); ++i) {
    const T diff = lhs[i] - T();
    accum += diff * diff;
  }
  for (size_t i = num_shared_attrs; i < rhs.size(); ++i) {
    const T diff = T() - rhs[i];
    accum += diff * diff;
  }
  return accum;
}

// Implements Euclidean distance for real-valued dense feature vectors.
//
// When the arguments are dense feature vectors of distinct lengths, the
// computation substitutes the missing components of the shorter vector by
// default-constructed values of type 'T'.
template <typename T,
          typename Enabled = std::enable_if<std::is_floating_point<T>::value>>
T EuclideanDistance(const DenseFeatureVector<T>& lhs,
                    const DenseFeatureVector<T>& rhs) {
  return std::sqrt(SquaredEuclideanDistance<T>(lhs, rhs));
}

// Implements Euclidean distance for complex-valued dense feature vectors.
//
// When the arguments are dense feature vectors of distinct lengths, the
// computation substitutes the missing components of the shorter vector by
// default-constructed values of type 'T'.
template <typename T>
std::complex<T> EuclideanDistance(
    const DenseFeatureVector<std::complex<T>>& lhs,
    const DenseFeatureVector<std::complex<T>>& rhs) {
  return std::sqrt<std::complex<T>>(SquaredEuclideanDistance(lhs, rhs));
}

}  // namespace rct
}  // namespace indexing

#endif  // INDEXING_RCT_CONTRIB_METRICS_H_
