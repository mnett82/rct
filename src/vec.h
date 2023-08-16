#ifndef DATA_H_
#define DATA_H_

#include <string_view>
#include <utility>
#include <vector>

#include "DistData.h"

// Represents a densly-packed real-valued vector.
//
// Euclidean distance is currently built in and should be made into a parameter.
class Vec final : public DistData {
 public:
  using Repr = std::vector<float>;

  ~Vec() override = default;

  explicit Vec(Repr repr) : repr_(std::move(repr)) {}

  float distanceTo(DistData* other) override;

 private:
  Repr repr_;
};

std::vector<Vec> LoadVecsFromHDF5(std::string_view hdf5_file,
                                  std::string_view dset_name);

#endif  // DATA_H_