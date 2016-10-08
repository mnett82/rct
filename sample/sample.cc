#include <cstddef>
#include <random>
#include <vector>

#include <glog/logging.h>

#include "indexing/rct/public/dense_vec_data.h"

using indexing::rct::DenseVecData;

std::vector<DenseVecData<float>*> CreateRandomData() {
  constexpr size_t kNumPoints = 1000;
  constexpr size_t kDimension = 100;

  std::vector<DenseVecData<float>*> data;
  for (size_t i = 0; i < kNumPoints; ++i) {
    std::mt19937 random;
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);

    auto* vec_data = new DenseVecData<float>(kDimension);
    for (size_t j = 0; j < kDimension; ++j) {
      (*vec_data)[j] = unit(random);
    }
    data.push_back(vec_data);
  }
  return data;
}

int main(int argc, char** argv) {
  // TODO: See why command line flags are not parsed properly.
  FLAGS_logtostderr = 1;
  FLAGS_v = 3;
  google::InitGoogleLogging(argv[0]);

  VLOG(1) << "Generating random data...";
  std::vector<DenseVecData<float>*> data = CreateRandomData();
  
  return 0;
}
