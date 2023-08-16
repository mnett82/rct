#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "rct.h"
#include "vec.h"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " hdf5_file" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Vec> vecs = LoadVecsFromHDF5(argv[1], "train");

  auto rct = std::make_unique<RCT>();
  rct->setVerbosity(2);
  rct->setCoverageParameter(8.0);
  rct->setSampleRate(pow(vecs.size(), 1.0 / 3.0));

  std::vector<DistData*> ptrs;
  ptrs.reserve(vecs.size());
  for (Vec& vec : vecs) {
    ptrs.push_back(&vec);
  }

  rct->build(ptrs.data(), ptrs.size());

  rct->saveToFile(argv[1]);

  return 0;
}