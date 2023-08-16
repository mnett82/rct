#include "vec.h"

#include <H5Cpp.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "DistData.h"

using namespace H5;

float Vec::distanceTo(DistData* const other) {
  const Repr& other_repr = static_cast<Vec*>(other)->repr_;

  float s = 0.0f;
  for (std::size_t i = 0; i < repr_.size(); ++i) {
    s += std::pow(repr_[i] - other_repr[i], 2.0f);
  }
  return std::sqrt(s);
}

std::vector<Vec> LoadVecsFromHDF5(std::string_view hdf5_file,
                                  std::string_view dset_name) {
  std::cout << "Opening file '" << hdf5_file << "'...\n";
  H5File file(std::string(hdf5_file), H5F_ACC_RDONLY);

  std::cout << "Loading data set '" << dset_name << "'...\n";
  DataSet dset = file.openDataSet(std::string(dset_name));
  assert(dset.getTypeClass() == H5T_FLOAT);

  DataSpace dspace = dset.getSpace();
  hsize_t dims[2];
  hsize_t rank = dspace.getSimpleExtentDims(dims, NULL);
  assert(rank == 2);
  assert(dims[0] > 0);
  assert(dims[1] > 0);
  std::cout << "Data set has rank " << rank << " and dimensions " << dims[0]
            << "," << dims[1] << "\n";

  hsize_t mdims[1] = {dims[1]};
  DataSpace mspace(1, mdims);

  hsize_t dcount[2] = {1, dims[1]};
  hsize_t doffset[2] = {0, 0};
  hsize_t mcount[1] = {dims[1]};
  hsize_t moffset[1] = {0};
  std::vector<Vec> vecs;
  vecs.reserve(dims[0]);
  mspace.selectHyperslab(H5S_SELECT_SET, mcount, moffset);
  for (hsize_t i = 0; i < dims[0]; ++i) {
    std::vector<float> values(dims[1]);
    doffset[0] = i;
    dspace.selectHyperslab(H5S_SELECT_SET, dcount, doffset);
    dset.read(values.data(), PredType::NATIVE_FLOAT, mspace, dspace);
    vecs.emplace_back(Vec(std::move(values)));
  }

  file.close();
  return vecs;
}