#ifndef INDEXING_RCT_DIST_DATA_H_
#define INDEXING_RCT_DIST_DATA_H_

namespace indexing {
namespace rct {

// Defines an abstract distance metric.
//
// A distance metric on any type T implements a non-negative distance function 
// between pairs of Ts.
class DistData {
public:
  virtual ~DistData() {}
  virtual double distanceTo(const DistData&) = 0;
};

class MetricData {

}



class DistData {
 public:
  virtual ~DistData() {}
  virtual float distanceTo(DistData*) = 0;
};

}  // namespace rct
}  // namespace indexing

#endif  // INDEXING_RCT_DIST_DATA_H_
