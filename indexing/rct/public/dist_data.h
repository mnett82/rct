#ifndef INDEXING_RCT_DIST_DATA_H_
#define INDEXING_RCT_DIST_DATA_H_

class DistData {
 public:
  virtual ~DistData() {}
  virtual float distanceTo(DistData*) = 0;
};

#endif  // INDEXING_RCT_DIST_DATA_H_
