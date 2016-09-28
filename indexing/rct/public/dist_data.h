#ifndef INDEXING_RCT_DIST_DATA_H_
#define INDEXING_RCT_DIST_DATA_H_

namespace indexing {
namespace rct {

class DistData {
 public:
  virtual ~DistData() {}
  virtual float distanceTo(DistData*) = 0;
};

}  // namespace rct
}  // namespace indexing

#endif  // INDEXING_RCT_DIST_DATA_H_
