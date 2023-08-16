#ifndef __DISTDATA_H
#define __DISTDATA_H

class DistData {
 public:
  virtual ~DistData() = default;
  virtual float distanceTo(DistData*) = 0;
};

#endif /* __DISTDATA_H */
