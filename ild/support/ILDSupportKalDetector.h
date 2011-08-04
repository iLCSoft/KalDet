#ifndef __ILDSUPPORTDETECTOR__
#define __ILDSUPPORTDETECTOR__

#include "kaltest/TVKalDetector.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDSupportKalDetector : public TVKalDetector {
public:

  /** Initialize the support structures from GEAR */
  ILDSupportKalDetector( const gear::GearMgr& gearMgr );
  
  ~ILDSupportKalDetector();
  
 private:

};

#endif
