#ifndef __ILDTPCDETECTOR__
#define __ILDTPCDETECTOR__

#include "kaltest/TVKalDetector.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDTPCKalDetector : public TVKalDetector {
public:
  
  /** Initialize the TPC from GEAR */
  ILDTPCKalDetector( const gear::GearMgr& gearMgr );
  
  ~ILDTPCKalDetector();
  
private:
  
};

#endif
