#ifndef __ILDVXDKALDETECTOR__
#define __ILDVXDKALDETECTOR__

#include "kaltest/TVKalDetector.h"

#include "TMath.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDVXDKalDetector : public TVKalDetector {
    
 public:

  ILDVXDKalDetector( const gear::GearMgr& gearMgr );

  ~ILDVXDKalDetector();
  
};



#endif
