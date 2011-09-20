#ifndef __ILDSUPPORTDETECTOR__
#define __ILDSUPPORTDETECTOR__

#include "kaltest/TVKalDetector.h"

class TNode;

namespace gear{
  class GearMgr ;
}

class ILDCylinderMeasLayer;

class ILDSupportKalDetector : public TVKalDetector {
public:

  /** Initialize the support structures from GEAR */
  ILDSupportKalDetector( const gear::GearMgr& gearMgr );
  
  ~ILDSupportKalDetector();

  ILDCylinderMeasLayer* getIPLayer() { return _ipLayer; }
	
 private:

	ILDCylinderMeasLayer* _ipLayer;
	
};

#endif
