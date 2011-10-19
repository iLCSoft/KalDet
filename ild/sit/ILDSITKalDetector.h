#ifndef __ILDSITKALDETECTOR__
#define __ILDSITKALDETECTOR__

#include "kaltest/TVKalDetector.h"

#include "TMath.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDSITKalDetector : public TVKalDetector {
  
public:
  
  ILDSITKalDetector( const gear::GearMgr& gearMgr );
  
  ~ILDSITKalDetector();
  
private:
  
  void setupGearGeom( const gear::GearMgr& gearMgr ) ;
  
  int _nLayers ;
  double _bZ ;
  
  struct SIT_Layer {
    int nLadders;
    double phi0;
    double dphi;
    double senRMin;
    double supRMin;
    double length;
    double width;
    double offset;
    double senThickness;
    double supThickness;
  };
  std::vector<SIT_Layer> _SITgeo;
  
};



#endif
