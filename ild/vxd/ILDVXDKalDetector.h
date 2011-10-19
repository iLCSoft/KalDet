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
  
private:
  
  void setupGearGeom( const gear::GearMgr& gearMgr ) ;
  
  int _nLayers ;
  double _bZ ;
  
  struct VXD_Layer {
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
  std::vector<VXD_Layer> _VXDgeo;
  
};



#endif
