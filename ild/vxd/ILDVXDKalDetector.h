#ifndef __ILDVXDKALDETECTOR__
#define __ILDVXDKALDETECTOR__

/** Ladder based VXD to be used for ILD DBD studies 
 *
 * @author S.Aplin DESY
 */

#include "kaltest/TVKalDetector.h"

#include "TMath.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDVXDKalDetector : public TVKalDetector {
  
public:
  
  ILDVXDKalDetector( const gear::GearMgr& gearMgr );
  
  
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
  
  
   struct VXD_Cryostat {
     double alRadius;
     double alThickness;
     double alInnerR;
     double alZEndCap;
     double alHalfZ;
     double shellRadius;
     double shellThickness;
     double shellInnerR;
     double shellZEndCap;
     double shelllHalfZ;

     bool   exists;
   };

  VXD_Cryostat _vxd_Cryostat;
  
};



#endif
