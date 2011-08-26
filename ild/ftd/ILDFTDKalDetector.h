#ifndef __ILDFTDDETECTOR__
#define __ILDFTDDETECTOR__

#include "kaltest/TVKalDetector.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDFTDKalDetector : public TVKalDetector {
public:

  /** Initialize the FTD from GEAR */
  ILDFTDKalDetector( const gear::GearMgr& gearMgr );
  
  ~ILDFTDKalDetector() {} ;
  
 private:

  void setupGearGeom( const gear::GearMgr& gearMgr ) ;

  int _nDisks ;
  double _bZ ;

  struct FTD_Disk {
    int nPetals;
    double phi0;
    double dphi;
    double alpha;
    double rInner;
    double height;
    double innerBaseLength;
    double outerBaseLength;
    double senThickness;
    double supThickness;
    double zPos;

  };
  std::vector<FTD_Disk> _FTDgeo;

};

#endif
