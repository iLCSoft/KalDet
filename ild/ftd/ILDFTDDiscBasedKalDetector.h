#ifndef __ILDFTDDISCBASEDDETECTOR__
#define __ILDFTDDISCBASEDDETECTOR__

#include "kaltest/TVKalDetector.h"

class TNode;

namespace gear{
  class GearMgr ;
}


class ILDFTDDiscBasedKalDetector : public TVKalDetector {
public:

  /** Initialize the FTD from GEAR */
  ILDFTDDiscBasedKalDetector( const gear::GearMgr& gearMgr );
  
  ~ILDFTDDiscBasedKalDetector() {} ;
  
 private:

  void setupGearGeom( const gear::GearMgr& gearMgr ) ;

  int _nDisks ;
  double _bZ ;

  struct FTD_Disk {
    double rInner;
    double rOuter;
    double senThickness;
    double supThickness;
    double zPos;

  };
  std::vector<FTD_Disk> _FTDgeo;

};

#endif
