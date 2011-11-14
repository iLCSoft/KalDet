#ifndef __ILDFTDDETECTOR__
#define __ILDFTDDETECTOR__

/** Petal based FTD to be used for ILD DBD studies 
 * WARNING: Still very experimental
 *
 * @author S.Aplin DESY
 */

#include "kaltest/TVKalDetector.h"

class TNode;
class TVector3;

namespace gear{
  class GearMgr ;
}


class ILDFTDKalDetector : public TVKalDetector {
public:
  
  /** Initialize the FTD from GEAR */
  ILDFTDKalDetector( const gear::GearMgr& gearMgr );
  
  
private:
  
  struct FTD_Petal {
    
    int    ipetal;
    double phi;
    double alpha;
    double rInner;
    double height;
    double innerBaseLength;
    double outerBaseLength;
    double senThickness;
    double supThickness;
    double senZPos;
    bool faces_ip;
    
  };
  
  
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
    
    double senZPos_even_petal1;
    double senZPos_even_petal2;
    double senZPos_even_petal3;
    double senZPos_even_petal4;
    
    double supZPos_even_petal1;
    double supZPos_even_petal2;
    double supZPos_even_petal3;
    double supZPos_even_petal4;
    
    double senZPos_odd_petal1;
    double senZPos_odd_petal2;
    double senZPos_odd_petal3;
    double senZPos_odd_petal4;
    
    double supZPos_odd_petal1;
    double supZPos_odd_petal2;
    double supZPos_odd_petal3;
    double supZPos_odd_petal4;
    
    
    
  };
  
  void build_turbine_design();
  void build_staggered_design();
  
  //void create_petal(TVector3 measurement_plane_centre, FTD_Petal petal, int CellID);
  void create_segmented_disk_layers(int idisk, int nsegments, bool even_petals, double phi0, double zpos, bool front);
  
  
  void setupGearGeom( const gear::GearMgr& gearMgr ) ;
  
  int _nDisks ;
  double _bZ ;
  
  bool _is_staggered_design;
  
  std::vector<FTD_Disk> _FTDgeo;
  
};

#endif
