#ifndef ILDMEASUREMENTSURFACESTOREFILLER_H
#define ILDMEASUREMENTSURFACESTOREFILLER_H

#include "gear/gearsurf/MeasurementSurfaceStore.h"

#include <vector>

namespace gear{
  class GearMgr;
  class ZPlanarParameters;
  class FTDParameters;
}

using namespace GearSurfaces;

class ILDMeasurementSurfaceStoreFiller : protected MeasurementSurfaceStoreFiller{
  
  public:
  
  ILDMeasurementSurfaceStoreFiller(gear::GearMgr* gear_mgr) : _gear_mgr(gear_mgr) {

    get_gear_parameters();
    
  }
  
  protected:
  
  void fill_store( std::vector<MeasurementSurface*>& surface_list ) const;
  
  private:
  
  /** adds MeasurementSufaces to the store
   * @param param: the ZPlanarParameters pointer of the detector, of which the measurement surfaces shall be added
   * 
   * @param det_id: the detector id (as in ILDConf)
   */
  void storeZPlanar( const gear::ZPlanarParameters* param , int det_id, std::vector<MeasurementSurface*>& surface_list ) const;
  
  void storeFTD( const gear::FTDParameters* param, std::vector<MeasurementSurface*>& surface_list ) const;
  
  gear::GearMgr* _gear_mgr;
  
  void get_gear_parameters();
  
#define HARDCODEDGEAR 1
#ifdef HARDCODEDGEAR
  
  
  /** the strip angles for every layer */
  std::vector< double > _VTXStripAngles;
  std::vector< double > _SITStripAngles;
  std::vector< double > _SETStripAngles;
  
  /** the strip angles for every layer and sensor */
  std::vector< std::vector< double > > _FTDStripAngles;
  
  unsigned _nVTXLayers;
  unsigned _nSITLayers;
  unsigned _nFTDLayers;
  unsigned _nSETLayers;
  
  const gear::ZPlanarParameters* _paramVXD;
  const gear::ZPlanarParameters* _paramSIT;
  const gear::ZPlanarParameters* _paramSET;
  const gear::FTDParameters* _paramFTD;
  
#endif
  
};

#endif

