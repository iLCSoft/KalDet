
#include "ILDVXDKalDetector.h"

#include "MaterialDataBase.h"

#include "ILDParallelPlanarMeasLayer.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include <gear/GEAR.h>
#include "gear/BField.h"
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include "gearimpl/Util.h"

#include "TMath.h"

#include "math.h"
#include <sstream>

#include "streamlog/streamlog.h"

ILDVXDKalDetector::ILDVXDKalDetector( const gear::GearMgr& gearMgr )
: TVKalDetector(300) // SJA:FIXME initial size, 300 looks reasonable for ILD, though this would be better stored as a const somewhere
{
  
  streamlog_out(DEBUG1) << "ILDVXDKalDetector building VXD detector using GEAR " << std::endl ;
  
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon");
  TMaterial & carbon    = *MaterialDataBase::Instance().getMaterial("carbon");
  
  
  this->setupGearGeom(gearMgr) ;
  
  //--The Ladder structure (realistic ladder)--
  int nLadders;
  
  Bool_t active = true;
  Bool_t dummy  = false;
  
  static const double eps = 1e-6; 
  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  
  for (int layer=0; layer<_nLayers; ++layer) {
    
    nLadders = _VXDgeo[layer].nLadders ;
    
    double phi0 = _VXDgeo[layer].phi0 ;
    
    double ladder_distance = _VXDgeo[layer].supRMin ;
    double ladder_thickness = _VXDgeo[layer].supThickness ;
    
    double sensitive_distance = _VXDgeo[layer].senRMin ;
    double sensitive_thickness = _VXDgeo[layer].senThickness ;
    
    double width = _VXDgeo[layer].width ;
    double length = _VXDgeo[layer].length;
    double offset = _VXDgeo[layer].offset;
    
    double pos_xi_nonoverlap_width = (2.0 * (( width / 2.0 ) - fabs(offset))); 
    
    double currPhi;
    double dphi = _VXDgeo[layer].dphi ;
    
    static const double z_offset = 0.0; // all VXD planes are centred at zero 
    
    for (int ladder=0; ladder<nLadders; ++ladder) {
      
      currPhi = phi0 + (dphi * ladder);
      
      encoder.reset() ;  // reset to 0
      
      encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::VXD ;
      encoder[lcio::ILDCellID0::side] = 0 ;
      encoder[lcio::ILDCellID0::layer]  = layer ;
      encoder[lcio::ILDCellID0::module] = ladder ;
      encoder[lcio::ILDCellID0::sensor] = 0 ;
      
      int CellID = encoder.lowWord() ;
      
      
      if(layer%2 == 0 ){ // overlap section of ladder0 is defined after the last ladder,
        
        
        double sen_front_sorting_policy         = sensitive_distance  + (4 * ladder+0) * eps ;
        double measurement_plane_sorting_policy = sensitive_distance  + (4 * ladder+1) * eps ;
        double sen_back_sorting_policy          = sensitive_distance  + (4 * ladder+2) * eps ;
        double sup_back_sorting_policy          = ladder_distance     + (4 * ladder+4) * eps ;
        
        
        if(ladder==0){   // bacause overlap section of ladder0 is further outer than the last ladder.
          
          streamlog_out(DEBUG0) << "ILDVXDKalDetector add surface with CellID = "
          << CellID
          << std::endl ;
          
          // non overlapping region
          // air - sensitive boundary
          Add(new ILDParallelPlanarMeasLayer(air, silicon, sensitive_distance, currPhi, _bZ, sen_front_sorting_policy, pos_xi_nonoverlap_width, length, 0.0, z_offset, offset, dummy,-1,"VXDSenFront" )) ;
          
          // measurement plane defined as the middle of the sensitive volume 
          Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, pos_xi_nonoverlap_width, length, 0.0, z_offset, offset, active, CellID, "VXDMeasLayer" )) ;          
          
          // sensitive - support boundary 
          Add(new ILDParallelPlanarMeasLayer(silicon, carbon, sensitive_distance+sensitive_thickness, currPhi, _bZ, sen_back_sorting_policy, pos_xi_nonoverlap_width, length, 0.0, z_offset, offset, dummy,-1,"VXDSenSuppportIntf" )) ; 
          
          // support - air boundary
          Add(new ILDParallelPlanarMeasLayer(carbon, air, ladder_distance+ladder_thickness, currPhi, _bZ, sup_back_sorting_policy, pos_xi_nonoverlap_width, length, 0.0, z_offset, offset, dummy,-1,"VXDSupRear" )) ;           
          
          
          // overlapping region
          double overlap_region_width  = width - pos_xi_nonoverlap_width ;
          double overlap_region_offset = -(overlap_region_width/2.0) - (pos_xi_nonoverlap_width)/2.0 ;

          
          // overlap sorting policy uses nLadders as the overlapping "ladder" is the order i.e. there will now be nLadders+1 
          double overlap_front_sorting_policy                = sensitive_distance + (4* nLadders+0) * eps;
          double overlap_measurement_plane_sorting_policy    = sensitive_distance + (4* nLadders+1) * eps;
          double overlap_back_sorting_policy                 = sensitive_distance + (4* nLadders+2) * eps;
          double overlap_sup_back_sorting_policy             = sensitive_distance + (4* nLadders+3) * eps;
          
          streamlog_out(DEBUG0) << "ILDVXDKalDetector add surface with CellID = "
          << CellID
          << std::endl ;
          
          // air - sensitive boundary
          Add(new ILDParallelPlanarMeasLayer(air, silicon, sensitive_distance, currPhi, _bZ, overlap_front_sorting_policy, overlap_region_width, length, overlap_region_offset, z_offset, offset, dummy,-1,"VXDSenFront")) ;
          
          // measurement plane defined as the middle of the sensitive volume
          Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, overlap_measurement_plane_sorting_policy, overlap_region_width, length, overlap_region_offset, z_offset, offset, active, CellID, "VXDMeasLayer" )) ;
          
          
          // sensitive - support boundary 
          Add(new ILDParallelPlanarMeasLayer(silicon, carbon, sensitive_distance+sensitive_thickness, currPhi, _bZ, overlap_back_sorting_policy, overlap_region_width, length, overlap_region_offset, z_offset, offset, dummy,-1,"VXDSenSuppportIntf")) ; 
          
          // support - air boundary
          Add(new ILDParallelPlanarMeasLayer(carbon, air, ladder_distance+ladder_thickness, currPhi, _bZ, overlap_sup_back_sorting_policy, overlap_region_width, length, overlap_region_offset, z_offset, offset, dummy,-1,"VXDSupRear")) ; 
          
        }
        else{
          
          streamlog_out(DEBUG0) << "ILDVXDKalDetector (ILDParallelPlanarMeasLayer) add surface with CellID = "
          << CellID
          << std::endl ;                                        
          
          // air - sensitive boundary
          Add(new ILDParallelPlanarMeasLayer(air, silicon, sensitive_distance, currPhi, _bZ, sen_front_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSenFront")) ;
          
          // measurement plane defined as the middle of the sensitive volume
          Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, length, offset, z_offset, offset, active, CellID, "VXDMeaslayer" )) ;
          
          // sensitive - support boundary 
          Add(new ILDParallelPlanarMeasLayer(silicon, carbon, sensitive_distance+sensitive_thickness, currPhi, _bZ, sen_back_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSenSuppportIntf" )) ; 
          
          // support - air boundary
          Add(new ILDParallelPlanarMeasLayer(carbon, air, ladder_distance+ladder_thickness, currPhi, _bZ, sup_back_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSupRear" )) ; 
          
        }        
      }
      else{ // counting from 0, odd numbered layers are placed with the support closer to the IP than the sensitive
                                                                
        double sup_forward_sorting_policy        = ladder_distance + (4 * ladder+0) * eps ;
        double sup_back_sorting_policy           = ladder_distance + (4 * ladder+1) * eps ;
        double measurement_plane_sorting_policy  = ladder_distance + (4 * ladder+2) * eps ;
        double sen_back_sorting_policy           = ladder_distance + (4 * ladder+3) * eps ;
        
        streamlog_out(DEBUG0) << "ILDVXDKalDetector (ILDPlanarMeasLayer) add surface with CellID = "
        << CellID
        << std::endl ;
        
        // air - support boundary
        Add(new ILDParallelPlanarMeasLayer(air, carbon, ladder_distance, currPhi, _bZ, sup_forward_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSupFront" )) ; 
        
        // support - sensitive boundary 
        Add(new ILDParallelPlanarMeasLayer(carbon, silicon, (ladder_distance+ladder_thickness), currPhi, _bZ, sup_back_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSenSuppportIntf")) ; 
        
        // measurement plane defined as the middle of the sensitive volume
        Add(new ILDParallelPlanarMeasLayer(silicon, silicon, (sensitive_distance+sensitive_thickness*0.5), currPhi, _bZ, measurement_plane_sorting_policy, width, length, offset, z_offset, offset, active, CellID, "VXDMeaslayer")) ; 
        
        // sensitive air - sensitive boundary
        Add(new ILDParallelPlanarMeasLayer(silicon, air, (sensitive_distance+sensitive_thickness), currPhi, _bZ, sen_back_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSenRear")) ;
        
        
      }
    }
  }
  
  SetOwner();                   
}




void ILDVXDKalDetector::setupGearGeom( const gear::GearMgr& gearMgr ){
  
  const gear::VXDParameters& pVXDDetMain = gearMgr.getVXDParameters();
  const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();
  
  _bZ = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  _nLayers = pVXDLayerLayout.getNLayers(); 
  _VXDgeo.resize(_nLayers);
  
  //SJA:FIXME: for now the support is taken as the same size the sensitive
  //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
  //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
  //           for a significant distance 
  
  for( int layer=0; layer < _nLayers; ++layer){
    _VXDgeo[layer].nLadders = pVXDLayerLayout.getNLadders(layer); 
    _VXDgeo[layer].phi0 = pVXDLayerLayout.getPhi0(layer); 
    _VXDgeo[layer].dphi = 2*M_PI / _VXDgeo[layer].nLadders; 
    _VXDgeo[layer].senRMin = pVXDLayerLayout.getSensitiveDistance(layer); 
    _VXDgeo[layer].supRMin = pVXDLayerLayout.getLadderDistance(layer); 
    _VXDgeo[layer].length = pVXDLayerLayout.getSensitiveLength(layer) * 2.0 ; // note: gear for historical reasons uses the halflength 
    _VXDgeo[layer].width = pVXDLayerLayout.getSensitiveWidth(layer); 
    _VXDgeo[layer].offset = pVXDLayerLayout.getSensitiveOffset(layer); 
    _VXDgeo[layer].senThickness = pVXDLayerLayout.getSensitiveThickness(layer); 
    _VXDgeo[layer].supThickness = pVXDLayerLayout.getLadderThickness(layer); 
  }
  
  
  
  
  
}
