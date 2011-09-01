
#include "ILDVXDKalDetector.h"

#include "MaterialDataBase.h"

#include "ILDPlanarMeasLayer.h"
#include "ILDPlanarHit.h"

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
   
  streamlog_out(DEBUG4) << "ILDVXDKalDetector building VXD detector using GEAR " << std::endl ;

  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon");
  TMaterial & carbon   = *MaterialDataBase::Instance().getMaterial("carbon");


  this->setupGearGeom(gearMgr) ;

  //--The Ladder structure (realistic ladder)--
  int nLadders;
  
  Bool_t active = true;
  Bool_t dummy  = false;
  
  static const double eps = 1e-6; 

  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
 
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
    
    double pos_xi_nonoverlap_width = sensitive_distance*tan( M_PI / nLadders);
    double ladder_xi_min = offset - width/2.0 ;
    double ladder_xi_max = offset + width/2.0 ;

    double currPhi;
    double dphi = _VXDgeo[layer].dphi ;
    double cosphi, sinphi;

    for (int ladder=0; ladder<nLadders; ++ladder) {
      
      
      currPhi = phi0 + (dphi * ladder);
      cosphi = cos(currPhi);
      sinphi = sin(currPhi);

      TVector3 normal( cosphi, sinphi, 0) ;
      
      encoder.reset() ;  // reset to 0

      encoder[ILDCellID0::subdet] = ILDDetID::VXD ;
      encoder[ILDCellID0::side] = 0 ;
      encoder[ILDCellID0::layer]  = layer ;
      encoder[ILDCellID0::module] = ladder ;
      encoder[ILDCellID0::sensor] = 0 ;

      int layerID = encoder.lowWord() ;
      
      
      if(layer%2 == 0 ){ // overlap section of ladder0 is defined after the last ladder,

      TVector3 sen_front_face_centre( sensitive_distance*cosphi, sensitive_distance*sinphi, 0) ; 
      TVector3 measurement_plane_centre( (sensitive_distance+sensitive_thickness*0.5)*cosphi, (sensitive_distance+sensitive_thickness*0.5)*sinphi, 0) ; 
      TVector3 sen_back_face_centre( (sensitive_distance+sensitive_thickness)*cosphi, (sensitive_distance+sensitive_thickness)*sinphi, 0) ; 
      TVector3 sup_back_face_centre( (ladder_distance+ladder_thickness)*cosphi, (ladder_distance+ladder_thickness)*sinphi, 0) ; 
      
      double sen_front_sorting_policy         = sensitive_distance  + (4 * ladder+0) * eps ;
      double measurement_plane_sorting_policy = sensitive_distance  + (4 * ladder+1) * eps ;
      double sen_back_sorting_policy          = sensitive_distance  + (4 * ladder+2) * eps ;
      double sup_back_sorting_policy          = ladder_distance     + (4 * ladder+4) * eps ;


	if(ladder==0){   // bacause overlap section of ladder0 is further outer than the last ladder.
	  
	  
	  // non overlapping region
	  // air - sensitive boundary
	  Add(new ILDPlanarMeasLayer(air, silicon, sen_front_face_centre, normal, _bZ, sen_front_sorting_policy, pos_xi_nonoverlap_width - ladder_xi_min, length, (pos_xi_nonoverlap_width + ladder_xi_min)/2, dummy )) ;

	  // measurement plane defined as the middle of the sensitive volume 
	  Add(new ILDPlanarMeasLayer(silicon, silicon, measurement_plane_centre, normal, _bZ, measurement_plane_sorting_policy, pos_xi_nonoverlap_width - ladder_xi_min, length, (pos_xi_nonoverlap_width + ladder_xi_min)/2, active, layerID )) ;
	  streamlog_out(DEBUG3) << "ILDVXDKalDetector add surface with layerID = "
			    << layerID
			    << std::endl ;


	  // sensitive - support boundary 
	  Add(new ILDPlanarMeasLayer(silicon, carbon, sen_back_face_centre, normal, _bZ, sen_back_sorting_policy, pos_xi_nonoverlap_width - ladder_xi_min, length, (pos_xi_nonoverlap_width + ladder_xi_min)/2, dummy)) ;

	  // support - air boundary
 	  Add(new ILDPlanarMeasLayer(carbon, air, sup_back_face_centre, normal, _bZ, sup_back_sorting_policy, pos_xi_nonoverlap_width - ladder_xi_min, length, (pos_xi_nonoverlap_width + ladder_xi_min)/2, dummy)) ; 

	  // overlapping region
	  double overlap_region_width  = width + ladder_xi_min - pos_xi_nonoverlap_width ;
	  double overlap_region_offset = (ladder_xi_max + pos_xi_nonoverlap_width)/2 ;
	  
	  // overlap sorting policy uses nLadders as the overlapping "ladder" is the order i.e. there will now be nLadders+1 
	  double overlap_front_sorting_policy                = sensitive_distance + (4* nLadders+0) * eps;
	  double overlap_measurement_plane_sorting_policy    = sensitive_distance + (4* nLadders+1) * eps;
	  double overlap_back_sorting_policy                 = sensitive_distance + (4* nLadders+2) * eps;
	  double overlap_sup_back_sorting_policy             = sensitive_distance + (4* nLadders+3) * eps;

	  // air - sensitive boundary
	  Add(new ILDPlanarMeasLayer(air, silicon, sen_front_face_centre, normal, _bZ, overlap_front_sorting_policy, overlap_region_width, length, overlap_region_offset, dummy)) ;

	  // measurement plane defined as the middle of the sensitive volume
	  Add(new ILDPlanarMeasLayer(silicon, silicon, measurement_plane_centre, normal, _bZ, overlap_measurement_plane_sorting_policy, overlap_region_width, length, overlap_region_offset, active, layerID )) ;
	  streamlog_out(DEBUG3) << "ILDVXDKalDetector add surface with layerID = "
				<< layerID
				<< std::endl ;


	  // sensitive - support boundary 
	  Add(new ILDPlanarMeasLayer(silicon, air, sen_back_face_centre, normal, _bZ, overlap_back_sorting_policy, overlap_region_width, length, overlap_region_offset, dummy)) ; 
	  
	  // support - air boundary
	  Add(new ILDPlanarMeasLayer(carbon, air, sup_back_face_centre, normal, _bZ, overlap_sup_back_sorting_policy, overlap_region_width, length, overlap_region_offset, dummy)) ; 

	}
	else{

	  // air - sensitive boundary
	  Add(new ILDPlanarMeasLayer(air, silicon, sen_front_face_centre, normal, _bZ, sen_front_sorting_policy, width, length, offset, dummy )) ;

	  // measurement plane defined as the middle of the sensitive volume
	  Add(new ILDPlanarMeasLayer(silicon, silicon, measurement_plane_centre, normal, _bZ, measurement_plane_sorting_policy, width, length, offset, active, layerID )) ;
	  streamlog_out(DEBUG3) << "ILDVXDKalDetector add surface with layerID = "
				<< layerID
				<< std::endl ;

	  // sensitive - support boundary 
	  Add(new ILDPlanarMeasLayer(silicon, carbon, sen_back_face_centre, normal, _bZ, sen_back_sorting_policy, width, length, offset, dummy )) ; 

	  // support - air boundary
	  Add(new ILDPlanarMeasLayer(carbon, air, sup_back_face_centre, normal, _bZ, sup_back_sorting_policy, width, length, offset, dummy )) ; 

	}	 
      }
      else{ // counting from 0, odd numbered layers are placed with the support closer to the IP than the sensitive

	TVector3 sup_front_face_centre( ladder_distance*cosphi, ladder_distance*sinphi, 0) ; 
	TVector3 sup_back_face_centre( (ladder_distance+ladder_thickness)*cosphi, (ladder_distance+ladder_thickness)*sinphi, 0) ; 
	TVector3 measurement_plane_centre( (sensitive_distance+sensitive_thickness*0.5)*cosphi, (sensitive_distance+sensitive_thickness*0.5)*sinphi, 0) ; 
	TVector3 sen_back_face_centre( (sensitive_distance+sensitive_thickness)*cosphi, (sensitive_distance+sensitive_thickness)*sinphi, 0) ; 
	
	double sup_forward_sorting_policy        = ladder_distance + (4 * ladder+0) * eps ;
	double sup_back_sorting_policy           = ladder_distance + (4 * ladder+1) * eps ;
	double measurement_plane_sorting_policy  = ladder_distance + (4 * ladder+2) * eps ;
	double sen_back_sorting_policy           = ladder_distance + (4 * ladder+3) * eps ;
	

	// air - support boundary
	Add(new ILDPlanarMeasLayer(air, carbon, sup_front_face_centre, normal, _bZ, sup_forward_sorting_policy, width, length, offset, dummy )) ; 

	// support - sensitive boundary 
	Add(new ILDPlanarMeasLayer(carbon, silicon, sup_back_face_centre, normal, _bZ, sup_back_sorting_policy, width, length, offset, dummy)) ; 

	// measurement plane defined as the middle of the sensitive volume
	Add(new ILDPlanarMeasLayer(silicon, silicon, measurement_plane_centre, normal, _bZ, measurement_plane_sorting_policy, width, length, offset, active, layerID )) ; 
	streamlog_out(DEBUG3) << "ILDVXDKalDetector add surface with layerID = "
			      << layerID
			      << std::endl ;

	// sensitive air - sensitive boundary
	Add(new ILDPlanarMeasLayer(silicon, air, sen_back_face_centre, normal, _bZ, sen_back_sorting_policy, width, length, offset, dummy)) ;


      }
    }
  }

   SetOwner();			
}

ILDVXDKalDetector::~ILDVXDKalDetector()
{
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
    _VXDgeo[layer].length = pVXDLayerLayout.getSensitiveLength(layer); 
    _VXDgeo[layer].width = pVXDLayerLayout.getSensitiveWidth(layer); 
    _VXDgeo[layer].offset = pVXDLayerLayout.getSensitiveOffset(layer); 
    _VXDgeo[layer].senThickness = pVXDLayerLayout.getSensitiveThickness(layer); 
    _VXDgeo[layer].supThickness = pVXDLayerLayout.getLadderThickness(layer); 
    }





}
