
#include "ILDVXDKalDetector.h"

#include "MaterialDataBase.h"

#include "ILDPlanarMeasLayer.h"
#include "ILDPlanarHit.h"

#include <ILDDetectorIDs.h>

#include <gear/GEAR.h>
#include "gear/BField.h"
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include "gearimpl/Util.h"

#include "TMath.h"

#include "math.h"
#include <sstream>

ILDVXDKalDetector::ILDVXDKalDetector( const gear::GearMgr& gearMgr )
                : TVKalDetector(6) // SJA:FIXME initial size, 6 looks reasonable for ILD, though this would be better stored as a const somewhere
{
   

  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon");
  TMaterial & carbon   = *MaterialDataBase::Instance().getMaterial("carbon");

  const gear::VXDParameters& pVXDDetMain = gearMgr.getVXDParameters();
  const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();

  const Double_t bz = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;

  //--Get additional parameters from GEAR--
  const int nLayersVTX = pVXDLayerLayout.getNLayers();

  //--The Ladder structure (realistic ladder)--
  int nLadders;
  
  Bool_t active = true;
  Bool_t dummy  = false;
  
  static const Double_t eps = 1e-6; 

  for (int layer=0; layer<nLayersVTX; ++layer) {

    nLadders = pVXDLayerLayout.getNLadders(layer);
    
    float ladder_phi0 = float(pVXDLayerLayout.getPhi0(layer));

    float ladder_distance = float(pVXDLayerLayout.getLadderDistance(layer));
    float ladder_thickness = float(pVXDLayerLayout.getLadderThickness(layer));
    float ladder_width = float(pVXDLayerLayout.getLadderWidth(layer));
    float ladder_length = float (pVXDLayerLayout.getLadderLength(layer));
    float ladder_offset = float (pVXDLayerLayout.getLadderOffset(layer));
    
    float sensitive_distance = float(pVXDLayerLayout.getSensitiveDistance(layer));
    float sensitive_thickness = float(pVXDLayerLayout.getSensitiveThickness(layer));
    float sensitive_width = float(pVXDLayerLayout.getSensitiveWidth(layer));
    float sensitive_length = float(pVXDLayerLayout.getSensitiveLength(layer));
    float sensitive_offset = float (pVXDLayerLayout.getSensitiveOffset(layer));
    
    Double_t pos_xi_nonoverlap_width = sensitive_distance*tan( M_PI / nLadders);
    Double_t ladder_xi_min = ladder_offset - sensitive_width/2.0 ;
    Double_t ladder_xi_max = ladder_offset + sensitive_width/2.0 ;

    float currPhi;
    float angleLadders = 2*M_PI / nLadders;
    float cosphi, sinphi;

    for (int ladder=0; ladder<nLadders; ++ladder) {
      
      currPhi = ladder_phi0 + (angleLadders * ladder);
      cosphi = cos(currPhi);
      sinphi = sin(currPhi);

      //Add(new EXVTXMeasLayer(air, si, xc, normal, rmin + 2*ladder*eps, plyhwidth - sximin, 2*hlength, (plyhwidth + sximin)/2, sigmaxi, sigmazeta, inner,ss.str().data()));

      TVector3 normal( cosphi, sinphi, 0) ;

      TVector3 sen_front_face_centre( sensitive_distance*cosphi, sensitive_distance*sinphi, 0) ; 
      TVector3 sen_back_face_centre( (sensitive_distance+sensitive_thickness)*cosphi, (sensitive_distance+sensitive_thickness)*sinphi, 0) ; 
      
      Double_t sen_front_sorting_policy = sensitive_distance + (2 * ladder) * eps ;
      Double_t sen_back_sorting_policy = sensitive_distance + (2 * ladder+1) * eps ;

      int layerID = ILDDetectorIDs::DetID::VXD * ILDDetectorIDs::DetID::Factor  + layer * 100 + ladder;

      if(layer%2 == 0 ){ // overlap section of ladder0 is defined after the last ladder,
	if(ladder==0){   // bacause overlap section of ladder0 is further outer than the last ladder.

	  // non overlapping region
	  Add(new ILDPlanarMeasLayer(air, silicon, sen_front_face_centre, normal, bz, sen_front_sorting_policy, pos_xi_nonoverlap_width - ladder_xi_min, sensitive_length, (pos_xi_nonoverlap_width + ladder_xi_min)/2, active, layerID )) ;
	  Add(new ILDPlanarMeasLayer(silicon, air, sen_back_face_centre, normal, bz, sen_back_sorting_policy, pos_xi_nonoverlap_width - ladder_xi_min, sensitive_length, (pos_xi_nonoverlap_width + ladder_xi_min)/2, dummy)) ; // back side declared not sensitive	  

	  // overlapping region
	  Double_t overlap_region_width  = sensitive_width + ladder_xi_min - pos_xi_nonoverlap_width ;
	  Double_t overlap_region_offset = (ladder_xi_max + pos_xi_nonoverlap_width)/2 ;
	  
	  Double_t overlap_front_sorting_policy = sensitive_distance + 2*nLadders*eps;
	  Double_t overlap_back_sorting_policy  = sensitive_distance + (2*nLadders+1)*eps;

	  Add(new ILDPlanarMeasLayer(air, silicon, sen_front_face_centre, normal, bz, overlap_front_sorting_policy, overlap_region_width, sensitive_length, overlap_region_offset, active, layerID )) ;
	  Add(new ILDPlanarMeasLayer(silicon, air, sen_back_face_centre, normal, bz, overlap_back_sorting_policy, overlap_region_width, sensitive_length, overlap_region_offset, dummy)) ; // back side declared not sensitive	  
	  
	}
	else{
	  
	  Add(new ILDPlanarMeasLayer(air, silicon, sen_front_face_centre, normal, bz, sen_front_sorting_policy, sensitive_width, sensitive_length, sensitive_offset, active, layerID )) ;
	  Add(new ILDPlanarMeasLayer(silicon, air, sen_back_face_centre, normal, bz, sen_back_sorting_policy, sensitive_width, sensitive_length, sensitive_offset, dummy )) ; // back side declared not sensitive

	}	 
      }
      else{ //SJA::FIXME: in mokka the hits will be placed on the surface facing the IP so the inner surface is defined as active here 

	Add(new ILDPlanarMeasLayer(air, silicon, sen_front_face_centre, normal, bz, sen_front_sorting_policy, sensitive_width, sensitive_length, sensitive_offset, active, layerID )) ;
	Add(new ILDPlanarMeasLayer(silicon, air, sen_back_face_centre, normal, bz, sen_back_sorting_policy, sensitive_width, sensitive_length, sensitive_offset, dummy )) ; // back side declared not sensitive

      }
    }
  }

   SetOwner();			
}

ILDVXDKalDetector::~ILDVXDKalDetector()
{
}
