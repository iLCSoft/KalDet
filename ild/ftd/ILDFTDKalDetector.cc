
#include "ILDFTDKalDetector.h"

#include "MaterialDataBase.h"

#include <sstream>
#include <iomanip>

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gearimpl/Util.h"
#include "gear/FTDLayerLayout.h"

#include "ILDRotatedTrapMeaslayer.h"
#include "ILDPlanarHit.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "streamlog/streamlog.h"

#include "TVector3.h"

ILDFTDKalDetector::ILDFTDKalDetector( const gear::GearMgr& gearMgr ) : 
  TVKalDetector(300), _nDisks(0) // SJA:FIXME initial size, 300 looks reasonable for ILD, though this would be better stored as a const somewhere
{

  streamlog_out(DEBUG4) << "ILDFTDKalDetector building FTD detector using GEAR " << std::endl ;

  setupGearGeom( gearMgr ) ; 

  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air") ;
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon") ;
  TMaterial & carbon    = *MaterialDataBase::Instance().getMaterial("carbon") ;

  Bool_t active = true ;
  Bool_t dummy  = false ;
  
  std::string name = "FTD" ;
  
  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 


  for (int idisk = 0; idisk < _nDisks; ++idisk) {

    static const  Double_t eps1 = 1e-6;
    static const  Double_t eps2 = 1e-4;

    int npetals =  _FTDgeo[idisk].nPetals ;
    double phi0 =  _FTDgeo[idisk].phi0 ;
    double alpha  = _FTDgeo[idisk].alpha ; 
    double height = _FTDgeo[idisk].height ;
    double rInner = _FTDgeo[idisk].rInner ;
    double innerBaseLength = _FTDgeo[idisk].innerBaseLength ;
    double outerBaseLength = _FTDgeo[idisk].outerBaseLength ;
    double senThickness = _FTDgeo[idisk].senThickness ;
    double supThickness = _FTDgeo[idisk].supThickness ;
    double zPos = _FTDgeo[idisk].zPos;

    double dphi = M_PI / npetals;
   
    for(int ipet=0; ipet<npetals; ++ipet){    

      encoder.reset() ;  // reset to 0
      
      encoder[ILDCellID0::subdet] = ILDDetID::FTD ;
      encoder[ILDCellID0::side] = 1 ;
      encoder[ILDCellID0::layer]  = idisk ;
      encoder[ILDCellID0::module] = ipet ;
      encoder[ILDCellID0::sensor] = 0 ;

      int CELL_ID_FWD = encoder.lowWord() ;

      encoder[ILDCellID0::side] = -1 ;

      int CELL_ID_BWD = encoder.lowWord() ;

      double cosphi = cos( phi0 + ipet * dphi ) ;
      double sinphi = sin( phi0 + ipet * dphi ) ; 
      double sinalpha = sin( alpha ) ;
      double cosalpha = cos( alpha ) ;

      TVector3 sen_front_face_centre_fwd( cosphi * rInner + height*0.5, sinphi * rInner + height*0.5, +zPos );         // for +z  

      TVector3 sen_rear_face_centre_fwd( sen_front_face_centre_fwd.X(), 
					 sen_front_face_centre_fwd.Y(), 
					 sen_front_face_centre_fwd.Z() + senThickness ); 

      TVector3 sup_rear_face_centre_fwd( sen_rear_face_centre_fwd.X(), 
					 sen_rear_face_centre_fwd.Y(), 
					 sen_rear_face_centre_fwd.Z() + supThickness ); 



      // note this is the petal facing the one in +z, not the rotated one 
      TVector3 sen_front_face_centre_bwd( cosphi * rInner + height*0.5, sinphi * rInner + height*0.5, -zPos );         // for -z  

      TVector3 sen_rear_face_centre_bwd( sen_front_face_centre_bwd.X(), 
					 sen_front_face_centre_bwd.Y(), 
					 sen_front_face_centre_bwd.Z() - senThickness ); 

      TVector3 sup_rear_face_centre_bwd( sen_rear_face_centre_bwd.X(), 
					 sen_rear_face_centre_bwd.Y(), 
					 sen_rear_face_centre_bwd.Z() - supThickness ); 


      TVector3 normalF( sinphi*sinalpha, -cosphi*sinalpha, cosalpha );
      TVector3 normalB( -normalF );

      double dist_to_IP =  sen_front_face_centre_fwd.Mag() ;

      // +z
      // air - sensitive boundary
      Add(new ILDRotatedTrapMeaslayer( air, silicon, sen_front_face_centre_fwd, normalF, _bZ, dist_to_IP+3*ipet*eps1,
				       height, innerBaseLength, outerBaseLength, alpha, 0, active , CELL_ID_FWD) );
      streamlog_out(DEBUG3) << "ILDFTDKalDetector add surface with layerID = "
			    << CELL_ID_FWD
			    << std::endl ;
      // sensitive - support boundary 
      Add(new ILDRotatedTrapMeaslayer( silicon, carbon, sen_rear_face_centre_fwd, normalF, _bZ, dist_to_IP+(3*ipet+1)*eps1,
				       height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );

      // support - air boundary
      Add(new ILDRotatedTrapMeaslayer( carbon, air, sup_rear_face_centre_fwd, normalF, _bZ, dist_to_IP+(3*ipet+2)*eps1,
				       height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );

      
      // -z
      // air - sensitive boundary
      Add(new ILDRotatedTrapMeaslayer( air, silicon, sen_front_face_centre_bwd, normalB, _bZ, dist_to_IP+3*ipet*eps1+eps2,
				       height, innerBaseLength, outerBaseLength, alpha, 0, active , CELL_ID_BWD) );
      streamlog_out(DEBUG3) << "ILDFTDKalDetector add surface with layerID = "
			    << CELL_ID_BWD
			    << std::endl ;
      // sensitive - support boundary 
      Add(new ILDRotatedTrapMeaslayer( silicon, carbon, sen_rear_face_centre_bwd, normalB, _bZ, dist_to_IP+(3*ipet+1)*eps1+eps2,
				       height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );

      // support - air boundary
      Add(new ILDRotatedTrapMeaslayer( carbon, air, sup_rear_face_centre_bwd, normalB, _bZ, dist_to_IP+(3*ipet+2)*eps1+eps2,
				       height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );

      
    }


  }

  SetOwner();

}
  


void ILDFTDKalDetector::setupGearGeom( const gear::GearMgr& gearMgr ){

  const gear::FTDParameters& ftdParams = gearMgr.getFTDParameters() ;
  const gear::FTDLayerLayout& ftdlayers = ftdParams.getFTDLayerLayout() ;
  
  _bZ = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;

  _nDisks = 1 ; // just do the first disk for now 
  _FTDgeo.resize(_nDisks);

  //SJA:FIXME: for now the support is taken as the same size the sensitive
  //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
  //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
  //           for a significant distance 
  
  for(int disk=0; disk< _nDisks; ++disk){
    
    // numbers taken from the ILD_01 gear file for the sensitive part 
    _FTDgeo[disk].nPetals = ftdlayers.getNLadders(disk) ;    
    _FTDgeo[disk].dphi = 2*M_PI /  _FTDgeo[disk].nPetals ;
    _FTDgeo[disk].phi0 = ftdlayers.getPhi(disk) + 0.5 * _FTDgeo[disk].dphi ; // SJA:FIXME: phi0 is not well defined, the gear file has it as the phi angle of the trailing edge of the first petal
    _FTDgeo[disk].alpha = ftdlayers.getAlpha(disk) ;
    _FTDgeo[disk].rInner = ftdlayers.getSensitiveRinner(disk) ;
    _FTDgeo[disk].height = ftdlayers.getSensitiveWidth(disk) ;
    _FTDgeo[disk].innerBaseLength =  ftdlayers.getSensitiveLengthMin(disk) ;
    _FTDgeo[disk].outerBaseLength =  ftdlayers.getSensitiveLengthMax(disk) ;
    _FTDgeo[disk].senThickness =  ftdlayers.getSensitiveThickness(disk) ;
    _FTDgeo[disk].supThickness =  ftdlayers.getLadderThickness(disk) ;
    _FTDgeo[disk].zPos = ftdlayers.getZposition(disk) ;


  }


}
