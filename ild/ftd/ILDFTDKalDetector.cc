
#include "ILDFTDKalDetector.h"

#include "MaterialDataBase.h"

#include <sstream>
#include <iomanip>

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gearimpl/Util.h"
#include "gear/FTDLayerLayout.h"

#include "ILDRotatedTrapMeaslayer.h"
#include "ILDSegmentedDiscMeasLayer.h"
#include "ILDPlanarHit.h"
#include "ILDDiscMeasLayer.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "streamlog/streamlog.h"

#include "TVector3.h"

ILDFTDKalDetector::ILDFTDKalDetector( const gear::GearMgr& gearMgr ) : 
TVKalDetector(300), _nDisks(0), _is_staggered_design(true) // SJA:FIXME initial size, 300 looks reasonable for ILD, though this would be better stored as a const somewhere
{
  
  streamlog_out(DEBUG1) << "ILDFTDKalDetector building FTD detector using GEAR " << std::endl ;
  
  setupGearGeom( gearMgr ) ; 
  
  if( _is_staggered_design ){
    streamlog_out(DEBUG1) << "ILDFTDKalDetector use staggered design" << std::endl ;
    this->build_staggered_design();
  }
  else{
    streamlog_out(DEBUG1) << "ILDFTDKalDetector use tilted design " << std::endl ;
    this->build_turbine_design(); 
  }
  
  SetOwner();
  
}


void ILDFTDKalDetector::create_segmented_disk_layers( int idisk, int nsegments, bool even_petals, double phi0, double zpos , bool front){
  
  Bool_t active = true ;
  Bool_t dummy  = false ;
  
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air") ;
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon") ;
  TMaterial & carbon    = *MaterialDataBase::Instance().getMaterial("carbon") ;
  
  double senThickness = _FTDgeo[idisk].senThickness ;
  double supThickness = _FTDgeo[idisk].supThickness ;
  double innerBaseLength = _FTDgeo[idisk].innerBaseLength ;
  double outerBaseLength = _FTDgeo[idisk].outerBaseLength ;
  double height = _FTDgeo[idisk].height ;
  double rInner = _FTDgeo[idisk].rInner ;
  
  
  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  encoder.reset() ;  // reset to 0
  
  encoder[ILDCellID0::subdet] = ILDDetID::FTD ;
  encoder[ILDCellID0::side] = zpos > 0 ? 1 : -1 ;
  encoder[ILDCellID0::layer]  = idisk ;
  
  int start_index = even_petals ? 0 : 1 ;
  std::vector<int> module_ids;
  
  for (int i=0; i<nsegments; ++i) {
    encoder[ILDCellID0::module] = start_index + i*2 ;
    
    if (front) {
      encoder[ILDCellID0::sensor] = 1 ;
    }
    else{
      encoder[ILDCellID0::sensor] = 3 ;
    }
    
    module_ids.push_back(encoder.lowWord());
    
    if (front) {
      encoder[ILDCellID0::sensor] = 2 ;
    }
    else{
      encoder[ILDCellID0::sensor] = 4 ;
    }
    
    module_ids.push_back(encoder.lowWord());
    
  }
  
  // create segmented disk 
  
  // front face of sensitive  
  double z = zpos - 0.5*(senThickness) ;  
  //  double z = zpos ;
  double sort_policy = fabs(z) ;
  int zsign = z > 0 ? +1 : -1 ;
  
  // if this is the negative z disk add epsilon to the policy
  if( z < 0 ) sort_policy += 1.0e-06 ; 
  
  streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add front face of sensitive at " << z << std::endl;
  Add( new ILDSegmentedDiscMeasLayer(air, silicon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, dummy) );
  
  
  // measurement plane
  z += zsign*0.5*senThickness;  
  sort_policy = z ;
  if( z < 0 ) sort_policy += 1.0e-06 ;
  streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add measurement plane at " << z << " number of module_ids = " << module_ids.size() << std::endl;
  Add( new ILDSegmentedDiscMeasLayer(silicon, silicon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, active, module_ids));
  
  
  // interface between sensitive and support
  z += zsign*0.5*senThickness;   
  sort_policy = z ;
  if( z < 0 ) sort_policy += 1.0e-06 ;
  
  streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add interface between sensitive and support at " << z << std::endl;
  Add( new ILDSegmentedDiscMeasLayer(silicon, carbon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, dummy));
  
  
  // rear face of support
  z += zsign*supThickness;   
  sort_policy = z ;
  if( z < 0 ) sort_policy += 1.0e-06 ;
  
  streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add rear face of support at " << z << std::endl;
  Add( new ILDSegmentedDiscMeasLayer(silicon, carbon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, dummy));
  
  
}


void ILDFTDKalDetector::build_staggered_design() {
  
  streamlog_out(DEBUG) << "ILDFTDKalDetector::build_staggered_design " << std::endl;
  
  
  std::string name = "FTD" ;
  
  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  
  double z_of_last_disc = 0.0 ;
  
  for (int idisk = 0; idisk < _nDisks; ++idisk) {
    
    streamlog_out(DEBUG) << "ILDFTDKalDetector::build_staggered_design build disk " << idisk << std::endl;
    
    int npetals =  _FTDgeo[idisk].nPetals ;
    double phi0 =  _FTDgeo[idisk].phi0 ;
    double senThickness = _FTDgeo[idisk].senThickness ;
    double supThickness = _FTDgeo[idisk].supThickness ;
    
    double senZPos_even_front = _FTDgeo[idisk].senZPos_even_petal1;
    double supZPos_even_front = _FTDgeo[idisk].supZPos_even_petal1;
    
    double senZPos_even_back = _FTDgeo[idisk].senZPos_even_petal3;
    double supZPos_even_back = _FTDgeo[idisk].supZPos_even_petal3;
    
    double senZPos_odd_front = _FTDgeo[idisk].senZPos_odd_petal1;
    double supZPos_odd_front = _FTDgeo[idisk].supZPos_odd_petal1;
    
    double senZPos_odd_back = _FTDgeo[idisk].senZPos_odd_petal3;
    double supZPos_odd_back = _FTDgeo[idisk].supZPos_odd_petal3;
    
    
    // for this design we are assumming that the sensitive and support are share a common boundary
    // here we check if this is the case and exit if not
    double sepration = fabs( senZPos_even_front - supZPos_even_front ) - ( 0.5*senThickness + 0.5*supThickness ) ;
    if( sepration > 1.0e-04 /* 0.1 microns */ ) {
      streamlog_out(ERROR) << "ILDFTDKalDetector design assumes that the sensitive and support are share a common boundary. Separation found to be: "
      << sepration << " microns. exit(1) called" 
      << std::endl ;
      exit(1);
    }
    
    sepration = fabs( senZPos_odd_front - supZPos_odd_front ) - ( 0.5*senThickness + 0.5*supThickness ) ;
    if( sepration > 1.0e-04 /* 0.1 microns */ ) {
      streamlog_out(ERROR) << "ILDFTDKalDetector design assumes that the sensitive and support are share a common boundary. Separation found to be: "
      << sepration << " microns. exit(1) called" 
      << std::endl ;
      exit(1);
    }
    
    // check that the number of petals is divisible by 2
    int nsegments = npetals/2.0;
    
    // even segments forward
    this->create_segmented_disk_layers(idisk, nsegments, true, phi0,  senZPos_even_front, true);
    
    // even segments backwards
    this->create_segmented_disk_layers(idisk, nsegments, true, phi0, -senZPos_even_front, true);
    
    
    // even segments forward
    this->create_segmented_disk_layers(idisk, nsegments, true, phi0,  senZPos_even_back, false);
    
    // even segments backwards
    this->create_segmented_disk_layers(idisk, nsegments, true, phi0, -senZPos_even_back, false);
    
    
    // odd segements 
    // update phi0 by the angular distance of one petal
    phi0 += 2.0 * M_PI / npetals; 
    
    // odd segments forward
    this->create_segmented_disk_layers(idisk, nsegments, false, phi0,  senZPos_odd_front, true);
    
    // odd segments backward
    this->create_segmented_disk_layers(idisk, nsegments, false, phi0, -senZPos_odd_front, true);
    
    // odd segments forward
    this->create_segmented_disk_layers(idisk, nsegments, false, phi0,  senZPos_odd_back, false);
    
    // odd segments backward
    this->create_segmented_disk_layers(idisk, nsegments, false, phi0, -senZPos_odd_back, false);
    
    
    TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air") ;
    
    Bool_t dummy  = false ;
    
    // place air discs to help transport during track extrapolation
    if( idisk != 0 ){
      
      // place the disc half way between the two discs 
      double z = z_of_last_disc + (supZPos_even_front - z_of_last_disc) * 0.5 ;
      
      TVector3 xc_fwd(0.0, 0.0, z) ;
      TVector3 normal_fwd(xc_fwd) ;
      normal_fwd.SetMag(1.0) ;
      
      static const  Double_t eps = 1e-6;
      
      double height = _FTDgeo[idisk].height ;
      double rInner = _FTDgeo[idisk].rInner ;
      
      
      Add(new ILDDiscMeasLayer( air, air, xc_fwd, normal_fwd, _bZ, z,
                               rInner, rInner+height, dummy ) );
      
      TVector3 xc_bwd(0.0, 0.0, -z) ;
      TVector3 normal_bwd(xc_bwd) ;
      normal_bwd.SetMag(1.0) ;
      
      Add(new ILDDiscMeasLayer( air, air, xc_bwd, normal_bwd, _bZ, z+eps,
                               rInner, rInner+height, dummy ) );
      
      
      // save the position of this disc for the next loop
      z_of_last_disc = supZPos_even_front ;   
      
    }
    
  }
  
}


void ILDFTDKalDetector::build_turbine_design() {
  
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air") ;
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon") ;
  TMaterial & carbon    = *MaterialDataBase::Instance().getMaterial("carbon") ;
  
  Bool_t active = true ;
  Bool_t dummy  = false ;
  
  std::string name = "FTD" ;
  
  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  
  
  double z_of_last_disc = 0.0 ;
  
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
    double senZPos = _FTDgeo[idisk].senZPos_even_petal1;
    double supZPos = _FTDgeo[idisk].supZPos_even_petal1;
    
    // for this design we are assumming that the sensitive and support are share a common boundary
    // here we check if this is the case and exit if not
    double sepration = fabs( senZPos - supZPos ) - ( 0.5*senThickness + 0.5*supThickness ) ;
    if( sepration > 1.0e-04 /* 0.1 microns */ ) {
      streamlog_out(ERROR) << "ILDFTDKalDetector design assumes that the sensitive and support are share a common boundary. Separation found to be: "
      << sepration << " microns. exit(1) called" 
      << std::endl ;
      exit(1);
    }
    
    
    double dphi = M_PI / npetals;
    
    for(int ipet=0; ipet<npetals; ++ipet){    
      
      encoder.reset() ;  // reset to 0
      
      encoder[ILDCellID0::subdet] = ILDDetID::FTD ;
      encoder[ILDCellID0::side] = 1 ;
      encoder[ILDCellID0::layer]  = idisk ;
      encoder[ILDCellID0::module] = ipet ;
      encoder[ILDCellID0::sensor] = 0 ;
      
      int CellID_FWD = encoder.lowWord() ;
      
      encoder[ILDCellID0::side] = -1 ;
      
      int CellID_BWD = encoder.lowWord() ;
      
      double cosphi = cos( phi0 + ipet * dphi ) ;
      double sinphi = sin( phi0 + ipet * dphi ) ; 
      double sinalpha = sin( alpha ) ;
      double cosalpha = cos( alpha ) ;
      
      TVector3 sen_front_face_centre_fwd( cosphi * rInner + height*0.5, sinphi * rInner + height*0.5, +senZPos - senThickness*0.5 );         // for +z  
      
      TVector3 measurement_plane_centre_fwd( sen_front_face_centre_fwd.X(), 
                                            sen_front_face_centre_fwd.Y(), 
                                            sen_front_face_centre_fwd.Z() + senThickness*0.5 ); 
      
      TVector3 sen_rear_face_centre_fwd( sen_front_face_centre_fwd.X(), 
                                        sen_front_face_centre_fwd.Y(), 
                                        sen_front_face_centre_fwd.Z() + senThickness ); 
      
      TVector3 sup_rear_face_centre_fwd( sen_rear_face_centre_fwd.X(), 
                                        sen_rear_face_centre_fwd.Y(), 
                                        sen_rear_face_centre_fwd.Z() + supThickness ); 
      
      
      
      // note this is the petal facing the one in +z, not the rotated one 
      TVector3 sen_front_face_centre_bwd( cosphi * rInner + height*0.5, sinphi * rInner + height*0.5, -senZPos + senThickness*0.5 );         // for -z  
      
      TVector3 measurement_plane_centre_bwd( sen_front_face_centre_bwd.X(), 
                                            sen_front_face_centre_bwd.Y(), 
                                            sen_front_face_centre_bwd.Z() - senThickness*0.5 ); 
      
      TVector3 sen_rear_face_centre_bwd( sen_front_face_centre_bwd.X(), 
                                        sen_front_face_centre_bwd.Y(), 
                                        sen_front_face_centre_bwd.Z() - senThickness ); 
      
      TVector3 sup_rear_face_centre_bwd( sen_rear_face_centre_bwd.X(), 
                                        sen_rear_face_centre_bwd.Y(), 
                                        sen_rear_face_centre_bwd.Z() - supThickness ); 
      
      
      TVector3 normalF( sinphi*sinalpha, -cosphi*sinalpha, cosalpha );
      TVector3 normalB( -normalF );
      
      double dist_to_IP =  sen_front_face_centre_fwd.Mag() ;
      
      if(ipet == 0 ){ // need to split the first petal due to overlap
                      // sorting policy is modified for the half which is further from the IP than the last petal in the disc. 
                      // sorting policy = dist_to_IP+3*npetals*eps1 
                      // note it is npetals used here not ipetal 
                      // also the int side argument is set to 1 and -1 as opposed to 0 for the other petals 
        
        // +z
        // air - sensitive boundary
        Add(new ILDRotatedTrapMeaslayer( air, silicon, sen_front_face_centre_fwd, normalF, _bZ, dist_to_IP+4*ipet*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 1, dummy) );
        
        // measurement plane defined as the middle of the sensitive volume
        Add(new ILDRotatedTrapMeaslayer( silicon, silicon, measurement_plane_centre_fwd, normalF, _bZ, dist_to_IP+(4*ipet+1)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 1, active , CellID_FWD) );
        streamlog_out(DEBUG0) << "ILDFTDKalDetector add surface with CellID = "
        << CellID_FWD
        << std::endl ;
        
        // sensitive - support boundary 
        Add(new ILDRotatedTrapMeaslayer( silicon, carbon, sen_rear_face_centre_fwd, normalF, _bZ, dist_to_IP+(4*ipet+2)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 1, dummy ) );
        
        // support - air boundary
        Add(new ILDRotatedTrapMeaslayer( carbon, air, sup_rear_face_centre_fwd, normalF, _bZ, dist_to_IP+(4*ipet+3)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 1, dummy ) );
        
        // +z
        // air - sensitive boundary
        Add(new ILDRotatedTrapMeaslayer( air, silicon, sen_front_face_centre_fwd, normalF, _bZ, dist_to_IP+4*npetals*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, -1, dummy) );
        
        // measurement plane defined as the middle of the sensitive volume
        Add(new ILDRotatedTrapMeaslayer( silicon, silicon, measurement_plane_centre_fwd, normalF, _bZ, dist_to_IP+(4*npetals+1)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, -1, active , CellID_FWD) );
        streamlog_out(DEBUG0) << "ILDFTDKalDetector add surface with CellID = "
        << CellID_FWD
        << std::endl ;
        
        // sensitive - support boundary 
        Add(new ILDRotatedTrapMeaslayer( silicon, carbon, sen_rear_face_centre_fwd, normalF, _bZ, dist_to_IP+(4*npetals+2)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, -1, dummy ) );
        
        // support - air boundary
        Add(new ILDRotatedTrapMeaslayer( carbon, air, sup_rear_face_centre_fwd, normalF, _bZ, dist_to_IP+(4*npetals+3)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, -1, dummy ) );
        
        // -z
        // air - sensitive boundary
        Add(new ILDRotatedTrapMeaslayer( air, silicon, sen_front_face_centre_bwd, normalB, _bZ, dist_to_IP+4*ipet*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, 1, dummy) );
        
        // measurement plane defined as the middle of the sensitive volume
        Add(new ILDRotatedTrapMeaslayer( silicon, silicon, measurement_plane_centre_bwd, normalB, _bZ, dist_to_IP+(4*ipet+1)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 1, active , CellID_BWD) );
        streamlog_out(DEBUG0) << "ILDFTDKalDetector add surface with CellID = "
        << CellID_BWD
        << std::endl ;
        
        // sensitive - support boundary 
        Add(new ILDRotatedTrapMeaslayer( silicon, carbon, sen_rear_face_centre_bwd, normalB, _bZ, dist_to_IP+(4*ipet+2)*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, 1, dummy ) );
        
        // support - air boundary
        Add(new ILDRotatedTrapMeaslayer( carbon, air, sup_rear_face_centre_bwd, normalB, _bZ, dist_to_IP+(4*ipet+3)*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, 1, dummy ) );
        
        
        // -z
        // air - sensitive boundary
        Add(new ILDRotatedTrapMeaslayer( air, silicon, sen_front_face_centre_bwd, normalB, _bZ, dist_to_IP+4*npetals*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, -1, dummy) );
        
        // measurement plane defined as the middle of the sensitive volume
        Add(new ILDRotatedTrapMeaslayer( silicon, silicon, measurement_plane_centre_bwd, normalB, _bZ, dist_to_IP+(4*ipet+1)*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, -1, active , CellID_BWD) );
        streamlog_out(DEBUG0) << "ILDFTDKalDetector add surface with CellID = "
        << CellID_BWD
        << std::endl ;
        
        // sensitive - support boundary 
        Add(new ILDRotatedTrapMeaslayer( silicon, carbon, sen_rear_face_centre_bwd, normalB, _bZ, dist_to_IP+(4*npetals+2)*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, -1, dummy ) );
        
        // support - air boundary
        Add(new ILDRotatedTrapMeaslayer( carbon, air, sup_rear_face_centre_bwd, normalB, _bZ, dist_to_IP+(4*npetals+3)*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, -1, dummy ) );
        
        
        
      }
      else{
        
        // +z
        // air - sensitive boundary
        Add(new ILDRotatedTrapMeaslayer( air, silicon, sen_front_face_centre_fwd, normalF, _bZ, dist_to_IP+4*ipet*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );
        
        // measurement plane defined as the middle of the sensitive volume
        Add(new ILDRotatedTrapMeaslayer( silicon, silicon, measurement_plane_centre_fwd, normalF, _bZ, dist_to_IP+(4*ipet+1)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 0, active , CellID_FWD) );
        streamlog_out(DEBUG0) << "ILDFTDKalDetector add surface with CellID = "
        << CellID_FWD
        << std::endl ;
        
        // sensitive - support boundary 
        Add(new ILDRotatedTrapMeaslayer( silicon, carbon, sen_rear_face_centre_fwd, normalF, _bZ, dist_to_IP+(4*ipet+2)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );
        
        // support - air boundary
        Add(new ILDRotatedTrapMeaslayer( carbon, air, sup_rear_face_centre_fwd, normalF, _bZ, dist_to_IP+(4*ipet+3)*eps1,
                                        height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );
        
        
        // -z
        // air - sensitive boundary
        Add(new ILDRotatedTrapMeaslayer( air, silicon, sen_front_face_centre_bwd, normalB, _bZ, dist_to_IP+4*ipet*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, 0, dummy) );
        
        // measurement plane defined as the middle of the sensitive volume
        Add(new ILDRotatedTrapMeaslayer( silicon, silicon, measurement_plane_centre_bwd, normalB, _bZ, dist_to_IP+(4*ipet+1)*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, 0, active , CellID_BWD) );
        streamlog_out(DEBUG0) << "ILDFTDKalDetector add surface with CellID = "
        << CellID_FWD
        << std::endl ;
        
        // sensitive - support boundary 
        Add(new ILDRotatedTrapMeaslayer( silicon, carbon, sen_rear_face_centre_bwd, normalB, _bZ, dist_to_IP+(4*ipet+2)*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );
        
        // support - air boundary
        Add(new ILDRotatedTrapMeaslayer( carbon, air, sup_rear_face_centre_bwd, normalB, _bZ, dist_to_IP+(4*ipet+3)*eps1+eps2,
                                        height, innerBaseLength, outerBaseLength, alpha, 0, dummy ) );
        
      }
      
    }
    
    // place air discs to help transport during track extrapolation
    if( idisk != 0 ){
      
      // place the disc half way between the two discs 
      double z = z_of_last_disc + (senZPos - z_of_last_disc) * 0.5 ;
      
      TVector3 xc_fwd(0.0, 0.0, z) ;
      TVector3 normal_fwd(xc_fwd) ;
      normal_fwd.SetMag(1.0) ;
      
      Add(new ILDDiscMeasLayer( air, air, xc_fwd, normal_fwd, _bZ, z,
                               rInner, rInner+height, dummy ) );
      
      TVector3 xc_bwd(0.0, 0.0, -z) ;
      TVector3 normal_bwd(xc_bwd) ;
      normal_bwd.SetMag(1.0) ;
      
      Add(new ILDDiscMeasLayer( air, air, xc_bwd, normal_bwd, _bZ, z+eps2,
                               rInner, rInner+height, dummy ) );
      
      
      // save the position of this disc for the next loop
      z_of_last_disc = senZPos ;   
    }
    
  }
  
}


void ILDFTDKalDetector::setupGearGeom( const gear::GearMgr& gearMgr ){
  
  const gear::FTDParameters& ftdParams = gearMgr.getFTDParameters() ;
  const gear::FTDLayerLayout& ftdlayers = ftdParams.getFTDLayerLayout() ;
  
  _bZ = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  _nDisks = ftdlayers.getNLayers() ; // just do the first disk for now 
  _FTDgeo.resize(_nDisks);
  
  //SJA:FIXME: for now the support is taken as the same size the sensitive
  //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
  //           which if traversed in the reverse direction to the next boundary then the track be propagated through carbon
  //           for a significant distance 
  
  for(int disk=0; disk< _nDisks; ++disk){
    
    // numbers taken from the ILD_01 gear file for the sensitive part 
    _FTDgeo[disk].nPetals = ftdlayers.getNPetals(disk) ;    
    _FTDgeo[disk].dphi = 2*M_PI /  _FTDgeo[disk].nPetals ;
    _FTDgeo[disk].phi0 = ftdlayers.getPhi0(disk) ;
    _FTDgeo[disk].alpha = ftdlayers.getAlpha(disk) ;
    _FTDgeo[disk].rInner = ftdlayers.getSensitiveRinner(disk) ;
    _FTDgeo[disk].height = ftdlayers.getSensitiveWidth(disk) ;
    _FTDgeo[disk].innerBaseLength =  ftdlayers.getSensitiveLengthMin(disk) ;
    _FTDgeo[disk].outerBaseLength =  ftdlayers.getSensitiveLengthMax(disk) ;
    _FTDgeo[disk].senThickness =  ftdlayers.getSensitiveThickness(disk) ;
    _FTDgeo[disk].supThickness =  ftdlayers.getSupportThickness(disk) ;
    
    _FTDgeo[disk].senZPos_even_petal1 = ftdlayers.getSensitiveZposition(disk, 0, 1) ; 
    _FTDgeo[disk].senZPos_even_petal2 = ftdlayers.getSensitiveZposition(disk, 0, 2) ; 
    _FTDgeo[disk].senZPos_even_petal3 = ftdlayers.getSensitiveZposition(disk, 0, 3) ; 
    _FTDgeo[disk].senZPos_even_petal4 = ftdlayers.getSensitiveZposition(disk, 0, 4) ; 
    
    // currently the design assumes that the petal on the same side are at the same z
    assert(_FTDgeo[disk].senZPos_even_petal1==_FTDgeo[disk].senZPos_even_petal2);
    assert(_FTDgeo[disk].senZPos_even_petal3==_FTDgeo[disk].senZPos_even_petal4);
    
    _FTDgeo[disk].senZPos_odd_petal1 = ftdlayers.getSensitiveZposition(disk, 1, 1) ; 
    _FTDgeo[disk].senZPos_odd_petal2 = ftdlayers.getSensitiveZposition(disk, 1, 2) ; 
    _FTDgeo[disk].senZPos_odd_petal3 = ftdlayers.getSensitiveZposition(disk, 1, 3) ; 
    _FTDgeo[disk].senZPos_odd_petal4 = ftdlayers.getSensitiveZposition(disk, 1, 4) ; 
    
    // currently the design assumes that the petal on the same side are at the same z
    assert(_FTDgeo[disk].senZPos_odd_petal1==_FTDgeo[disk].senZPos_odd_petal2);
    assert(_FTDgeo[disk].senZPos_odd_petal3==_FTDgeo[disk].senZPos_odd_petal4);
    
    _FTDgeo[disk].supZPos_even_petal1 = ftdlayers.getSensitiveZposition(disk, 0, 1) ; 
    _FTDgeo[disk].supZPos_even_petal2 = ftdlayers.getSensitiveZposition(disk, 0, 2) ; 
    _FTDgeo[disk].supZPos_even_petal3 = ftdlayers.getSensitiveZposition(disk, 0, 3) ; 
    _FTDgeo[disk].supZPos_even_petal4 = ftdlayers.getSensitiveZposition(disk, 0, 4) ; 
    
    assert(_FTDgeo[disk].supZPos_even_petal1==_FTDgeo[disk].supZPos_even_petal2);
    assert(_FTDgeo[disk].supZPos_even_petal3==_FTDgeo[disk].supZPos_even_petal4);
    
    _FTDgeo[disk].supZPos_odd_petal1 = ftdlayers.getSensitiveZposition(disk, 1, 1) ; 
    _FTDgeo[disk].supZPos_odd_petal2 = ftdlayers.getSensitiveZposition(disk, 1, 2) ; 
    _FTDgeo[disk].supZPos_odd_petal3 = ftdlayers.getSensitiveZposition(disk, 1, 3) ; 
    _FTDgeo[disk].supZPos_odd_petal4 = ftdlayers.getSensitiveZposition(disk, 1, 4) ; 
    
    assert(_FTDgeo[disk].supZPos_odd_petal1==_FTDgeo[disk].supZPos_odd_petal2);
    assert(_FTDgeo[disk].supZPos_odd_petal3==_FTDgeo[disk].supZPos_odd_petal4);
    
    
    
    
    // rough check to see if the petal is rotated
    if( fabs(_FTDgeo[disk].alpha) > 1.0e-08  ) { 
      _is_staggered_design = false;
    } else {
      _is_staggered_design = true; 
    }
    
  }
  
  
}
