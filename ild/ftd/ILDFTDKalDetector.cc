
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
TVKalDetector(300), _nDisks(0) // SJA:FIXME initial size, 300 looks reasonable for ILD, though this would be better stored as a const somewhere
{
  
  streamlog_out(DEBUG1) << "ILDFTDKalDetector building FTD detector using GEAR " << std::endl ;
  
  setupGearGeom( gearMgr ) ; 
  
  this->build_staggered_design();

  
  SetOwner();
  
}




void ILDFTDKalDetector::build_staggered_design() {
  
  streamlog_out(DEBUG) << "ILDFTDKalDetector::build_staggered_design " << std::endl;
  
  
  std::string name = "FTD" ;
  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  
  
  
  for (int idisk = 0; idisk < _nDisks; ++idisk) {
    
    streamlog_out(DEBUG) << "ILDFTDKalDetector::build_staggered_design build disk " << idisk << std::endl;
    
    int npetals =  _FTDgeo[idisk].nPetals ;
    double phi0 =  _FTDgeo[idisk].phi0 ;
    
    double senZPos_even_front = _FTDgeo[idisk].senZPos_even_front;
    double senZPos_odd_front = _FTDgeo[idisk].senZPos_odd_front;
    
    
    // check that the number of petals is divisible by 2
    int nsegments = npetals/2.0;
    
    // even segments forward
    this->create_segmented_disk_layers(idisk, nsegments, true, phi0,  senZPos_even_front );
    
    // even segments backwards
    this->create_segmented_disk_layers(idisk, nsegments, true, phi0, -senZPos_even_front );
    
 
    // odd segements 
    // update phi0 by the angular distance of one petal
    phi0 += 2.0 * M_PI / npetals; 
    
    // odd segments forward
    this->create_segmented_disk_layers(idisk, nsegments, false, phi0,  senZPos_odd_front );
    
    // odd segments backward
    this->create_segmented_disk_layers(idisk, nsegments, false, phi0, -senZPos_odd_front );
    
   
    
    
    // make the air disks
    
    TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air") ;
    
    Bool_t dummy  = false ;
    
    // place air discs to help transport during track extrapolation
    if( idisk < _nDisks-1 ){
      
      // place the disc half way between the two discs 
      double z = ( _FTDgeo[idisk].senZPos_even_front + _FTDgeo[idisk+1].senZPos_even_front ) * 0.5 ;
      
      TVector3 xc_fwd(0.0, 0.0, z) ;
      TVector3 normal_fwd(xc_fwd) ;
      normal_fwd.SetMag(1.0) ;
      
      double eps1 = 1.0e-04 ; // offset for disk number 
      double eps2 = 1.0e-05 ; // odd or even 
      double eps4 = 1.0e-08 ; // offset for forwards and backwards
      
      double height = _FTDgeo[idisk].height ;
      double rInner = _FTDgeo[idisk].rInner ;
      
      
      double sort_policy = rInner+height + eps1 * idisk + eps2 * 2 ; // we need to be after the even and odd
      
      Add(new ILDDiscMeasLayer( air, air, xc_fwd, normal_fwd, _bZ, sort_policy,
                                rInner, rInner+height, dummy,-1, "FTDAirSupportDiscFront" ) );
      
      TVector3 xc_bwd(0.0, 0.0, -z) ;
      TVector3 normal_bwd(xc_bwd) ;
      normal_bwd.SetMag(1.0) ;
      
      // offset needed for rear disks 
      sort_policy += eps4 ;
      
      
      Add(new ILDDiscMeasLayer( air, air, xc_bwd, normal_bwd, _bZ, sort_policy,
                                rInner, rInner+height, dummy,-1, "FTDAirSupportDiscRear" ) );
      
      
      
    }
    
  }
  
}






void ILDFTDKalDetector::create_segmented_disk_layers( int idisk, int nsegments, bool even_petals, double phi0, double zpos ){
  
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
  bool isDoubleSided = _FTDgeo[idisk].isDoubleSided ;
  int nSensors = _FTDgeo[idisk].nSensors ;
  int zsign = zpos > 0 ? +1 : -1 ;
  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  encoder.reset() ;  // reset to 0
  
  encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::FTD ;
  encoder[lcio::ILDCellID0::side] = zsign ;
  encoder[lcio::ILDCellID0::layer]  = idisk ;
  
  int start_index = even_petals ? 0 : 1 ;
  std::vector<int> sensors_front;
  std::vector<int> sensors_back;
  std::vector<int> module_ids_front;
  std::vector<int> module_ids_back;
  
  if( isDoubleSided ){ // sensors on front and back    double supZPos_odd_front = _FTDgeo[idisk].supZPos_odd;
    
    // first half is on the front, second half on the back, sensors start with sensor number 1
    for( int iSensor=1; iSensor <= nSensors/2; iSensor++ ) sensors_front.push_back( iSensor );
    for( int iSensor=nSensors/2 + 1; iSensor <= nSensors; iSensor++ ) sensors_back.push_back( iSensor );
    
  }
  else{ // only sensors on the front
    
    for( int iSensor=1; iSensor <= nSensors; iSensor++ ) sensors_front.push_back( iSensor );
    
  }

  for (int i=0; i<nsegments; ++i) {
    
    encoder[lcio::ILDCellID0::module] = start_index + i*2 ;
    
    for( unsigned j=0; j<sensors_front.size(); j++ ){
      
      encoder[lcio::ILDCellID0::sensor] = sensors_front[j];
      module_ids_front.push_back( encoder.lowWord() );
      
    }
    
    for( unsigned j=0; j<sensors_back.size(); j++ ){
      
      encoder[lcio::ILDCellID0::sensor] = sensors_back[j];
      module_ids_back.push_back( encoder.lowWord() );
      
    }
    
  }
  
  // create segmented disk 
  
  // front face of sensitive  
  double z = zpos - zsign*0.5*(senThickness) ;  
  //  double sort_policy = fabs(z) ;
  double eps1 = 1.0e-04 ; // disk  
  double eps2 = 1.0e-05 ; // odd or even 
  double eps3 = 1.0e-06 ; // layer in disk 
  double eps4 = 1.0e-08 ; // forward or backwards
  
  double sort_policy = rInner+height + eps1 * idisk + eps3 * 1 ;
  if ( ! even_petals ) sort_policy += eps2;
  
  // if this is the negative z disk add epsilon to the policy
  if( z < 0 ) sort_policy += eps4 ; 
  const char *name1 = z > 0 ? "FTDSenFrontPositiveZ" : "FTDSenFrontNegativeZ";  
  streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add front face of sensitive at " << z << std::endl;
  Add( new ILDSegmentedDiscMeasLayer(air, silicon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, dummy,name1) );
  
  
  // measurement plane

  z += zsign*0.5*senThickness;  
  //sort_policy = fabs(z) ;
  sort_policy = rInner+height + eps1 * idisk + eps3 * 2 ;
  if( z < 0 ) sort_policy += eps4 ;
  if ( ! even_petals ) { 
    sort_policy += eps2;
 
    const char *name2 = z > 0 ? "FTDMeasLayerFrontPositiveZOdd" : "FTDMeasLayerFrontNegativeZOdd";
    streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add measurement plane at " << z << " number of module_ids = " << module_ids_front.size() << std::endl;
    Add( new ILDSegmentedDiscMeasLayer(silicon, silicon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, active, module_ids_front,name2));
  }
  else{
    const char *name2 = z > 0 ? "FTDMeasLayerFrontPositiveZEven" : "FTDMeasLayerFrontNegativeZEven";
    streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add measurement plane at " << z << " number of module_ids = " << module_ids_front.size() << std::endl;
    Add( new ILDSegmentedDiscMeasLayer(silicon, silicon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, active, module_ids_front,name2));
  }
  
  // interface between sensitive and support
  z += zsign*0.5*senThickness;   
  //  sort_policy = fabs(z) ;
  sort_policy = rInner+height + eps1 * idisk + eps3 * 3 ;
  if( z < 0 ) sort_policy += eps4 ;
  if ( ! even_petals ) sort_policy += eps2;
  
  const char *name3 = z > 0 ? "FTDSenSupportIntfPositiveZ" : "FTDSenSupportIntfNegativeZ";
  
  streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add interface between sensitive and support at " << z << std::endl;
  Add( new ILDSegmentedDiscMeasLayer(silicon, carbon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, dummy,name3));
  
  if( isDoubleSided ){
    
    // interface between support and sensitive
    z += zsign*supThickness;   
    //  sort_policy = fabs(z) ;
    sort_policy = rInner+height + eps1 * idisk + eps3 * 4 ;
    if( z < 0 ) sort_policy += eps4 ;
    if ( ! even_petals ) sort_policy += eps2;
    
    const char *name4 = z > 0 ? "FTDSupportSenIntfPositiveZ" : "FTDSupportSenIntfNegativeZ";
    
    streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add interface between support and sensitive at " << z << std::endl;
    Add( new ILDSegmentedDiscMeasLayer(carbon, silicon , _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, dummy,name4));
    
    
    // measurement plane at the back
    z += zsign + 0.5*senThickness;   
    //  sort_policy = fabs(z) ;
    sort_policy = rInner+height + eps1 * idisk + eps3 * 5 ;
    if( z < 0 ) sort_policy += eps4 ;
    if ( ! even_petals ){ 
      
      sort_policy += eps2;
      
      const char *name5 = z > 0 ? "FTDMeasLayerBackPositiveZOdd" : "FTDMeasLayerBackNegativeZOdd";
      streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add measurement plane at " << z << " number of module_ids = " << module_ids_back.size() << std::endl;
      Add( new ILDSegmentedDiscMeasLayer(silicon, silicon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, active, module_ids_back,name5));
    }
    else{
      const char *name5 = z > 0 ? "FTDMeasLayerBackPositiveZEven" : "FTDMeasLayerBackNegativeZEven";
      streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add measurement plane at " << z << " number of module_ids = " << module_ids_back.size() << std::endl;
      Add( new ILDSegmentedDiscMeasLayer(silicon, silicon, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, active, module_ids_back,name5));
    }
    
    // rear face of sensitive
    z += zsign + 0.5*senThickness;  
    //  sort_policy = fabs(z) ;
    sort_policy = rInner+height + eps1 * idisk + eps3 * 4 ;
    if( z < 0 ) sort_policy += eps4 ;
    if ( ! even_petals ) sort_policy += eps2;
    
    const char *name6 = z > 0 ? "FTDSenRearPositiveZ" : "FTDSenRearNegativeZ";
    
    streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add rear face of sensitive at " << z << std::endl;
    Add( new ILDSegmentedDiscMeasLayer(silicon, air, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, dummy,name6));
    
    
    
  }
  else{
    
    // rear face of support
    z += zsign*supThickness;   
    //  sort_policy = fabs(z) ;
    sort_policy = rInner+height + eps1 * idisk + eps3 * 4 ;
    if( z < 0 ) sort_policy += eps4 ;
    if ( ! even_petals ) sort_policy += eps2;
    
    const char *name4 = z > 0 ? "FTDSupRearPositiveZ" : "FTDSupRearNegativeZ";
    
    streamlog_out(DEBUG) << "ILDFTDKalDetector::create_segmented_disk_layers add rear face of support at " << z << std::endl;
    Add( new ILDSegmentedDiscMeasLayer(carbon, air, _bZ, sort_policy, nsegments, z, phi0, rInner, height, innerBaseLength, outerBaseLength, dummy,name4));
    
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
    
    _FTDgeo[disk].senZPos_even_front = ftdlayers.getSensitiveZposition(disk, 0, 1) ;
    _FTDgeo[disk].senZPos_odd_front = ftdlayers.getSensitiveZposition(disk, 1, 1) ;
    
    _FTDgeo[disk].isDoubleSided = ftdlayers.isDoubleSided( disk );
    _FTDgeo[disk].nSensors = ftdlayers.getNSensors( disk );
    
    
  }
  
  
}
