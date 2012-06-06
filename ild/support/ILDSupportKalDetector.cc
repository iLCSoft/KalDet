
#include "ILDSupportKalDetector.h"
#include "ILDCylinderMeasLayer.h"
#include "ILDConeMeasLayer.h"
#include "ILDPolygonBarrelMeasLayer.h"
#include "ILDDiscMeasLayer.h"

#include "TMath.h"
#include "TTUBE.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "MaterialDataBase.h"

#include <sstream>
#include <cmath>

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gearimpl/Util.h"
#include "gear/CalorimeterParameters.h"
#include "streamlog/streamlog.h"


ILDSupportKalDetector::ILDSupportKalDetector( const gear::GearMgr& gearMgr ) : 
TVKalDetector(10) 
{
  
  streamlog_out(DEBUG1) << "ILDSupportKalDetector building beampipe using GEAR " << std::endl ;
  
  const gear::GearParameters& pBeamPipe = gearMgr.getGearParameters("BeamPipe");
  const Double_t bz = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  // the beampipe is supposed to be a chain of cut cones (cut, that means the spike is cut off, also called a cone frustum).
  // as they are connected, RStart[i] == REnd[i-1]. With this all we need are the z values and the radii at the place.
  const std::vector<double> z = pBeamPipe.getDoubleVals("Z");
  const std::vector<double> rInner = pBeamPipe.getDoubleVals("RInner"); //inner radius of the cone
  const std::vector<double> rOuter = pBeamPipe.getDoubleVals("ROuter"); //outer radius of the cone
  
  MaterialDataBase::Instance().registerForService(gearMgr);
  TMaterial & beam      = *MaterialDataBase::Instance().getMaterial("beam");
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  //  TMaterial & aluminium = *MaterialDataBase::Instance().getMaterial("aluminium");
  TMaterial & beryllium = *MaterialDataBase::Instance().getMaterial("beryllium");
  
  Bool_t dummy  = false;
  
  
  // add beam pipe cones
  for( unsigned i=0; i<z.size()-1; i++){
     
     double zStart = z[i];
     double zEnd = z[i+1];
     double rInnerStart = rInner[i];
     double rInnerEnd = rInner[i+1];
//     double rOuterStart = rOuter[i];
//     double rOuterEnd = rOuter[i+1];
     
     std::stringstream sname;
     sname << "BeamPipeCone" << i;
     std::string name = sname.str();
      
     double epsilon = 0.001;
     if( fabs( zEnd-zStart ) > epsilon ){
     
////         Add( new ILDConeMeasLayer(beam, beryllium , zStart, rInnerStart, zEnd, rInnerEnd, bz, dummy,-1, name.c_str() ) );
////         Add( new ILDConeMeasLayer(beam, beryllium , -zStart, rInnerStart, -zEnd, rInnerEnd, bz, dummy,-1, name.c_str() ) );
         streamlog_out( DEBUG0 )   << " *** adding inner " << name << " Measurement layer using CellID: [ beampipe ] at"
         << " z1 = +-" << zStart
         << " z2 = +-" << zEnd
         << " r1 = " << rInnerStart
         << " r2 = " << rInnerEnd 
         << " X0_in = " << beam.GetRadLength() << "  X0_out = " <<  beryllium.GetRadLength()    
         << std::endl ;  
            
         
////         Add( new ILDConeMeasLayer(beryllium , air , zStart, rOuterStart, zEnd, rOuterEnd, bz, dummy,-1, name.c_str() ) );
////         Add( new ILDConeMeasLayer(beryllium , air , -zStart, rOuterStart, -zEnd, rOuterEnd, bz, dummy,-1, name.c_str() ) );
         streamlog_out( DEBUG0 )   << " *** adding outer " << name << " Measurement layer using CellID: [ beampipe ] at"
         << " z1 = +-" << zStart
         << " z2 = +-" << zEnd
         << " r1 = " << rInnerStart
         << " r2 = " << rInnerEnd 
         << " X0_in = " << beryllium.GetRadLength() << "  X0_out = " <<  air.GetRadLength()    
         << std::endl ;  
    
     }
     
  }
   
  
  // add vacuum layer 1mm inside the beam pipe to assist propagation to the IP
  // therefore make a cylinder that is 1mm smaller than the lowest RInner value of the cones
  // and make it so long that it gets as long as the beamtube
  const double rvacuum = *min_element(rInner.begin(), rInner.end()) - 1.0; 
  const double zHalfVacuum = z.back();
  
  _ipLayer = new ILDCylinderMeasLayer(beam, beam , rvacuum , zHalfVacuum, bz, dummy,-1,"IPLayer" );
  
  Add( _ipLayer );
  streamlog_out( DEBUG0 )   << " *** adding " << "IPLayer" << " Measurement layer using CellID: [ beampipe ] at R = " << rvacuum
  << " zHalf = " << zHalfVacuum << " X0_in = " << beam.GetRadLength() << "  X0_out = " <<  beam.GetRadLength()    
  << std::endl ;  
  
  
  // add calo bounding inner surface as ILDPolygonBarrelMeasLayer and two planes at the inner z of the calo endcaps with ILDDiscMeasLayer

  // SJA:FIXME: OK for now we will just use a set of planes as there is not yet an implementation to propagate to a multilayer
  
  const gear::CalorimeterParameters& ecalB = gearMgr.getEcalBarrelParameters();
  const gear::CalorimeterParameters& ecalE = gearMgr.getEcalEndcapParameters();
  
  if (ecalB.getSymmetryOrder()!=8) {
    streamlog_out(ERROR) << "ILDSupportKalDetector::ILDSupportKalDetector ECal barrel is not eightfold symmetry: exit(1) called from " << __FILE__ << "   line " << __LINE__ << std::endl; 
    exit(1);

  }
  
  double phi0  = ecalB.getPhi0();
  double r_min_ecal_bar = ecalB.getExtent()[0];
  double z_max_ecal_bar = ecalB.getExtent()[3];

  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  encoder.reset() ;  // reset to 0
  
  encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::ECAL ;
  encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::barrel;
  encoder[lcio::ILDCellID0::layer]  = 0 ;
  
  std::vector<int> module_ids;
  
  for (int i=0; i<8; ++i) {

    encoder[lcio::ILDCellID0::module] = i;
    module_ids.push_back(encoder.lowWord());

    double segment_dphi = 2.0*M_PI / 8; 
    double phi = i*segment_dphi+phi0;
    double width = 2.0*(r_min_ecal_bar*tan(segment_dphi*0.5));
    double length = 2.0*z_max_ecal_bar;
    Add ( new ILDParallelPlanarMeasLayer(air,air,r_min_ecal_bar,phi,bz,r_min_ecal_bar+i*1.0e-06,width,length,0.0,0.0,0.0,true,encoder.lowWord(),"ECalBarrelFace"));    
    
  }
  
  //  Add ( new ILDPolygonBarrelMeasLayer(air,air,bz,r_min_ecal_bar,r_min_ecal_bar,z_max_ecal_bar*0.5,8,0.0,phi0,module_ids,"ECalBarrelFace"));  
  
  streamlog_out( DEBUG0 )   << " *** adding ECalBarrelFace Measurement layer at Rmin = " << r_min_ecal_bar << std::endl ;
  
  double r_max_ecal_ecap = ecalE.getExtent()[1];
  double z_min_ecal_ecap = ecalE.getExtent()[2];

  encoder[lcio::ILDCellID0::module] = 0;
  
  encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::fwd;
  
  TVector3 front_face_centre_fwd( 0.0, 0.0, z_min_ecal_ecap); // for +z  
  TVector3 front_face_normal_fwd(front_face_centre_fwd);
  front_face_normal_fwd.SetMag(1.0);
  
  Add (new ILDDiscMeasLayer(air, air, front_face_centre_fwd,front_face_normal_fwd, bz, z_min_ecal_ecap, 0., r_max_ecal_ecap/cos(M_PI/8.0), true, encoder.lowWord(),"ECalEndcapFace+Z") );
  
  streamlog_out( DEBUG0 )   << " *** adding ECalEndcapFace+Z Measurement layer at Zmin = " << z_min_ecal_ecap << " and Rmax = " << r_max_ecal_ecap/cos(M_PI/8.0) << std::endl ;
  
  encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::bwd;
  TVector3 front_face_centre_bwd( -front_face_centre_fwd ); // for -z  
  TVector3 front_face_normal_bwd(front_face_centre_bwd);
  front_face_normal_bwd.SetMag(1.0);
  
  Add (new ILDDiscMeasLayer(air, air, front_face_centre_bwd,front_face_normal_bwd, bz, z_min_ecal_ecap+1.0e-06, 0., r_max_ecal_ecap/cos(M_PI/8.0), true, encoder.lowWord(),"ECalEndcapFace-Z"));

  streamlog_out( DEBUG0 )   << " *** adding ECalEndcapFace-Z Measurement layer at Zmin = " << z_min_ecal_ecap << " and Rmax = " << r_max_ecal_ecap/cos(M_PI/8.0) << std::endl ;
  
  SetOwner();
}

