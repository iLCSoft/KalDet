
#include "ILDSupportKalDetector.h"
#include "ILDCylinderMeasLayer.h"
#include "ILDPolygonBarrelMeasLayer.h"
#include "ILDDiscMeasLayer.h"

#include "TMath.h"
#include "TTUBE.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "MaterialDataBase.h"

#include <sstream>
#include <iomanip>
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
  
  const Double_t rtub       = pBeamPipe.getDoubleVal("BeamPipeRadius") ; // inner r of the tube
  const Double_t thickness  = pBeamPipe.getDoubleVal("BeamPipeThickness") ;   // thickness of beampipe
  const Double_t halfZ      = pBeamPipe.getDoubleVal("BeamPipeHalfZ") ;
  
  TMaterial & beam      = *MaterialDataBase::Instance().getMaterial("beam");
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & aluminium = *MaterialDataBase::Instance().getMaterial("aluminium");
  
  Bool_t dummy  = false;
  
  std::string name = "Support";
  
  // add beam pipe
  Add( new ILDCylinderMeasLayer(beam, aluminium , rtub, halfZ, bz, dummy,-1,"BeamPipeInr" ) );
  streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: [ beampipe ] at R = " << rtub
  << " X0_in = " << air.GetRadLength() << "  X0_out = " <<  aluminium.GetRadLength()    
  << std::endl ;  
  
  Add( new ILDCylinderMeasLayer(aluminium , air, rtub+thickness, halfZ, bz, dummy,-1,"BeamPipeOtr"  ) );
  streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: [ beampipe ] at R = " << rtub+thickness
  << " X0_in = " << aluminium.GetRadLength() << "  X0_out = " <<  air.GetRadLength()    
  << std::endl ;  
  
  
  // add vacuum layer 1mm inside the beam pipe to assist propagation to the IP
  const Double_t rvacuum = rtub - 1.0; 
  
  _ipLayer = new ILDCylinderMeasLayer(beam, beam , rvacuum , halfZ, bz, dummy,-1,"IPLayer" );
  
  Add( _ipLayer );
  streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: [ beampipe ] at R = " << rtub
  << " X0_in = " << air.GetRadLength() << "  X0_out = " <<  aluminium.GetRadLength()    
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
    double width = 2*(r_min_ecal_bar*tan(segment_dphi*0.5));
    Add ( new ILDParallelPlanarMeasLayer(air,air,r_min_ecal_bar,phi,bz,r_min_ecal_bar+i*1.0e-06,width,z_max_ecal_bar,0.0,true,encoder.lowWord(),"ECalBarrelFace"));    
    
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

