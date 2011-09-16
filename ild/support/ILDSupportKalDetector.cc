
#include "ILDSupportKalDetector.h"
#include "ILDCylinderMeasLayer.h"

#include "TMath.h"
#include "TTUBE.h"

#include "MaterialDataBase.h"

#include <sstream>
#include <iomanip>

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gearimpl/Util.h"

#include "streamlog/streamlog.h"


ILDSupportKalDetector::ILDSupportKalDetector( const gear::GearMgr& gearMgr ) : 
  TVKalDetector(10) 
{
  
  streamlog_out(DEBUG1) << "ILDSupportKalDetector building beampipe using GEAR " << std::endl ;

  const gear::GearParameters& pBeamPipe = gearMgr.getGearParameters("BeamPipe");
  const Double_t bz = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;

  static const Double_t rtub       = pBeamPipe.getDoubleVal("BeamPipeRadius") ; // inner r of the tube
  static const Double_t thickness  = pBeamPipe.getDoubleVal("BeamPipeThickness") ;   // thickness of beampipe
  static const Double_t halfZ      = pBeamPipe.getDoubleVal("BeamPipeHalfZ") ;

  TMaterial & beam      = *MaterialDataBase::Instance().getMaterial("beam");
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & aluminium = *MaterialDataBase::Instance().getMaterial("aluminium");
  
  Bool_t dummy  = false;
   
  std::string name = "Support";

 // add beam pipe
  Add( new ILDCylinderMeasLayer(beam, aluminium , rtub, halfZ, bz, dummy ) );
  streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using layerID: [ beampipe ] at R = " << rtub
			    << " X0_in = " << air.GetRadLength() << "  X0_out = " <<  aluminium.GetRadLength()    
			    << std::endl ;  

  Add( new ILDCylinderMeasLayer(aluminium , air, rtub+thickness, halfZ, bz, dummy ) );
  streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using layerID: [ beampipe ] at R = " << rtub+thickness
			    << " X0_in = " << aluminium.GetRadLength() << "  X0_out = " <<  air.GetRadLength()    
			    << std::endl ;  

  
  SetOwner();
}

ILDSupportKalDetector::~ILDSupportKalDetector()
{
}

