
#include "ILDMeasurementSurfaceStoreFiller.h"

#include "UTIL/ILDConf.h"

#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>
#include <gear/FTDLayerLayout.h>
#include <gear/FTDParameters.h>
#include <gear/GEAR.h>
#include <gear/GearMgr.h>

#include "streamlog/streamlog.h"


#include "gear/gearsurf/CartesianCoordinateSystem.h"
#include "gear/gearsurf/MeasurementSurface.h"
#include "gear/gearsurf/BoundaryRectangle.h"
#include "gear/gearsurf/BoundaryTrapezoid.h"

void ILDMeasurementSurfaceStoreFiller::get_gear_parameters(const gear::GearMgr& gear_mgr) {

  _paramVXD = &(gear_mgr.getVXDParameters());
  _paramSIT = &(gear_mgr.getSITParameters());
  _paramSET = &(gear_mgr.getSETParameters());
  _paramFTD = &(gear_mgr.getFTDParameters());
  
  // set the strip angle values, that are at the moment not in gear, so it's hardcoded here
#ifdef HARDCODEDGEAR
  
  unsigned nVTXLayers = _paramVXD->getZPlanarLayerLayout().getNLayers();
  for( unsigned i=0; i<nVTXLayers; i++) _VTXStripAngles.push_back( 0. );
  
  unsigned nSITLayers = _paramSIT->getZPlanarLayerLayout().getNLayers();
  //  for( unsigned i=0; i<nSITLayers; i++) _SITStripAngles.push_back( pow(-1,i) * 5.* M_PI/180. ); // alternately + and - 5°

  for( unsigned i=0; i<nSITLayers; i++) _SITStripAngles.push_back( 0.0 ); 
  
  unsigned nSETLayers = _paramSET->getZPlanarLayerLayout().getNLayers();
  //  for( unsigned i=0; i<nSETLayers; i++) _SETStripAngles.push_back( pow(-1,i) * 5.* M_PI/180. ); // alternately + and - 5°
  for( unsigned i=0; i<nSETLayers; i++) _SETStripAngles.push_back( 0.0 );
  
  const gear::FTDLayerLayout& ftdLayers = _paramFTD->getFTDLayerLayout() ;
  unsigned nFTDLayers = ftdLayers.getNLayers();

  for( unsigned layer=0; layer<nFTDLayers; layer++ ){
    
    std::vector< double > angles;
    
    unsigned nSensors = 4; //ftdLayers.getNSensors( layer );
    bool isDoubleSided = true; //ftdLayers.isDoubleSided( layer );
    
    if( isDoubleSided){ // it is a strip detector
      
      
      for( unsigned sensor=1; sensor <= nSensors; sensor++ ){
        
//        if ( sensor <= nSensors/2 ) angles.push_back( 5.* M_PI/180. );   // the first half of the sensors is in front with one angle,
//        else                        angles.push_back( -5.* M_PI/180. );  // the other is in the back with the opposite angle 
        
        angles.push_back(0.0);
        
      }
      
    }
    else{ // a pixel detector, so no stripAngle
      
      for( unsigned sensor=1; sensor <= nSensors; sensor++ ) angles.push_back(0.);
      
    }
    
    _FTDStripAngles.push_back( angles );
    
  }

  
}
#endif



void ILDMeasurementSurfaceStoreFiller::getMeasurementSurfaces( std::vector<MeasurementSurface*>& surface_list ) const {
  
  
  this->storeZPlanar( _paramVXD , UTIL::ILDDetID::VXD, surface_list );
  this->storeZPlanar( _paramSIT , UTIL::ILDDetID::SIT, surface_list );
  this->storeZPlanar( _paramSET , UTIL::ILDDetID::SET, surface_list );
  this->storeFTD( _paramFTD , surface_list);
  
  
}


void ILDMeasurementSurfaceStoreFiller::storeZPlanar( const gear::ZPlanarParameters* param , int det_id, std::vector<MeasurementSurface*>& surface_list  ) const {
  
  
#ifdef HARDCODEDGEAR
  std::vector< double > angles;
  if( det_id == UTIL::ILDDetID::VXD ) angles = _VTXStripAngles;
  else if ( det_id == UTIL::ILDDetID::SIT ) angles = _SITStripAngles;
  else if ( det_id == UTIL::ILDDetID::SET ) angles = _SETStripAngles;
  else return;
#endif
  
  const gear::ZPlanarLayerLayout& layerLayout = param->getZPlanarLayerLayout();
  
  unsigned nLayers = layerLayout.getNLayers();
  
  for( unsigned layerNumber = 0; layerNumber < nLayers; layerNumber++ ){
    
    unsigned nLadders = layerLayout.getNLadders( layerNumber );
    double ladder_r            = layerLayout.getSensitiveDistance(layerNumber); // the distance of the ladders from (0,0,0)
    double sensitive_offset    = layerLayout.getSensitiveOffset(layerNumber); // the offset, see ZPlanarLayerLayout.h for more details
    double deltaPhi            = ( 2 * M_PI ) / nLadders ; // the phi difference between two ladders
    double phi0                = layerLayout.getPhi0( layerNumber );
    double sensitive_length  = layerLayout.getSensitiveLength(layerNumber) * 2.0 ; // note: gear for historical reasons uses the halflength 
    double sensitive_width  = layerLayout.getSensitiveWidth(layerNumber);
    double sensitive_thickness = layerLayout.getSensitiveThickness(layerNumber);
    
    double stripAngle = 0.; 
    // TODO: this should come from gear or somewhere else! And not be hardcoded
#ifdef HARDCODEDGEAR
    stripAngle = angles[layerNumber];       
#endif
    
    for( unsigned ladderNumber = 0; ladderNumber < nLadders; ladderNumber++ ){
      
      // determine the CellID0 for the  ladder
      UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
      cellID[ lcio::ILDCellID0::subdet ] = det_id ;
      cellID[ lcio::ILDCellID0::side   ] = 0 ;
      cellID[ lcio::ILDCellID0::layer  ] = layerNumber ;
      cellID[ lcio::ILDCellID0::module ] = ladderNumber ;
      cellID[ lcio::ILDCellID0::sensor ] = 0 ;
      int cellID0 = cellID.lowWord();
    
      
      // Let's start with the translation T: the new center of coordinates:
      // The center of the first ladder (when we ignore an offset and phi0 for now) is (R,0,0)
      // If we include the offset, the center gets shifted to (R,offset,0)
      CLHEP::Hep3Vector T( ladder_r + sensitive_thickness/2., sensitive_offset, 0 );
      // Now we have to take into account phi0 and that the number of the ladder.
      // Together the center is rotated by phi0 + ladderNumber*deltaPhi around the z axis
      CLHEP::HepRotation rot;
      rot.rotateZ( deltaPhi * ladderNumber + phi0 );
      
      T = rot * T;
      
      // Next, we want to determinte the rotation matrix R
      // We start with u,v,w alligned with x,y,z.
      // As u is perpendicular to the strip orientation it looks like this.
      //               y
      //               |     
      //            ---|---
      //            |  |  |
      //            |  |  |
      //      <--------|--------> x
      //            |  |  |
      //            |  |  |
      //            ---|---
      //               |
      // To get this sensor in place we have to do a few rotations:
      // First we'll rotate around the z axis. With a strip angle of 0,
      // we would just rotate by 90°, but with a strip angle by
      // 90°-stripAngle in clockwise direction. 
      CLHEP::HepRotation R;
      //           printRotation( R );        
      R.rotateZ( stripAngle - M_PI/2. );
      //           printRotation( R );
      
      // Next we rotate 90° clockwise around y, so the strip now points in z direction (if strip angle == 0)
      R.rotateY( -M_PI/2. );
      
      // Finally we have to get the ladder in place w.r. to its number and the resulting phi angle
      R.rotateZ( deltaPhi * ladderNumber + phi0 );
      
      CartesianCoordinateSystem* cartesian = new CartesianCoordinateSystem( T, R );
      
      BoundaryRectangle* b = new BoundaryRectangle( sensitive_width, sensitive_length, 1., -stripAngle );
      MeasurementSurface* ms = new MeasurementSurface( cellID0, cartesian, b );
      surface_list.push_back(ms);
      
    }
    
  }
  
}


void ILDMeasurementSurfaceStoreFiller::storeFTD( const gear::FTDParameters* param, std::vector<MeasurementSurface*>& surface_list ) const {
  
  const gear::FTDLayerLayout& ftdLayers = param->getFTDLayerLayout() ;
  unsigned nLayers = ftdLayers.getNLayers();
  
  UTIL::BitField64  cellID( UTIL::ILDCellID0::encoder_string );
  cellID[ lcio::ILDCellID0::subdet ] = UTIL::ILDDetID::FTD ;
  
  for( unsigned layer = 0; layer < nLayers; layer++ ){
    
    
    cellID[ lcio::ILDCellID0::layer  ] = layer ;
    
    unsigned nPetals = ftdLayers.getNPetals( layer );
    double deltaPhi  = ( 2 * M_PI ) / nPetals;
    double phi0      = ftdLayers.getPhi0( layer );
    double lengthMin = ftdLayers.getSensitiveLengthMin( layer );
    double lengthMax = ftdLayers.getSensitiveLengthMax( layer );
    double width     = ftdLayers.getSensitiveWidth( layer );
    unsigned nSensors = ftdLayers.getNSensors( layer ); 
    bool isDoubleSided = ftdLayers.isDoubleSided( layer );
    // Get some information about the petal (a trapezoid, see gear for in depth description)
    double distMin = ftdLayers.getSensitiveRinner(layer); // the xy-distance from (0,0) of the center of the inner base
    double distMax = distMin + width; // the xy-distance from (0,0) of the center of the outer base
    unsigned nSensorsOn1Side = nSensors;
    if( isDoubleSided ) nSensorsOn1Side = nSensorsOn1Side / 2;
    
    for( unsigned petal=0; petal< nPetals; petal++ ){
      
      
      cellID[ lcio::ILDCellID0::module ] = petal ;
      
      //TODO: should on back to back strip detectors of the FTD the w-vectors point in opposite directions or the same?
      for ( unsigned sensor = 1; sensor <=nSensors; sensor++ ){
        
        streamlog_out(DEBUG1) << "layer = " << layer << "\tpetal = " << petal << "\tsensor = " << sensor << "\n";
        
        double stripAngle = 0.; 
        
        // TODO: this should come from gear and not be hardcoded
#ifdef HARDCODEDGEAR
        stripAngle = _FTDStripAngles[layer][sensor];       
#endif
        
        cellID[ lcio::ILDCellID0::side   ] = -1 ;                    
        cellID[ lcio::ILDCellID0::sensor ] = sensor ;
        int cellID0 = cellID.lowWord();
        
        // We start with the Translation Vector T (=origin of new CoordinateSystem )
        // the first petal (if phi0 = 0) sits with its symmetry axis on the x axis.
        // It starts at x= rInner (lMin)nand reaches rInner+ width (lMax) on the x axis
        // So far the center of Mass will be on the x axis.
        //
        // Now we have to take into account how many sensors there are on the petal, and where they sit.
        double deltaX = (distMax-distMin) / double( nSensorsOn1Side ); //how much delta x one sensor covers
        
        // The first sensor is the outermost one. Let's calculate its center
        double xFirst = distMax - deltaX/2.;
        
        // From here, all we have to do is go down in steps of deltaX until we reach the right sensor.
        // We have to be careful: when there are sensors on the back, they start again at the top.
        // If we for example have 6 sensors on a doublesided petal: on the front it is from out to inside: 1,2,3
        // and on the back: 4,5,6. So 1 and 4 will be (apart from z) at the same position. So will 2 and 5, as well
        // as 3 and 6
        unsigned steps = (sensor-1)%nSensorsOn1Side; // the number of steps to go deltaX (in negative direction)
                                                     // In our example with 6 sensors, this gives sensor 1: 0 steps, 2:1, 3:2, 4:0, 5:1, 6:2
        
        double x = xFirst - steps * deltaX;
        double y = 0.;
        double z = -ftdLayers.getSensitiveZposition( layer, petal, sensor ); // we first calculate for the petals at -z 
                                                                             //(because the imagination of the rotation is easier, as w point the same direction as z)
                                                                             // later we from the result we can calculate R and T for +z pretty easy.
        
        // Calculate lengths of the inner and outer base of the trapezoid
        double baseOuter = lengthMax - steps * (lengthMax - lengthMin);
        double baseInner = baseOuter - (lengthMax - lengthMin) / nSensorsOn1Side;
        
        CLHEP::Hep3Vector T( x , y, z);
        //             printVector( T );
        // Now we only have to rotate the petal around the z-axis into its place
        CLHEP::HepRotation rot;
        rot.rotateZ( petal * deltaPhi + phi0 );
        T = rot * T;
        //             printVector( T );
        
        //On to the rotation matrix
        CLHEP::HepRotation R;
        //             printRotation( R );
        R.rotateZ( petal * deltaPhi + phi0 + stripAngle - M_PI/2. );
        //             printRotation( R );

        BoundaryTrapezoid* b1 = new BoundaryTrapezoid( baseInner, baseOuter, deltaX, 1., -stripAngle);
        CartesianCoordinateSystem* cartesian = new CartesianCoordinateSystem( T, R );
        MeasurementSurface* ms = new MeasurementSurface( cellID0, cartesian, b1);
        surface_list.push_back(ms);
        
        
        // Once more for the other side
        cellID[ lcio::ILDCellID0::side ] = 1 ;   
        cellID0 = cellID.lowWord();
        
        T.setZ( -T.z() ); // switch to -z
        
        // R is pretty much the same as the strip orientation will be the same,
        // but as (by chosen definition) w should point towards the IP,
        // we have to flip u around. So we acutally have to rotate 180° around v
        // So first we get the vector v
        CLHEP::Hep3Vector v = R*CLHEP::Hep3Vector(0,1,0);
        // Then we rotate around it
        R.rotate( M_PI , v );
        //             printRotation( R );

        BoundaryTrapezoid* b2 = new BoundaryTrapezoid( baseInner, baseOuter, deltaX, 1., stripAngle);
        CartesianCoordinateSystem* cartesian2 = new CartesianCoordinateSystem( T, R );
        MeasurementSurface* ms2 = new MeasurementSurface( cellID0, cartesian2, b2);
        surface_list.push_back(ms2);
        
      
        
        //TODO: do we need alpha in here?
        
        
      }
      
    }      
    
  }
  
}







