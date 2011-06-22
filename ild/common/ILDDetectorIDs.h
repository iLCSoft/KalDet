#ifndef INCLUDE_ILDDetectorIDs
#define INCLUDE_ILDDetectorIDs 1

#include <lcio.h>
#include <LCRTRelations.h>

namespace ILDDetectorIDs{

  // helper class to assign additional parameters to TrackerHits
  struct HitInfoStruct{
  HitInfoStruct() :layerID(-1) {}
    int layerID ;
  } ;
  struct HitInfo : lcio::LCOwnedExtension<HitInfo, HitInfoStruct> {} ;
  
  /** Enums for identifying detectors
   */
  struct DetID {
    static const int unknown = 0 ;
    static const int VXD =  1 ;
    static const int SIT =  2 ;
    static const int TPC =  3 ;
    static const int SET =  4 ;
    static const int ETD =  5 ;
    static const int FTD =  6 ; // add new detecors here
    static const int IP =   7 ; // add new detecors here
    static const int Size = 8 ;
    static const int Factor = 10000 ; 
  } ;
  

}


#endif
