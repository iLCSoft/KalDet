#ifndef ILDCYLINDERMEASLAYER_H
#define ILDCYLINDERMEASLAYER_H

/** ILDCylinderMeasLayer: User defined KalTest measurement layer class 
 *
 * @author S.Aplin DESY
 */


#include "ILDVMeasLayer.h"

class ILDCylinderMeasLayer : public ILDVMeasLayer, public TCylinder {
  
public:
  
  /** Constructor Taking inner and outer materials, radius and half length, B-Field, whether the layer is sensitive, Cell ID, and an optional name */
  ILDCylinderMeasLayer(TMaterial &min,
                       TMaterial &mout,
                       Double_t   r0,
                       Double_t   lhalf,
                       Double_t   Bz,
                       Bool_t     is_active,
                       Int_t      CellID = -1,
                       const Char_t    *name = "ILDCylinderMeasL") 
  : ILDVMeasLayer(min, mout, Bz, is_active, CellID, name),
  TCylinder(r0, lhalf)
  { /* no op */ }
  
  

  // Parent's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVector3   &xv)   const;
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv)   const 
  
  { return this->XvToMv(xv); }  


  /** Local to Global coordinates */
  virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;
  
  
  /** Calculate Projector Matrix */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)    const;
  
  /** Convert LCIO Tracker Hit to an ILDCylinderHit  */
  virtual ILDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const ;
  
private:
  
};
#endif
