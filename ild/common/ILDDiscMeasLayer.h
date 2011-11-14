#ifndef __ILDDISCMEASLAYER__
#define __ILDDISCMEASLAYER__

/** ILDDiscMeasLayer: User defined KalTest Disc measurement layer class used with ILDPLanarTrackHit. 
 *
 * @author S.Aplin DESY
 */

#include "TVector3.h"
#include "TKalMatrix.h"
#include "TPlane.h"
#include "ILDVMeasLayer.h"
#include "KalTrackDim.h"
#include "TMath.h"
#include <sstream>
class TVTrackHit;


class ILDDiscMeasLayer : public ILDVMeasLayer, public TPlane {
  
public:
  
  /** Constructor Taking inner and outer materials, center and normal to the plane, B-Field, sorting policy, min and max r, whether the layer is sensitive, Cell ID, and an optional name */
  
  ILDDiscMeasLayer(TMaterial &min,
                   TMaterial &mout,
                   const TVector3  &center,
                   const TVector3  &normal,
                   double   Bz,
                   double   SortingPolicy,
                   double   rMin,
                   double   rMax,
                   Bool_t     is_active,
                   Int_t      CellID = -1,
                   const Char_t    *name = "ILDDiscMeasL")
  : ILDVMeasLayer(min, mout, Bz, is_active, CellID, name),
  TPlane(center, normal),
  _sortingPolicy(SortingPolicy), _rMin(rMin), _rMax(rMax)
  { /* no op */ }
  
  
  
  // Parrent's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv) const 
  { return this->XvToMv(xv); }
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
  
  /** Local to Global coordinates */  
  virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
  
  /** Calculate Projector Matrix */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const;
  
  /** Convert LCIO Tracker Hit to an ILDPLanarTrackHit  */
  virtual ILDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const ;
  
  /** Check if global point is on surface  */
  inline virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  /** Get sorting policy for this plane  */
  double GetSortingPolicy() const { return _sortingPolicy; }
  
private:
  double _sortingPolicy;
  double _rMin;
  double _rMax;
  
};

#endif
