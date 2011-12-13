#ifndef ILDPLANARHIT_H
#define ILDPLANARHIT_H

/** ILDPlanarStripHit: User defined KalTest hit class using u coordinate, which provides coordinate vector as defined by the MeasLayer 
 *  
 * @author S.Aplin DESY
 */

#include "KalTrackDim.h"

#include "ILDVTrackHit.h"


class ILDPlanarStripHit : public ILDVTrackHit {
  
public:
  
  /** Constructor Taking u and v coordinates and associated measurement layer, with bfield */
  ILDPlanarStripHit(const TVMeasLayer &ms,
               Double_t       *x,
               Double_t       *dx,
               Double_t        bfield,
               EVENT::TrackerHit* trkhit) 
  : ILDVTrackHit(ms, x, dx, bfield, 1,trkhit)
  { /* no op */ } 
  
  // TVTrackHit's pure virtuals that must be implemented
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;
  
  /** Print Debug information */
  virtual void       DebugPrint(Option_t *opt = "")           const;
  
  
private:
  
  
};
#endif
