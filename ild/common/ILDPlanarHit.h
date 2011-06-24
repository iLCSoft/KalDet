#ifndef ILDPLANARHIT_H
#define ILDPLANARHIT_H

#include "KalTrackDim.h"

#include "ILDVTrackHit.h"

#include "ILDPlanarMeasLayer.h"


class ILDPlanarHit : public ILDVTrackHit {

public:
  
  
  ILDPlanarHit(const TVMeasLayer &ms,
	       Double_t       *x,
	       Double_t       *dx,
	       Double_t        bfield) ;

  
  virtual ~ILDPlanarHit();
  
  virtual TKalMatrix XvToMv (const TVector3 &xv, Double_t t0) const;
  virtual void       DebugPrint(Option_t *opt = "")           const;
  
  
private:


};
#endif
