#ifndef ILDVTrackHIT_H
#define ILDVTrackHIT_H

/** ILDVMeasLayer:  Virtual hit class used by ILD[X]Hit Classes, which should provide coordinate vector as defined by the MeasLayer
 *
 * @author S.Aplin DESY
 */


#include "kaltest/TVTrackHit.h"

#include "ILDVMeasLayer.h"


class ILDVTrackHit : public TVTrackHit {
  
public:
  
   /** Constructor Taking coordinates and associated measurement layer, with bfield and number of measurement dimentions*/
  ILDVTrackHit(const TVMeasLayer &ms, Double_t *x, Double_t *dx, 
               Double_t bfield , Int_t dim) 
  : TVTrackHit(ms, x, dx, bfield, dim) 
  { /* no op */ }
  
  
private:
  
  
};
#endif
