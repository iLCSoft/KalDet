
#include "ILDVTrackHit.h"


ILDVTrackHit::ILDVTrackHit(const TVMeasLayer &ms, Double_t *x, Double_t *dx, 
                           Double_t bfield, Int_t dim)
: TVTrackHit(ms, x, dx, bfield, dim) 
{
}

ILDVTrackHit::~ILDVTrackHit()
{
}
