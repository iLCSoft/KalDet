
#include "ILDPlanarHit.h"
#include "ILDPlanarMeasLayer.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

using namespace std;


//_____________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------

      
ILDPlanarHit::ILDPlanarHit(const TVMeasLayer &ms,
			   Double_t       *x,
			   Double_t       *dx, 
			   Double_t        bfield)
        : ILDVTrackHit(ms, x, dx, bfield, 2)
{
}

ILDPlanarHit::~ILDPlanarHit()
{
}

//_____________________________________________________________________
//  ----------------------------------
//  Implementation of public methods  
//  ----------------------------------

TKalMatrix ILDPlanarHit::XvToMv(const TVector3 &xv, Double_t /*t0*/) const
{
   return dynamic_cast<const ILDPlanarMeasLayer &>(GetMeasLayer()).XvToMv(xv);
}

void ILDPlanarHit::DebugPrint(Option_t *) const
{
   cerr << "------------------- Site Info -------------------------" << endl;

   for (Int_t i=0; i<GetDimension(); i++) {
      Double_t x  = (*this)(i,0);
      Double_t dx = (*this)(i,1);
      cerr << " x[" << i << "] = " << setw(8) << setprecision(5) << x
           << "    "
           << "dx[" << i << "] = " << setw(6) << setprecision(2) << dx
           << setprecision(7)
           << resetiosflags(ios::showpoint)
           << endl;
   }
   cerr << "-------------------------------------------------------" << endl;
}
