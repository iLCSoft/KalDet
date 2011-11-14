
#include "ILDCylinderHit.h"
#include "ILDCylinderMeasLayer.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;
using std::ios;
using std::resetiosflags;

//_________________________________________________________________________
//  --------------------------------
//  Implementation of public methods
//  --------------------------------
//


/** Global to Local coordinates */

TKalMatrix ILDCylinderHit::XvToMv(const TVector3 &xv, Double_t t0) const
{
  
  // SJA:FIXME the dynamic_cast should be checked in the constructor and then a static cast could be used 
  
  
  const ILDCylinderMeasLayer &ms
  = dynamic_cast<const ILDCylinderMeasLayer &>(this->GetMeasLayer());
  
  
  TKalMatrix h    = ms.XvToMv(*this, xv); 
  Double_t   r    = ms.GetR();
  Double_t   phih = h(0, 0) / r;
  Double_t   phim = (*this)(0, 0) / r;
  Double_t   dphi = phih - phim;
  
  static Double_t kPi    = TMath::Pi();
  static Double_t kTwoPi = 2 * kPi;
  
  while (dphi < -kPi) dphi += kTwoPi;
  while (dphi >  kPi) dphi -= kTwoPi;
  
  h(0, 0)  = r * (phim + dphi);
  
  return h;
}

/** Print Debug information */

void ILDCylinderHit::DebugPrint(Option_t *) const
{
  cerr << "------------------- Site Info -------------------------" << endl;
  
  for (Int_t i = 0; i < GetDimension(); i++) {
    Double_t x  = (*this)(i, 0);
    Double_t dx = (*this)(i, 1);
    cerr << " x[" << i << "] = " << setw(8) << setprecision(5) << x
    << "    "
    << "dx[" << i << "] = " << setw(6) << setprecision(2) << dx
    << setprecision(7)
    << resetiosflags(ios::showpoint)
    << endl;
  }
  cerr << " r of ILDCylinderMeasLayer = " << setw(8)
  << static_cast<const ILDCylinderMeasLayer &>(GetMeasLayer()).GetR() << endl;
  cerr << "-------------------------------------------------------"  << endl;
}

