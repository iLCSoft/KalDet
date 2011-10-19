//*************************************************************************
//* ================
//*  ILDCylinderHit Class
//* ================
//*
//* (Description)
//*   User defined hit class
//*   provides coordinate vector as defined by the MeasLayer
//* (Requires)
//*     TVTrackHit
//* (Provides)
//*     class ILDCylinderHit
//*
//*************************************************************************
//
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
//  --------------
//  Ctors and Dtor
//  --------------
//

ILDCylinderHit::ILDCylinderHit(const TVMeasLayer &ms, Double_t *x, Double_t *dx, 
                               Double_t bfield ) 
: ILDVTrackHit(ms, x, dx, bfield, 2) 
{
  
  //SJA:FIXME is there any way to check the size of the Double_t supplied?
  
}


ILDCylinderHit::~ILDCylinderHit()
{
}

//_________________________________________________________________________
//  --------------------------------
//  Implementation of public methods
//  --------------------------------
//
TKalMatrix ILDCylinderHit::XvToMv(const TVector3 &xv, Double_t t0) const
{
  const ILDCylinderMeasLayer &ms
  = dynamic_cast<const ILDCylinderMeasLayer &>(GetMeasLayer());
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

