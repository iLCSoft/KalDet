//*************************************************************************
//* ================
//*  EXTPCHit Class
//* ================
//*
//* (Description)
//*   User defined hit class
//*   provides coordinate vector as defined by the MeasLayer
//* (Requires)
//*     TVTrackHit
//* (Provides)
//*     class EXTPCHit
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/tpc/EXTPCHit.cxx
//*   2009/11/23  K.Ikematsu   Added GetSortingPolicy() method
//*   2009/11/23  K.Ikematsu   Added data member *fHitPtr
//*                            as a pointer to raw Hit object
//*
//* $Id: EXTPCHit.cxx,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "EXTPCHit.h"
#include "EXTPCMeasLayer.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

using namespace std;

ClassImp(EXTPCHit)

//_________________________________________________________________________
//  --------------
//  Ctors and Dtor
//  --------------
//
EXTPCHit::EXTPCHit(Int_t m)
        : TVTrackHit(m),
          fSide(0),
          fVdrift(0)
{
}

EXTPCHit::EXTPCHit(const EXTPCMeasLayer &ms,
                         Double_t       *x,
                         Double_t       *dx,
                         Int_t           side,
                         Double_t        v,
                   const TVector3       &xx,
                         Double_t        b,
                         Int_t           m)
        : TVTrackHit(ms, x, dx, b, m),
          fSide(side),
          fVdrift(v),
          fXXPtr(&xx)
{
}

EXTPCHit::EXTPCHit(const EXTPCMeasLayer &ms,
                         Double_t       *x,
                         Double_t       *dx,
                         Int_t           side,
                         Double_t        v,
                   const void           *hitp,
                         Double_t        b,
                         Int_t           m)
        : TVTrackHit(ms, x, dx, b, m),
          fSide(side),
          fVdrift(v),
          fHitPtr(hitp)
{
}

EXTPCHit::~EXTPCHit()
{
}

//_________________________________________________________________________
//  --------------------------------
//  Implementation of public methods
//  --------------------------------
//
TKalMatrix EXTPCHit::XvToMv(const TVector3 &xv, Double_t t0) const
{
  const EXTPCMeasLayer &ms
    = dynamic_cast<const EXTPCMeasLayer &>(GetMeasLayer());
  TKalMatrix h    = ms.XvToMv(xv, GetSide());
  Double_t   r    = ms.GetR();
  Double_t   phih = h(0, 0) / r;
  Double_t   phim = (*this)(0, 0) / r;
  Double_t   dphi = phih - phim;

  static Double_t kPi    = TMath::Pi();
  static Double_t kTwoPi = 2 * kPi;

  while (dphi < -kPi) dphi += kTwoPi;
  while (dphi >  kPi) dphi -= kTwoPi;

  h(0, 0)  = r * (phim + dphi);
  h(1, 0) += fVdrift * t0;

  return h;
}

void EXTPCHit::DebugPrint(Option_t *) const
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
  cerr << " r    = " << setw(8)
       << static_cast<const EXTPCMeasLayer &>(GetMeasLayer()).GetR() << endl;
  cerr << "-------------------------------------------------------"  << endl;
}

//_________________________________________________________________________
//  ----------------
//  Compare two Hits
//  ----------------
//
Int_t EXTPCHit::Compare(const TObject *obj) const
{
  Double_t me  = GetSortingPolicy();
  Double_t you = dynamic_cast<const EXTPCHit *>(obj)->GetSortingPolicy();
  return me < you ? -1 : (me > you ? +1 : 0);
}

Double_t EXTPCHit::GetSortingPolicy() const
{
  // Calculate "r" (Caution!! vaild only for GEAR coordinate system)
  return GetMeasLayer().HitToXv(*this).Mag();
}

Bool_t EXTPCHit::IsSortable() const
{
  return kTRUE;
}
