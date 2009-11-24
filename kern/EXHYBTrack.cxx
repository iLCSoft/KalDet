//*************************************************************************
//* ==================
//*  EXHYBTrack Class
//* ==================
//*
//* (Description)
//*   Hybrid track class for Kalman filter
//* (Requires)
//*     TKalTrack
//* (Provides)
//*     class EXHYBTrack
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/kern/EXHYBTrack.cxx
//*
//* $Id: EXHYBTrack.cxx,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "EXHYBTrack.h"
#include "TKalTrackSite.h"  // from KalTrackLib
#include "TVirtualPad.h"    // from ROOT
#include "TPolyMarker3D.h"  // from ROOT

//_________________________________________________________________________
//  ------------------------------
//  EXHYBTrack: Kalman Track class
//  ------------------------------
//
ClassImp(EXHYBTrack)

//_________________________________________________________________________
//  --------------
//  Utility Method
//  --------------
//_________________________________________________________________________
//  --------------------------------------
//  Draw: Drawing method for event display
//  --------------------------------------
//
void EXHYBTrack::Draw(Int_t color, const Char_t *opt)
{
  if (!gPad || !GetEntries()) return;
  gPad->cd();

  TPolyMarker3D *pm3dp = new TPolyMarker3D(this->GetEntries());
  pm3dp->SetBit(kCanDelete);
  pm3dp->SetMarkerColor(color);
  pm3dp->SetMarkerStyle(6);

  Int_t nhits = 0;
  TIter next(this);
  TKalTrackSite *sitep = 0;
  while ((sitep = static_cast<TKalTrackSite *>(next()))) {
    TVector3 pos = sitep->GetPivot();
    pm3dp->SetPoint(nhits, pos.X(), pos.Y(), pos.Z());
    nhits++;
  }
  pm3dp->Draw();
  gPad->Update();
}
