#ifndef EXTPCMEASLAYER_H
#define EXTPCMEASLAYER_H
//*************************************************************************
//* ======================
//*  EXTPCMeasLayer Class
//* ======================
//*
//* (Description)
//*   User defined measurement layer class
//* (Requires)
//*     EXVMeasLayer
//* (Provides)
//*     class EXTPCMeasLayer
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/tpc/EXTPCMeasLayer.h
//*
//* $Id: EXTPCMeasLayer.h,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TCylinder.h"
#include "EXVMeasLayer.h"
#include "KalTrackDim.h"

class TVTrackHit;

class EXTPCMeasLayer : public EXVMeasLayer, public TCylinder {

public:
  EXTPCMeasLayer(TMaterial &min,
                 TMaterial &mout,
                 Double_t   r0,
                 Double_t   lhalf,
                 Double_t   sigmax0,
                 Double_t   sigmax1,
                 Double_t   sigmaz0,
                 Double_t   sigmaz1,
                 Bool_t     type = EXVMeasLayer::kActive);
  EXTPCMeasLayer(TMaterial &min,
                 TMaterial &mout,
                 Double_t   r0,
                 TVector3   xc,
                 Double_t   phimin,
                 Double_t   phimax,
                 Double_t   lhalf,
                 Double_t   sigmax0,
                 Double_t   sigmax1,
                 Double_t   sigmaz0,
                 Double_t   sigmaz1,
                 Bool_t     type = EXVMeasLayer::kDummy,
                 Int_t      module = -1,
                 Int_t      layer = -1);
  virtual ~EXTPCMeasLayer();

  inline Int_t GetModuleID() const { return fModule; }
  inline Int_t GetLayerID () const { return fLayer;  }

  // Parent's pure virtuals that must be implemented

  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv)   const;
  virtual TKalMatrix XvToMv    (const TVector3   &xv,
                                      Int_t       side) const;
  virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                      TKalMatrix &H)    const;
  virtual void       ProcessHit(const TVector3   &xx,
                                      TObjArray  &hits);

  virtual Double_t   GetSortingPolicy() const;

  Double_t GetSigmaX(Double_t z) const;
  Double_t GetSigmaZ(Double_t z) const;

private:
  Double_t fPhiMin;   // minimum phi
  Double_t fPhiMax;   // maximum phi
  Double_t fSigmaX0;  // xy resolution
  Double_t fSigmaX1;  // xy resolution
  Double_t fSigmaZ0;  // z  resolution
  Double_t fSigmaZ1;  // z  resolution
  Int_t    fModule;   // module number
  Int_t    fLayer;    // layer  number

  ClassDef(EXTPCMeasLayer, 1)  // User defined measurement layer class
};
#endif
