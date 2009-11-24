#ifndef EXVKALDETECTOR_H
#define EXVKALDETECTOR_H
//*************************************************************************
//* ======================
//*  EXVKalDetector Class
//* ======================
//*
//* (Description)
//*   Abstract detector class for Kalman filter
//* (Requires)
//*     TVKalDetector
//* (Provides)
//*     class EXVKalDetector
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/kern/EXVKalDetector.h
//*
//* $Id: EXVKalDetector.h,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "TVector3.h"
#include "TVKalDetector.h"
#include "TAttDrawable.h"

class TVMeasLayer;
class TNode;

class EXVKalDetector : public TVKalDetector, public TAttDrawable {

public:
  EXVKalDetector(Int_t m = 100);
  virtual ~EXVKalDetector();

  inline virtual Bool_t IsPowerOn() const { return fIsPowerOn;   }
  inline virtual void   PowerOn  ()       { fIsPowerOn = kTRUE;  }
  inline virtual void   PowerOff ()       { fIsPowerOn = kFALSE; }

  static Double_t GetBfield (const TVector3 &xx = TVector3(0.))
                            { return fgBfield; }

  using TAttDrawable::Draw;
  virtual void Draw(Int_t color, const Char_t *opt = "");

  static void   SetNodePtr(TNode *nodep) { fgNodePtr = nodep; }
  static TNode *GetNodePtr();

private:
  Bool_t fIsPowerOn;           // power status
  static Double_t  fgBfield;   // magnetic field [kG]
  static TNode    *fgNodePtr;  // pointer to TNode

  ClassDef(EXVKalDetector, 1)  // Abstract detector class
};
#endif
