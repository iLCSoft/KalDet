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
//*                                         hybrid/kern/EXVKalDetector.cxx
//*   2010/11/17  K.Fujii      Changed unit system to (mm, nsec, T) 
//*                            from (cm, nsec, kG)
//*
//* $Id: EXVKalDetector.cxx,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "EXVKalDetector.h"
#include "EXVMeasLayer.h"
#include "TVKalDetector.h"
#include "TTUBE.h"
#include "TNode.h"
#include "TRotMatrix.h"
#include "TVirtualPad.h"

#ifdef OLD_UNIT_SYSTEM
Double_t  EXVKalDetector::fgBfield  = 10.; // [kG]
#else
Double_t  EXVKalDetector::fgBfield  = 1.; // [T]
#endif
TNode    *EXVKalDetector::fgNodePtr = 0;

ClassImp(EXVKalDetector)

EXVKalDetector::EXVKalDetector(Int_t m)
              : TVKalDetector(m),
                fIsPowerOn(kTRUE)
{
}

EXVKalDetector::~EXVKalDetector()
{
}

TNode *EXVKalDetector::GetNodePtr()
{
  if (! fgNodePtr) {
    new TRotMatrix("rotm", "rotm", 10., 80., 10., 80., 10., 80.);
#ifdef OLD_UNIT_SYSTEM
    new TTUBE("Det", "Det", "void", 210., 210., 260.);
#else
    new TTUBE("Det", "Det", "void", 2100., 2100., 2600.);
#endif
    fgNodePtr = new TNode("World", "World", "Det", 0., 0., 0., "rotm");
  }
  return fgNodePtr;
}

void EXVKalDetector::Draw(Int_t color, const Char_t *opt)
{
  if (! gPad) return;
  TNode *nodep = GetNodePtr();
  nodep->cd();
  TIter next(this);
  TObject *objp;
  while ((objp = next())) {
    TAttDrawable *dp = dynamic_cast<TAttDrawable *>(objp);
    if (dp) dp->Draw(color, opt);
  }
  nodep->Draw("pad same");
}
