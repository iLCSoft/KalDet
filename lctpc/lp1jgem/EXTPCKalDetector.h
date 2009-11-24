#ifndef EXTPCDETECTOR_H
#define EXTPCDETECTOR_H
//*************************************************************************
//* ========================
//*  EXTPCKalDetector Class
//* ========================
//*
//* (Description)
//*   User defined detector class
//* (Requires)
//*     EXVKalDetector
//* (Provides)
//*     class EXTPCKalDetector
//* (Update Recored)
//*   2009/11/23  K.Ikematsu   Derived from KalTest/examples/kaltest/
//*                                         hybrid/tpc/EXTPCKalDetector.h
//*
//* $Id: EXTPCKalDetector.h,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "EXVKalDetector.h"

class TNode;

class EXTPCKalDetector : public EXVKalDetector {

public:
  EXTPCKalDetector(Int_t m = 100);
  ~EXTPCKalDetector();

  static Double_t GetVdrift() { return fgVdrift; }

  using EXVKalDetector::Draw;
  void  Draw(Int_t color, const Char_t *opt = "");

private:
  TNode *fNodePtr;           // node pointer
  static Double_t fgVdrift;  // drift velocity

  ClassDef(EXTPCKalDetector, 1)  // User defined detector class
};
#endif
