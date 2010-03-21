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
//* $Id: EXTPCKalDetector.h,v 1.2 2010-03-21 21:23:48 fujiik Exp $
//*************************************************************************
//
#include "EXVKalDetector.h"

class TNode;

class EXTPCKalDetector : public EXVKalDetector {
private:
  EXTPCKalDetector(Int_t m = 100);

public:
  ~EXTPCKalDetector();
  static EXTPCKalDetector * GetInstance();

  static Double_t GetVdrift() { return fgVdrift; }

  using EXVKalDetector::Draw;
  void  Draw(Int_t color, const Char_t *opt = "");

private:
         TNode            * fNodePtr;  // node pointer
  static Double_t           fgVdrift;   // drift velocity
  static EXTPCKalDetector * fgInstance; //! singleton pointer

  ClassDef(EXTPCKalDetector, 1)  // User defined detector class
};
#endif
