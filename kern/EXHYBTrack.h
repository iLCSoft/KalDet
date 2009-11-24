#ifndef EXHYBTRACK_H
#define EXHYBTRACK_H
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
//*                                         hybrid/kern/EXHYBTrack.h
//*
//* $Id: EXHYBTrack.h,v 1.1.1.1 2009-11-24 00:13:59 ikematsu Exp $
//*************************************************************************
//
#include "TKalTrack.h"     // from KalTrackLib
#include "TAttDrawable.h"  // from Utils

//_________________________________________________________________________
//  ------------------------------
//  EXHYBTrack: Kalman Track class
//  ------------------------------

class EXHYBTrack : public TKalTrack, public TAttDrawable {

public:
  EXHYBTrack(Int_t n = 1) : TKalTrack(n) {}
  ~EXHYBTrack() {}

  using TAttDrawable::Draw;
  virtual void Draw(Int_t color, const Char_t *opt);

  ClassDef(EXHYBTrack, 1)  // Hybrid track class for Kalman filter
};
#endif
