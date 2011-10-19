#ifndef ILDVTrackHIT_H
#define ILDVTrackHIT_H
//*************************************************************************
//* ================
//*  ILDVTrackHit Class
//* ================
//*
//* (Description)
//*    Virtual hit class used by ILD[X]Hit Classes.
//*   provides coordinate vector as defined by the MeasLayer
//* (Requires)
//*     TVTrackHit
//* (Provides)
//*     class ILDVTrackHit
//*
//*************************************************************************
//

#include "kaltest/TVTrackHit.h"

#include "ILDVMeasLayer.h"


class ILDVTrackHit : public TVTrackHit {
  
public:
  
  ILDVTrackHit(const TVMeasLayer &ms, Double_t *x, Double_t *dx, 
               Double_t bfield , Int_t dim) ; 
  
  virtual ~ILDVTrackHit();
  
  
private:
  
  
};
#endif
