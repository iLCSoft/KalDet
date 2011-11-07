#ifndef __ILDPLANARMEASLAYER__
#define __ILDPLANARMEASLAYER__
//*************************************************************************
//* ===================
//*  ILDRotatedTrapMeaslayer Class
//* ===================
//*
//* (Description)
//*   Rotated Trapezoid Planar measurement layer class used with ILDPLanarTrackHit.
//* (Requires)
//*   ILDVMeasLayer
//* (Provides)
//*     class ILDRotatedTrapMeaslayer
//*
//*************************************************************************
//
#include "TVector3.h"
#include "TKalMatrix.h"
#include "TPlane.h"
#include "ILDVMeasLayer.h"
#include "KalTrackDim.h"
#include "TMath.h"
#include <sstream>
class TVTrackHit;

class ILDRotatedTrapMeaslayer : public ILDVMeasLayer, public TPlane {
public:
  // Ctors and Dtor
  
  ILDRotatedTrapMeaslayer(TMaterial &min,
                          TMaterial &mout,
                          const TVector3  &center,
                          const TVector3  &normal,
                          Double_t   Bz,
                          Double_t   SortingPolicy,
                          Double_t   height,
                          Double_t   innerBaseLength,
                          Double_t   outerBaseLength,
                          Double_t   alpha,
                          Int_t      half_petal,
                          Bool_t     is_active,
                          Int_t      layerID = -1,
                          const Char_t    *name = "ILDRotatedTrapMeasL");
  
  
  // Parrent's pure virtuals that must be implemented
  
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv) const;
  
  virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
  
  virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
  
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const;
  
  virtual ILDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const ;
  
  inline virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  Double_t GetSortingPolicy() const { return _sortingPolicy; }
  
private:
  Double_t _sortingPolicy ;
  Double_t _signZ ;
  Double_t _innerR ;
  Double_t _outerR ;
  Double_t _innerBaseLength ;
  Double_t _outerBaseLength ;
  Double_t _cosPhi ;  //** cos of the azimuthal angle of the petal 
  Double_t _sinPhi ;  //** sin of the azimuthal angle of the petal 
  Double_t _cosAlpha ; //** cos of the tilt angle of the petal 
  Double_t _sinAlpha ; //** sin of the tilt angle of the petal 
  Double_t _tanBeta ; //** tan of the openning angle of the petal
  
  // meaning of _halfPetal:
  //                  0 complete trapezoid
  //                 +1 positive half only, i.e. the side of the petal in which the transverse coordinate, mv(0,0), is positive 
  //                 -1 negative half only, i.e. the side of the petal in which the transverse coordinate, mv(0,0), is negative
  Int_t _halfPetal ;
  
  
  
};

#endif
