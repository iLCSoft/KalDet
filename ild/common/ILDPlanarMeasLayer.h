#ifndef __ILDPLANARMEASLAYER__
#define __ILDPLANARMEASLAYER__
//*************************************************************************
//* ===================
//*  ILDPlanarMeasLayer Class
//* ===================
//*
//* (Description)
//*   Planar measurement layer class used with ILDPLanarTrackHit.
//* (Requires)
//*   ILDVMeasLayer
//* (Provides)
//*     class ILDPlanarMeasLayer
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*   2011/06/17  D.Kamai           Modified to handle ladder structure.
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

class ILDPlanarMeasLayer : public ILDVMeasLayer, public TPlane {
public:
  // Ctors and Dtor
  
  ILDPlanarMeasLayer(TMaterial &min,
                     TMaterial &mout,
                     const TVector3  &center,
                     const TVector3  &normal,
                     Double_t   Bz,
                     Double_t   SortingPolicy,
                     Double_t   xiwidth,
                     Double_t   zetawidth,
                     Double_t   xioffset,
                     Bool_t     is_active,
                     Int_t      layerID = -1,
                     const Char_t    *name = "ILDPlanarMeasL");
  
  virtual ~ILDPlanarMeasLayer();
  
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
  
  virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  Double_t GetSortingPolicy() const { return fSortingPolicy; }
  Double_t GetXiwidth() const { return fXiwidth; }
  Double_t GetZetawidth() const { return fZetawidth; }
  Double_t GetXioffset() const { return fXioffset; }
  
private:
  Double_t fSortingPolicy;
  Double_t fXiwidth;
  Double_t fZetawidth;
  Double_t fXioffset;
  
};

#endif
