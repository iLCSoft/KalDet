#ifndef __ILDDISCMEASLAYER__
#define __ILDDISCMEASLAYER__
//*************************************************************************
//* ===================
//*  ILDDiscMeasLayer Class
//* ===================
//*
//* (Description)
//*   Disc measurement layer class used with ILDPLanarTrackHit.
//* (Requires)
//*   ILDVMeasLayer
//* (Provides)
//*     class ILDDiscMeasLayer
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

class ILDDiscMeasLayer : public ILDVMeasLayer, public TPlane {
public:
  // Ctors and Dtor
  
  ILDDiscMeasLayer(TMaterial &min,
                   TMaterial &mout,
                   const TVector3  &center,
                   const TVector3  &normal,
                   double   Bz,
                   double   SortingPolicy,
                   double   rMin,
                   double   rMax,
                   Bool_t     is_active,
                   Int_t      layerID = -1,
                   const Char_t    *name = "ILDDiscMeasL");
  
  
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
  
  double GetSortingPolicy() const { return _sortingPolicy; }
  
private:
  double _sortingPolicy;
  double _rMin;
  double _rMax;
  
};

#endif
