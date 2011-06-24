#ifndef ILDCYLINDERMEASLAYER_H
#define ILDCYLINDERMEASLAYER_H
//*************************************************************************
//* ======================
//*  ILDCylinderMeasLayer Class
//* ======================
//*
//* (Description)
//*   User defined measurement layer class
//* (Requires)
//*     ILDVMeasLayer
//* (Provides)
//*     class ILDCylinderMeasLayer
//*
//*************************************************************************
//

#include "ILDVMeasLayer.h"

class ILDCylinderMeasLayer : public ILDVMeasLayer, public TCylinder {

public:
  ILDCylinderMeasLayer(TMaterial &min,
		  TMaterial &mout,
		  Double_t   r0,
		  Double_t   lhalf,
		  Double_t   Bz,
		  Bool_t     is_active,
		  Int_t      layerID = -1,
		  const Char_t    *name = "ILDCylinderMeasL") ;
  
  virtual ~ILDCylinderMeasLayer();
  
  virtual TKalMatrix XvToMv    (const TVector3   &xv)   const;

  // Parent's pure virtuals that must be implemented

  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv)   const;

  virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;

  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                      TKalMatrix &H)    const;

private:
  
};
#endif
