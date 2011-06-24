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
//* (Update Recored)
//*
//*************************************************************************
//

#include "TKalTrack.h" 

#include "ILDCylinderMeasLayer.h"
#include "ILDCylinderHit.h"

#include "TMath.h"
#include <cmath>

ILDCylinderMeasLayer::ILDCylinderMeasLayer(TMaterial &min,
				 TMaterial &mout,
				 Double_t   r0,
				 Double_t   lhalf,
				 Double_t   Bz,
				 Bool_t is_active,
				 Int_t      layerID,
				 const Char_t    *name) 
  : ILDVMeasLayer(min, mout, Bz, is_active, layerID, name),
  TCylinder(r0, lhalf)
{
}


ILDCylinderMeasLayer::~ILDCylinderMeasLayer()
{
}

TKalMatrix ILDCylinderMeasLayer::XvToMv(const TVector3 &xv) const
{

  // Calculate hit coordinate information:
  //   mv(0, 0) = r * phi
  //     (1, 0) = drift distance

  // account for cylinder not centered at x=0.0, y=0.0
  TVector3 xxv = xv - GetXc();

  Double_t phi = TMath::ATan2(xxv.Y(), xxv.X());

  // bring phi back into +/- Pi range
  static Double_t kPi    = TMath::Pi();
  static Double_t kTwoPi = 2 * kPi;
  while (phi < -kPi) phi += kTwoPi;
  while (phi >  kPi) phi -= kTwoPi;

  TKalMatrix mv(kMdim, 1);

  mv(0, 0) = GetR() * phi;

  mv(1, 0) = xxv.Z();

  return mv;
}

TKalMatrix ILDCylinderMeasLayer::XvToMv(const TVTrackHit &vht, // Hit is not used here
                                  const TVector3   &xv) const
{
  return XvToMv(xv);
}

TVector3 ILDCylinderMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  const ILDCylinderHit &ht = dynamic_cast<const ILDCylinderHit &>(vht) ;

  Double_t phi = ht(0, 0) / GetR() ;
  Double_t z   = ht(1, 0);

  // account for cylinder not centered at x=0.0, y=0.0
  Double_t x   = GetR() * TMath::Cos(phi) + GetXc().X();
  Double_t y   = GetR() * TMath::Sin(phi) + GetXc().Y();

  return TVector3(x, y, z);
}

void ILDCylinderMeasLayer::CalcDhDa(const TVTrackHit &vht, // tracker hit not used here
                              const TVector3   &xxv,
                              const TKalMatrix &dxphiada,
                                    TKalMatrix &H) const
{

  //  const ILDCylinderHit &ht = dynamic_cast<const ILDCylinderHit &>(vht);

  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  //  where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //

  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5, sdim - 1);

  // account for cylinder not centered at x=0.0, y=0.0
  TVector3 xxvc = xxv - GetXc();

  Double_t xv   = xxvc.X();
  Double_t yv   = xxvc.Y();
  Double_t xxyy = xv * xv + yv * yv;

  // Set H = (@h/@a) = (@d/@a, @z/@a)^t

  for (Int_t i = 0; i < hdim; i++) {
    H(0, i)  = - (yv / xxyy) * dxphiada(0, i)
               + (xv / xxyy) * dxphiada(1, i);
    H(0, i) *= GetR();

    H(1, i)  = dxphiada(2, i);
  }

  if (sdim == 6) {
    H(0, sdim - 1) = 0.;
  }

}


