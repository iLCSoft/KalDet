/** User defined KalTest measurement layer class 
 *
 * @author S.Aplin DESY
 */


#include "TKalTrack.h" 

#include "ILDVTrackHit.h"
#include "ILDCylinderMeasLayer.h"
#include "ILDCylinderHit.h"


#include <lcio.h>
#include <EVENT/TrackerHit.h>

#include "streamlog/streamlog.h"

#include "TMath.h"
#include <cmath>


/** Global to Local coordinates */

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


/** Local to Global coordinates */

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


/** Calculate Projector Matrix */

void ILDCylinderMeasLayer::CalcDhDa(const TVTrackHit &vht, // tracker hit not used here
                                    const TVector3   &xxv,
                                    const TKalMatrix &dxphiada,
                                    TKalMatrix &H) const
{
  
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


/** Convert LCIO Tracker Hit to an ILDCylinderHit  */

ILDVTrackHit* ILDCylinderMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
  
  const TVector3 hit( trkhit->getPosition()[0], trkhit->getPosition()[1], trkhit->getPosition()[2]) ;
  
  // convert to layer coordinates       
  TKalMatrix h    = this->XvToMv(hit);
  
  Double_t  x[2] ;
  Double_t dx[2] ;
  
  x[0] = h(0, 0);
  x[1] = h(1, 0);
  
  // convert errors
  dx[0] = sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]) ;
  dx[1] = sqrt(trkhit->getCovMatrix()[5]) ; 
  
  bool hit_on_surface = IsOnSurface(hit);
  
  streamlog_out(DEBUG0) << "ILDCylinderMeasLayer::ConvertLCIOTrkHit ILDCylinderHit created" 
  << " R = " << hit.Perp()
  << " Layer R = " << this->GetR() 
  << " RPhi = "  <<  x[0]
  << " Z = "     <<  x[1]
  << " dRPhi = " << dx[0]
  << " dZ = "    << dx[1]
  << " x = " << trkhit->getPosition()[0]
  << " y = " << trkhit->getPosition()[1]
  << " z = " << trkhit->getPosition()[2]
  << " onSurface = " << hit_on_surface
  << std::endl ;  
  
  return hit_on_surface ? new ILDCylinderHit( *this , x, dx, this->GetBz()) : NULL; 
  
}

