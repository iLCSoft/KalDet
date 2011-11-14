
#include <iostream>

#include "ILDDiscMeasLayer.h"
#include "ILDPlanarHit.h"

#include "TVTrack.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRotMatrix.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TString.h"

#include <EVENT/TrackerHitPlane.h>

#include "streamlog/streamlog.h"


TKalMatrix ILDDiscMeasLayer::XvToMv(const TVector3 &xv) const
{
  
  // Calculate measurement vector (hit coordinates) from global coordinates:
  
  TKalMatrix mv(kMdim,1);
  
  mv(0,0)  = xv.X() ;
  
  
  mv(1,0)  = xv.Y() ;
  return mv;
  
}


TVector3 ILDDiscMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  const ILDPlanarHit &mv = dynamic_cast<const ILDPlanarHit &>(vht);
  
  double x =   mv(0,0) ;
  double y =   mv(1,0) ;
  
  double z = GetXc().Z() ;
  
  return TVector3(x,y,z);
}

void ILDDiscMeasLayer::CalcDhDa(const TVTrackHit &vht,
                                const TVector3   &xxv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const
{
  // Calculate
  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
  // where
  //        h(a) = (phi, z)^t: expected meas vector
  //        a = (drho, phi0, kappa, dz, tanl, t0)
  //
  
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);
  
  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
  
  for (Int_t i=0; i<hdim; i++) {
    
    H(0,i) = dxphiada(0,i);
    H(1,i) = dxphiada(1,i) ;
    
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
    H(1,sdim-1) = 0.;
  }
  
}



Bool_t ILDDiscMeasLayer::IsOnSurface(const TVector3 &xx) const
{
  
  //  std::cout << "IsOnSurface " << std::endl;  
  
  bool onSurface = false ;
  
  TKalMatrix mv = XvToMv(xx);
  
  // check whether the hit lies in the same plane as the surface
  if( TMath::Abs((xx.X()-GetXc().X())*GetNormal().X() + (xx.Y()-GetXc().Y())*GetNormal().Y() + (xx.Z()-GetXc().Z())*GetNormal().Z()) < 1e-4){
    // check whether the hit lies within the boundary of the surface 
    
    double r2 = mv(0,0) * mv(0,0) + mv(1,0) * mv(1,0) ;
    
    if(  r2 <= _rMax*_rMax && r2 >= _rMin*_rMin )
        { 
          onSurface = true ;
        }    
  }
  
  return onSurface;
  
}


ILDVTrackHit* ILDDiscMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
  EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;
  
  if( plane_hit == NULL )  return NULL; // SJA:FIXME: should be replaced with an exception  
  
  const TVector3 hit( plane_hit->getPosition()[0], plane_hit->getPosition()[1], plane_hit->getPosition()[2]) ;
  
  // convert to layer coordinates       
  TKalMatrix h    = this->XvToMv(hit);
  
  double  x[2] ;
  double dx[2] ;
  
  x[0] = h(0, 0);
  x[1] = h(1, 0);
  
  dx[0] = plane_hit->getdU() ;
  dx[1] = plane_hit->getdV() ;
  
  bool hit_on_surface = IsOnSurface(hit);
  
  streamlog_out(DEBUG0) << "ILDDiscMeasLayer::ConvertLCIOTrkHit ILDPlanarHit created" 
  << " u = "  <<  x[0]
  << " v = "  <<  x[1]
  << " du = " << dx[0]
  << " dv = " << dx[1]
  << " x = " << plane_hit->getPosition()[0]
  << " y = " << plane_hit->getPosition()[1]
  << " z = " << plane_hit->getPosition()[2]
  << " onSurface = " << hit_on_surface
  << std::endl ;
  
  return hit_on_surface ? new ILDPlanarHit( *this , x, dx, this->GetBz()) : NULL; 
  
}
