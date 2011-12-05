
#include "ILDParallelStripPlanarMeasLayer.h"
#include "ILDPlanarStripHit.h"
#include "TVTrack.h"
#include "TVTrackHit.h"

#include <EVENT/TrackerHitPlane.h>

#include "gearimpl/Vector3D.h"

#include "streamlog/streamlog.h"


TKalMatrix ILDParallelStripPlanarMeasLayer::XvToMv(const TVector3 &xv) const
{
  // Calculate hit coordinate information:
  //	mv(0,0) = u (transverse)

  
  TKalMatrix mv(1,1);
  mv(0,0) = (xv.Y() - GetXc().Y())*GetNormal().X()/GetNormal().Perp() - (xv.X() - GetXc().X())*GetNormal().Y()/GetNormal().Perp() ;

  return mv;

}

TKalMatrix ILDParallelStripPlanarMeasLayer::XvToMv(const TVTrackHit &,
                                 const TVector3   &xv) const
{
  return XvToMv(xv);
}

TVector3 ILDParallelStripPlanarMeasLayer::HitToXv(const TVTrackHit &vht) const
{

  double x = -vht(0,0) * this->GetNormal().Y() / this->GetNormal().Perp() + this->GetXc().X();
  double y =  vht(0,0) * this->GetNormal().X() / this->GetNormal().Perp() + this->GetXc().Y();

  double z   = this->GetXc().Z();
  
  return TVector3(x,y,z);

}

void ILDParallelStripPlanarMeasLayer::CalcDhDa(const TVTrackHit &vht,
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
    H(0,i) = (this->GetNormal().X() / this->GetNormal().Perp()) * dxphiada(1,i)
    -(this->GetNormal().Y() / this->GetNormal().Perp()) * dxphiada(0,i);   
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
  }
  
}

ILDVTrackHit* ILDParallelStripPlanarMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
  EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;

  if( plane_hit == NULL )  return NULL; // SJA:FIXME: should be replaced with an exception  
  
  gear::Vector3D U(1.0,plane_hit->getU()[1],plane_hit->getU()[0],gear::Vector3D::spherical);
  gear::Vector3D Z(0.0,0.0,1.0);
  
  const float eps = 1.0e-07;

  // U must be in the transverse plane
  if( fabs(U.dot(Z)) > eps ) {
    streamlog_out(ERROR) << "ILDParallelStripPlanarMeasLayer: TrackerHitPlane measurment vectors U is not in the global X-Y plane. exit(1) called from file " << __FILE__ << " and line " << __LINE__ << std::endl;
    exit(1);
  }

  if( plane_hit == NULL )  return NULL; // SJA:FIXME: should be replaced with an exception  
  
  // remember here the "position" of the hit in fact defines the origin of the plane it defines so u and v are per definition 0. 
  // this is still the case for a 1-dimentional measurement, and is then used to calculate the u coordinate according to the origin of the actual measurement plane.
  const TVector3 hit( plane_hit->getPosition()[0], plane_hit->getPosition()[1], plane_hit->getPosition()[2]) ;
  
  // convert to layer coordinates       
  TKalMatrix h    = this->XvToMv(hit);
  
  Double_t  x[1] ;
  Double_t dx[1] ;
  
  x[0] = h(0, 0);
  
  dx[0] = plane_hit->getdU() ;
  
  bool hit_on_surface = IsOnSurface(hit);
  
  streamlog_out(DEBUG0) << "ILDPlanarMeasLayer::ConvertLCIOTrkHit ILDPlanarHit created" 
  << " Layer R = " << this->GetXc().Mag() 
  << " Layer phi = " << this->GetXc().Phi() 
  << " u = "  <<  x[0]
  << " du = " << dx[0]
  << " x = " << plane_hit->getPosition()[0]
  << " y = " << plane_hit->getPosition()[1]
  << " z = " << plane_hit->getPosition()[2]
  << " r = " << sqrt( plane_hit->getPosition()[0]*plane_hit->getPosition()[0] + plane_hit->getPosition()[1]*plane_hit->getPosition()[1])
  << " onSurface = " << hit_on_surface
  << std::endl ;
  
  
  return hit_on_surface ? new ILDPlanarStripHit( *this , x, dx, this->GetBz()) : NULL; 
  
  
}


