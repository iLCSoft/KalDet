
#include "ILDParallelStripPlanarMeasLayer.h"
#include "ILDPlanarStripHit.h"
#include "TVTrack.h"
#include "TVTrackHit.h"

#include <EVENT/TrackerHitPlane.h>

#include "gearimpl/Vector3D.h"

#include "streamlog/streamlog.h"


TKalMatrix ILDParallelStripPlanarMeasLayer::XvToMv(const TVector3 &xv) const
{
  
  double cos_phi = GetNormal().X()/GetNormal().Perp();
  double sin_phi = GetNormal().Y()/GetNormal().Perp();
  
  // theta is the strip angle rotation in the plane of the wafer
  double cos_theta = cos(_stripAngle);
  double sin_theta = sin(_stripAngle);
  
  double delta_x = xv.X() - GetXc().X();
  double delta_y = xv.Y() - GetXc().Y();

  double delta_t = (delta_x * sin_phi - delta_y * cos_phi) ; 
  double delta_z = xv.Z() - GetXc().Z();


  TKalMatrix mv(ILDPlanarStripHit_DIM,1);
  
  mv(0,0) = ( (delta_t + fUOrigin) * cos_theta ) + delta_z * sin_theta ;
  
  if (ILDPlanarStripHit_DIM == 2) {
    mv(1,0) = delta_z * cos_theta - delta_t * sin_theta;
  }
  
  return mv;

}

TKalMatrix ILDParallelStripPlanarMeasLayer::XvToMv(const TVTrackHit &,
                                 const TVector3   &xv) const
{
  return XvToMv(xv);
}

TVector3 ILDParallelStripPlanarMeasLayer::HitToXv(const TVTrackHit &vht) const
{

  
  double cos_phi = GetNormal().X()/GetNormal().Perp();
  double sin_phi = GetNormal().Y()/GetNormal().Perp();
  
  // theta is the strip angle rotation in the plane of the wafer
  double cos_theta = cos(_stripAngle);
  double sin_theta = sin(_stripAngle);
  
  double t =  (vht(0,0) - fUOrigin ) * cos_theta ;
    
  double x =  t * sin_phi + this->GetXc().X();
  double y = -t * cos_phi + this->GetXc().Y();
  
  double dz = 0.0;
  
  if (ILDPlanarStripHit_DIM == 2) {
    dz = vht(1,0) * cos_theta ;
  }

  
  double z = vht(0,0) * sin_theta + this->GetXc().Z() + dz;

  this->IsOnSurface(TVector3(x,y,z));
  
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
  
  double cos_phi = GetNormal().X()/GetNormal().Perp();
  double sin_phi = GetNormal().Y()/GetNormal().Perp();

  // theta is the strip angle rotation in the plane of the wafer
  double cos_theta = cos(_stripAngle);
  double sin_theta = sin(_stripAngle);
  
  Int_t sdim = H.GetNcols();
  Int_t hdim = TMath::Max(5,sdim-1);
  
  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
  
  for (Int_t i=0; i<hdim; i++) {
    
    H(0,i) =  cos_theta * sin_phi * dxphiada(0,i) - cos_theta * cos_phi * dxphiada(1,i) + sin_theta * dxphiada(2,i) ;   

    if (ILDPlanarStripHit_DIM == 2) {
      H(1,i) = -sin_theta * sin_phi * dxphiada(0,i) + sin_theta * cos_phi * dxphiada(1,i) + cos_theta * dxphiada(2,i);
    }
    
  }
  if (sdim == 6) {
    H(0,sdim-1) = 0.;
    if (ILDPlanarStripHit_DIM == 2) {
      H(1,sdim-1) = 0.;
    }
  }  
}

ILDVTrackHit* ILDParallelStripPlanarMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
  EVENT::TrackerHitPlane* plane_hit = dynamic_cast<EVENT::TrackerHitPlane*>( trkhit ) ;

  if( plane_hit == NULL )  return NULL; // SJA:FIXME: should be replaced with an exception  
  
  gear::Vector3D U(1.0,plane_hit->getU()[1],plane_hit->getU()[0],gear::Vector3D::spherical);
  gear::Vector3D Z(0.0,0.0,1.0);
  
  const float eps = 1.0e-07;


  if( plane_hit == NULL )  return NULL; // SJA:FIXME: should be replaced with an exception  
  
  // remember here the "position" of the hit in fact defines the origin of the plane it defines so u and v are per definition 0. 
  // this is still the case for a 1-dimentional measurement, and is then used to calculate the u coordinate according to the origin of the actual measurement plane.
  const TVector3 hit( plane_hit->getPosition()[0], plane_hit->getPosition()[1], plane_hit->getPosition()[2]) ;
  
  // convert to layer coordinates       
  TKalMatrix h(ILDPlanarStripHit_DIM,1);    
  h = this->XvToMv(hit);
  
  Double_t  x[ILDPlanarStripHit_DIM] ;
  Double_t dx[ILDPlanarStripHit_DIM] ;
  
  x[0] = h(0, 0);
  if(ILDPlanarStripHit_DIM == 2) x[1] = h(1, 0);

  dx[0] = plane_hit->getdU() ;
  if(ILDPlanarStripHit_DIM == 2) dx[1] = plane_hit->getdV() ;
  
  bool hit_on_surface = IsOnSurface(hit);
  
  streamlog_out(DEBUG0) << "ILDParallelStripPlanarMeasLayer::ConvertLCIOTrkHit ILDPlanarStripHit created" 
  << " for CellID " << trkhit->getCellID0()
  << " Layer R = " << this->GetXc().Perp() 
  << " Layer phi = " << this->GetXc().Phi() 
  << " Layer z0 = " << this->GetXc().Z() 
  << " u = "  <<  x[0]
  //  << " v = "  <<  x[1]
  << " du = " << dx[0]
  //  << " dv = " << dx[1]
  << " x = " << plane_hit->getPosition()[0]
  << " y = " << plane_hit->getPosition()[1]
  << " z = " << plane_hit->getPosition()[2]
  << " r = " << sqrt( plane_hit->getPosition()[0]*plane_hit->getPosition()[0] + plane_hit->getPosition()[1]*plane_hit->getPosition()[1])
  << " onSurface = " << hit_on_surface
  << std::endl ;
  
  
  return hit_on_surface ? new ILDPlanarStripHit( *this , x, dx, this->GetBz(),trkhit) : NULL; 
  
  
}


