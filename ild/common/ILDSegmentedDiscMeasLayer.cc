
#include "ILDSegmentedDiscMeasLayer.h"
#include "ILDPlanarHit.h"

#include "TVTrack.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRotMatrix.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TString.h"

#include <EVENT/TrackerHitPlane.h>

#include <math.h>
#include <assert.h>
#include <algorithm>

#include "streamlog/streamlog.h"


ILDSegmentedDiscMeasLayer::ILDSegmentedDiscMeasLayer(TMaterial &min,
                                                     TMaterial &mout,
                                                     double   Bz,
                                                     double   sortingPolicy,
                                                     int      nsegments,
                                                     double   zpos,
                                                     double   phi0, // defined by the axis of symmerty of the first petal
                                                     double   trap_rmin,
                                                     double   trap_height,
                                                     double   trap_inner_base_length,
                                                     double   trap_outer_base_length,
                                                     bool     is_active,
                                                     std::vector<int>      CellIDs,
                                                     const Char_t    *name) : 
ILDVMeasLayer(min, mout, Bz, CellIDs, is_active, name),
TPlane(TVector3(0.,0.,zpos), TVector3(0.,0.,zpos)),
_sortingPolicy(sortingPolicy),_nsegments(nsegments),_trap_rmin(trap_rmin),_trap_height(trap_height),_trap_inner_base_length(trap_inner_base_length),_trap_outer_base_length(trap_outer_base_length)
{
  
  double h = _trap_rmin + _trap_height;
  double w = 0.5 * _trap_outer_base_length;
  _rmax = sqrt(h*h + w*w); 
  
  _trap_tan_beta = 0.5*(_trap_outer_base_length - _trap_inner_base_length) / _trap_height;
  
  _segment_dphi = 2.0*M_PI / _nsegments; 
  
  phi0 = angular_range_2PI(phi0);
  
  _start_phi = phi0 - 0.5*_segment_dphi;
  
  _start_phi = angular_range_2PI(_start_phi);
  
  // now check for constistency
  double phi_max = std::max( ((0.5 * _trap_outer_base_length) / h ) , ((0.5 * _trap_inner_base_length) / _trap_rmin));
  
  if( phi_max > _segment_dphi ) {
    streamlog_out(ERROR) << "ILDSegmentedDiscMeasLayer::ILDSegmentedDiscMeasLayer trapezoids overlaps: exit(1) called from " << __FILE__ << "   line " << __LINE__ << std::endl; 
    exit(1);
  }
  
}






ILDSegmentedDiscMeasLayer::ILDSegmentedDiscMeasLayer(TMaterial &min,
                                                     TMaterial &mout,
                                                     double   Bz,
                                                     double   sortingPolicy,
                                                     int      nsegments,
                                                     double   zpos,
                                                     double   phi0, // defined by the axis of symmerty of the first petal
                                                     double   trap_rmin,
                                                     double   trap_height,
                                                     double   trap_inner_base_length,
                                                     double   trap_outer_base_length,
                                                     bool     is_active,
                                                     const Char_t    *name) : 
ILDVMeasLayer(min, mout, Bz, is_active, -1, name),
TPlane(TVector3(0.,0.,zpos), TVector3(0.,0.,zpos)),
_sortingPolicy(sortingPolicy),_nsegments(nsegments),_trap_rmin(trap_rmin),_trap_height(trap_height),_trap_inner_base_length(trap_inner_base_length),_trap_outer_base_length(trap_outer_base_length)
{
  
  double h = _trap_rmin + _trap_height;
  double w = 0.5 * _trap_outer_base_length;
  _rmax = sqrt(h*h + w*w); 

  _trap_tan_beta = 0.5*(_trap_outer_base_length - _trap_inner_base_length) / _trap_height;
  
  _segment_dphi = 2.0*M_PI / _nsegments; 
  
  phi0 = angular_range_2PI(phi0);
  
  _start_phi = phi0 - 0.5*_segment_dphi;
 
  _start_phi = angular_range_2PI(_start_phi);
  
  // now check for constistency
  double phi_max = std::max( ((0.5 * _trap_outer_base_length) / h ) , ((0.5 * _trap_inner_base_length) / _trap_rmin));

  if( phi_max > _segment_dphi ) {
    streamlog_out(ERROR) << "ILDSegmentedDiscMeasLayer::ILDSegmentedDiscMeasLayer trapezoids overlaps: exit(1) called from " << __FILE__ << "   line " << __LINE__ << std::endl; 
    exit(1);
  }
  
}


TKalMatrix ILDSegmentedDiscMeasLayer::XvToMv(const TVector3 &xv) const
{
  
  // Calculate measurement vector (hit coordinates) from global coordinates:
  
  TKalMatrix mv(kMdim,1);
  
  mv(0,0)  = xv.X() ;
  
  
  mv(1,0)  = xv.Y() ;
  return mv;
  
}


TVector3 ILDSegmentedDiscMeasLayer::HitToXv(const TVTrackHit &vht) const
{
  const ILDPlanarHit &mv = dynamic_cast<const ILDPlanarHit &>(vht);
  
  //SJA:FIXME: in order to use real local coordinates we would have to get the CELLID from the ILDPlanarHit, this would tell us in which segment the hit was in 
  
  double x =   mv(0,0) ;
  double y =   mv(1,0) ;
  
  double z = this->GetXc().Z() ;
  
  return TVector3(x,y,z);
}

void ILDSegmentedDiscMeasLayer::CalcDhDa(const TVTrackHit &vht,
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


Int_t ILDSegmentedDiscMeasLayer::CalcXingPointWith(const TVTrack  &hel,
                                                    TVector3 &xx,
                                                    Double_t &phi,
                                                    Int_t     mode,
                                                    Double_t  eps) const{

  //  streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::CalcXingPointWith" << std::endl;
  
  // check that direction has one of the correct values
  if( !( mode == 0 || mode == 1 || mode == -1) ) return -1 ;
  
//  
//  TVector3 xx_n;
//  int cuts = TVSurface::CalcXingPointWith(hel, xx_n, phi, mode, eps);
//  streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::CalcXingPointWith from Newton: cuts = " << cuts << " x = " << xx_n.x() << " y = "<< xx_n.y() << " z = " << xx_n.z() << " r = " << xx_n.Perp() << " phi = " << xx_n.Phi() << " dphi = " <<  phi << std::endl;   
  
  // get helix parameters
  Double_t dr     = hel.GetDrho();
  Double_t phi0   = hel.GetPhi0(); //
  Double_t kappa  = hel.GetKappa();
  Double_t rho    = hel.GetRho();
  Double_t omega  = 1.0 / rho;
  Double_t r      = TMath::Abs(rho);
  Double_t z0     = hel.GetDz();
  Double_t tanl   = hel.GetTanLambda();
  
  TVector3 ref_point = hel.GetPivot();
  
  //
  // Check if charge is nonzero.
  //
  
  Int_t    chg = (Int_t)TMath::Sign(1.1,kappa);
  if (!chg) {
    streamlog_out(ERROR) << ">>>> Error >>>> ILDParallelPlanarMeasLayer::CalcXingPointWith" << std::endl
    << "      Kappa = 0 is invalid for a helix "          << std::endl;
    return -1;
  }
  
  const double sin_phi0 = sin(phi0); 
  const double cos_phi0 = cos(phi0); 
  
  const double x_pca = ref_point.x() + dr * cos_phi0 ; 
  const double y_pca = ref_point.y() + dr * sin_phi0 ; 
  const double z_pca = ref_point.z() + z0 ;
  
  const double z = this->GetXc().Z() ;
  // get path length to crossing point 

  const double s = ( z - z_pca ) / tanl ;
  
  const double delta_phi_half = (omega*s)/2.0 ;
  
  double x;
  double y;
  
  if( fabs(s) > FLT_MIN ){ // protect against starting on the plane
    x = x_pca - s * ( sin(delta_phi_half) / delta_phi_half ) *  sin( phi0 - delta_phi_half ) ;
    y = y_pca + s * ( sin(delta_phi_half) / delta_phi_half ) *  cos( phi0 - delta_phi_half ) ;
  }
  else{
    x = x_pca;
    y = y_pca;
  }
  
  
  // check if intersection with plane is within boundaries
  
  xx.SetXYZ(x, y, z);
  

  phi = s / r ;
 
//  streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::CalcXingPointWith                    :  x = " << xx.x() << " y = "<< xx.y() << " z = " << xx.z() << " r = " << xx.Perp() << " phi = " << xx.Phi() << " dphi = " <<  phi << std::endl; 
//  
//  streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::CalcXingPointWith :  xdiff = " << xx.x() - xx_n.x() << " ydiff = "<< xx.y() - xx_n.y() << " zdiff = " << xx.z() - xx_n.z() << std::endl;
  
  return (IsOnSurface(xx) ? 1 : 0);
  
  
}

Bool_t ILDSegmentedDiscMeasLayer::IsOnSurface(const TVector3 &xx) const
{
  
  //  std::cout << "IsOnSurface " << std::endl;  
  
  bool onSurface = false ;
  
  TKalMatrix mv = XvToMv(xx);
  
  // check whether the hit lies in the same plane as the surface, here we are resticted to planes perpendicular to z
  
    //  streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::IsOnSurface Xc.Z = " << GetXc().Z() << std::endl;
  
  if (TMath::Abs(xx.Z()-GetXc().Z()) < 1e-4) {
      //    streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::IsOnSurface z passed " << std::endl;
    double r2 = xx.Perp2();
    
    // quick check to see weather the hit lies inside the min max r 
    if(  r2 <= _rmax*_rmax && r2 >= _trap_rmin*_trap_rmin ) { 
        //      streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::IsOnSurface r2 passed " << std::endl;
      
      double phi_point = angular_range_2PI(xx.Phi());
      
      // get the angle in the local system 
      double gamma = angular_range_2PI(phi_point - _start_phi);
      
      // the angle local to the sector
      double local_phi = fmod(gamma, _segment_dphi) - (0.5*_segment_dphi) ;

      double dist_along_centre_line = xx.Perp()*cos(local_phi);
      
      // check if the point projected onto the centre line is within the limits of the trapezoid
      if( dist_along_centre_line > _trap_rmin && dist_along_centre_line < _trap_rmin + _trap_height){
          //        streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::IsOnSurface dist_along_centre_line passed " << std::endl;
        
        double adj = dist_along_centre_line - _trap_rmin ;
        
        double dist_from_centre_line = fabs(xx.Perp()*sin(local_phi));

        double opp = dist_from_centre_line - 0.5*_trap_inner_base_length;
        
        double tan_delta_angle = opp/adj; 
      
        // check if the point is within the angular limits of the trapezoid
        if (tan_delta_angle < _trap_tan_beta) {
            //          streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::IsOnSurface tan_delta_angle passed " << std::endl;
          onSurface = true ;

        }
      }
    }    
  }
  
  return onSurface;
  
}


ILDVTrackHit* ILDSegmentedDiscMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
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
  
  streamlog_out(DEBUG0) << "ILDSegmentedDiscMeasLayer::ConvertLCIOTrkHit ILDPlanarHit created" 
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

double ILDSegmentedDiscMeasLayer::angular_range_2PI( double phi ) const {
  
  //bring phi_point into range 0 < phi < +2PI
  while (phi < 0) {
    phi += 2.0 * M_PI;
  }
  while (phi >= 2.0*M_PI) {
    phi -= 2.0 * M_PI;
  }
  
  return phi;
  
}

