#include "LCTPCCylinderMeasLayer.h"
//#include "GearTPCHit.h"
//#include "GearTPCKalDetector.h"
//#include "GearTPCCylinderHit.h"

//#include "ILDCylinderMeasLayer.h"
//#include "ILDCylinderHit.h"

//#include <EVENT/TrackerHit.h>
//#include <EVENT/TrackerHitZCylinder.h>

#include <TMath.h>

#include <iostream>
#include <streamlog/streamlog.h>

#include <gear/GEAR.h>

namespace kaldet
{


LCTPCCylinderMeasLayer::LCTPCCylinderMeasLayer(
        TMaterial &min,
        TMaterial &mout,
        Double_t   r0,
        Double_t   lhalf,
        Double_t   x0,
        Double_t   y0,
        Double_t   z0,
        Double_t   Bz,
        Bool_t     is_active,
        Int_t      CellID,
        const Char_t  *name,

        Int_t      module,
        Int_t      row,
        Bool_t     isPerfect,
        Double_t   sigmaX0,
        Double_t   sigmaX1,
        Double_t   sigmaZ0,
        Double_t   sigmaZ1,
        Double_t   phiMin,
        Double_t   phiMax )
    : ILDCylinderMeasLayer(min, mout, r0, lhalf, x0, y0, z0, Bz, is_active, CellID, name)

//       : GearTPCMeasLayer(min, mout, module, isPerfect, isActive,
//                          sigmaX0, sigmaX1, sigmaZ0, sigmaZ1),
//         TCylinder(r0, lhalf, xc.X(), xc.Y(), xc.Z()),
//                fPhiMin(phimin),
//                fPhiMax(phimax)
{
  //FIXME: As the handling of cylinder segments is not defined yet we force the layer to be
  //a perfect cylinder
  fIsPerfect = true;

  // for a full cylinder phi min and max do not make sense. Give a warning if they are used.
  if (fIsPerfect)
  {
    if (phiMin!=-TMath::Pi())
    {
      streamlog_out(WARNING) << "LCTPCCylinderMeasLayer: Ignoring input parameter phimin."
			     << " The current implementation is a full cylinder." << std::endl;
      phiMin=-TMath::Pi();
    }
    if (phiMax!=TMath::Pi())
    {
      streamlog_out(WARNING) << "LCTPCCylinderMeasLayer: Ignoring input parameter phimax."
			     << " The current implementation is a full cylinder." << std::endl;
      phiMax=TMath::Pi();
    }
  }
}

LCTPCCylinderMeasLayer::~LCTPCCylinderMeasLayer()
{
}

TKalMatrix LCTPCCylinderMeasLayer::XvToMv(const TVector3 &xv) const
{
  TVector3 xxv = xv - GetXc();

  Double_t phi = TMath::ATan2(xxv.Y(), xxv.X());// - fPhiMin;

  static Double_t kPi    = TMath::Pi();
  static Double_t kTwoPi = 2 * kPi;
  while (phi < -kPi) phi += kTwoPi;
  while (phi >  kPi) phi -= kTwoPi;

  // Calculate hit coordinate information:
  //   mv(0, 0) = r * phi
  //     (1, 0) = drift distance

  TKalMatrix mv(kMdim, 1);
  mv(0, 0) = GetR() * phi;

  mv(1, 0) = xxv.Z();

  return mv;
}

TKalMatrix LCTPCCylinderMeasLayer::XvToMv(const TVTrackHit &vhit,
                                  const TVector3   &xv) const
{
  return XvToMv(xv);
}


//ILDVTrackHit* LCTPCCylinderMeasLayer::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
//
//  if ( ! trkhit) {
//    streamlog_out(ERROR) << "LCTPCCylinderMeasLayer::ConvertLCIOTrkHit trkhit pointer is NULL" << std::endl;
//    return NULL;
//  }
//
//  const TVector3 hit( trkhit->getPosition()[0], trkhit->getPosition()[1], trkhit->getPosition()[2]) ;
//  //SJA:FIXME: this assumes that the cylinder is centred at 0,0
//
//  // convert to layer coordinates
//  TKalMatrix h    = this->XvToMv(hit);
//
//  Double_t  x[2] ;
//  Double_t dx[2] ;
//
//  x[0] = h(0, 0);
//  x[1] = h(1, 0);
//
//
//  EVENT::TrackerHitZCylinder* cylinder_hit = dynamic_cast<EVENT::TrackerHitZCylinder*>( trkhit ) ;
//
//  if(cylinder_hit){
//    std::cout<<"Alex:: this one is called"<<std::endl;
//    // convert errors
//    dx[0] = cylinder_hit->getdRPhi();
//    dx[1] = cylinder_hit->getdZ();
//  }
//  else {
//    std::cout<<"Alex:: the different one is called"<<std::endl;
//    // convert errors
//    dx[0] = sqrt(trkhit->getCovMatrix()[0] + trkhit->getCovMatrix()[2]) ;
//    dx[1] = sqrt(trkhit->getCovMatrix()[5]) ;
//    // Alex::
//    // Some reasonably dummy values for tests
////    dx[0] = 3;
////    dx[1] = 0.01;
//  }
//
//
//  bool hit_on_surface = IsOnSurface(hit);
//
//  streamlog_out(DEBUG1) << "LCTPCCylinderMeasLayer::ConvertLCIOTrkHit ILDCylinderHit created"
//  << " R = " << hit.Perp()
//  << " Layer R = " << this->GetR()
//  << " RPhi = "  <<  x[0]
//  << " Z = "     <<  x[1]
//  << " dRPhi = " << dx[0]
//  << " dZ = "    << dx[1]
//  << " x = " << trkhit->getPosition()[0]
//  << " y = " << trkhit->getPosition()[1]
//  << " z = " << trkhit->getPosition()[2]
//  << " onSurface = " << hit_on_surface
//  << std::endl ;
//
//  return hit_on_surface ? new ILDCylinderHit( *this , x, dx, this->GetBz(), trkhit) : NULL;
//
//}


//TVector3 LCTPCCylinderMeasLayer::HitToXv(const TVTrackHit &vhit) const
//{
//  //  const EXTPCHit &hit = dynamic_cast<const EXTPCHit &>(vhit);
//
//
//  Double_t phi = vhit(0, 0) / GetR();// + fPhiMin;
//
//  Double_t z   = vhit(1, 0);
//
//  Double_t x   = GetR() * TMath::Cos(phi) + GetXc().X();
//  Double_t y   = GetR() * TMath::Sin(phi) + GetXc().Y();
//
//  return TVector3(x, y, z);
//}

//void LCTPCCylinderMeasLayer::CalcDhDa(const TVTrackHit & vhit,
//                              const TVector3   &xxv,
//                              const TKalMatrix &dxphiada,
//                                    TKalMatrix &H) const
//{
//  double vDrift;
//  try{
//    const GearTPCHit &hit = dynamic_cast<const GearTPCHit &>(vhit);
//    vDrift = hit.GetVdrift();
//  }
//  catch(std::bad_cast &)
//  {
//    streamlog_out(ERROR) << "LCTPCCylinderMeasLayer::CalcDhDa :"
//			 << "Cannot cast incoming TVTrackHit to GearTPCHit!"<< std::endl;
//    throw;
//  }
//
//  // Calculate
//  //    H = (@h/@a) = (@phi/@a, @z/@a)^t
//  //  where
//  //        h(a) = (phi, z)^t: expected meas vector
//  //        a = (drho, phi0, kappa, dz, tanl, t0)
//  //
//
//  Int_t sdim = H.GetNcols();
//  Int_t hdim = TMath::Max(5, sdim - 1);
//
//  TVector3 xxvc = xxv - GetXc();
//  Double_t xv   = xxvc.X();
//  Double_t yv   = xxvc.Y();
//  Double_t xxyy = xv * xv + yv * yv;
//
//  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
//
//  for (Int_t i = 0; i < hdim; i++) {
//    H(0, i)  = - (yv / xxyy) * dxphiada(0, i)
//               + (xv / xxyy) * dxphiada(1, i);
//    H(0, i) *= GetR();
//
//    H(1, i)  = dxphiada(2, i);
//  }
//
//  if (sdim == 6) {
//    H(0, sdim - 1) = 0.;
//
//    // KILLENB I don't understand what this does, but I removed vDrift, so I set this to 0.
//    H(1, sdim - 1) = - vDrift;
//  }
//}
//
//Double_t LCTPCCylinderMeasLayer::GetSortingPolicy() const
//{
//   // The sorting policy (copied from the header file):
//  // The layers are first sorted by radius + offset. This offset is only
//  // useful for segments of a cylinder, like the LP1.
//  // As offsets in this case can be positive or negative, but only make sense in one
//  // direction (you need a continuous number), we only allow offsets in X.
//  // This should not be too much of a problem, you should be able to rotate your coordinates
//  // so the offset is in X. If not you have to extend the sorting policy. (Please thake
//  // care not to reduce versatility when doing so. You might want to implement your own class?)
//  //
//  // For equal radii  + offset the layers are sorted by moduleID. As we have to squeeze this
//  // information into only one number, we multiply the radius + offset by 1e9 and add the moduleID.
//  // A double has a precision of 53 bits, which is 15.9 digits. So the radius can be up to 1e6.9 mm
//  // without causing the last digit of the the ModuleID to be cut, and for up to 1000 modules the
//  // layers can be distinguished down to 1 nm without the two numbers mixing, or down to 1 micron
//  // with up to 1.000.000 modules.
//
//  // The additional sorting by module is intended for cylinder segments. Here only one module/row
//  // per layer is allowed, so we just take the first entry in the set. In case of a perfect layer
//  // it does not matter because there should only be one layer at this radius, so the sort order
//  // should not be affected by adding an arbitrary module ID (as long as the module ID is < 1e6, as
//  // described above).
//
//  // check that the offset is onyl in X
//  if ( (GetXc().Y()!=0) || (GetXc().Z()!=0) )
//  {
//    throw gear::NotImplementedException("Offset is only allowed in X in the current implementation");
//  }
//
//  int moduleID = fModuleRows.begin()->first;
//  // give a warning in case of very large module IDs. The sorting policy might have to be adapted
//  if( moduleID > 1e6 )
//  {
//    streamlog_out(WARNING) << "LCTPCCylinderMeasLayer::GetSortingPolicy(): "
//			   << "Very large module ID " << moduleID << " found. "
//			   << "This might compromise the sorting policy."
//			   << std::endl;
//  }
//
//  // FIXME: Do we have to check for a max r-offsetX?
//  return (GetR() + GetXc().X())*1e9 + moduleID;
//}
//
//GearTPCHit * GearTPCCylinderMeasLayer::createHit(Double_t * meas,
//						 Double_t * dmeas,
//						 void * hitPointer,
//						 Double_t bField,
//						 Double_t vDrift,
//						 Int_t           m) const
//{
//  return new GearTPCCylinderHit(*this, meas, dmeas, hitPointer,
//				bField, vDrift, m);
//}


}//namespace kaldet
