#include "EXEventGen.h"
#include "EXVKalDetector.h"
#include "EXVMeasLayer.h"
#include "TPlane.h"
#include "TRandom.h"

//-----------------------------------
// Track Parameters
//-----------------------------------

#define __DR__     0.
#define __DZ__     0.

ClassImp(EXEventGen)

Double_t EXEventGen::fgT0 = 0.; // [nsec]

THelicalTrack EXEventGen::GenerateHelix(Double_t pt,
                                        Double_t cosmin,
                                        Double_t cosmax,
					Double_t phimin,
					Double_t phimax,
					TVector3 xv0)
{
   // ---------------------------
   //  Generate a helical track
   // ---------------------------

   Double_t dr  = __DR__;
   Double_t fi0 = gRandom->Uniform(phimin,phimax);
   Double_t cpa = 1. / pt;
   Double_t dz  = __DZ__;
   Double_t cs  = gRandom->Uniform(cosmin, cosmax);
   Double_t tnl = cs / TMath::Sqrt((1-cs)*(1+cs)); 
   Double_t x0  = xv0.X();
   Double_t y0  = xv0.Y();
   Double_t z0  = xv0.Z();

   Double_t b   = dynamic_cast<const EXVKalDetector &>
                 (dynamic_cast<EXVMeasLayer *>
                 (fCradlePtr->At(0))->GetParent(kFALSE)).GetBfield();

   return THelicalTrack(dr,fi0,cpa,dz,tnl,x0,y0,z0,b);
}

void EXEventGen::Swim(THelicalTrack &heltrk)
{
   // ---------------------------
   //  Swim track and Make hits
   // ---------------------------

   Double_t dfi       = -dynamic_cast<TVSurface *>(fCradlePtr->At(0))
                            ->GetSortingPolicy()
                         / heltrk.GetRho();

   Int_t    nlayers   = fCradlePtr->GetEntries();
   Int_t    dlyr      = 1;
   Double_t dfisum    = 0.;

   for (Int_t lyr = 0; lyr >= 0; lyr += dlyr) { // loop over layers
      // change direction if it starts looping back
      if (lyr == nlayers - 1) dlyr = -1;

      EXVMeasLayer &ml = *dynamic_cast<EXVMeasLayer *>(fCradlePtr->At(lyr));
      TVSurface    &ms = *dynamic_cast<TVSurface *>(fCradlePtr->At(lyr));
      TVector3 xx;
      Double_t dfis = dfi;
      if (!ms.CalcXingPointWith(heltrk,xx,dfi,1)
       || TMath::Abs(dfi) > TMath::Pi()
       || TMath::Abs(dfi + dfisum) > TMath::TwoPi()) {
         dfi = dfis;
         continue;
      }
      // should use the material behind the surface since dfi is measured 
      // from the last point to the current surface
      // Bool_t   dir    = dlyr < 0 ? kTRUE : kFALSE;

      dfisum += dfi;

      heltrk.MoveTo(xx,dfi);	// move pivot to current hit

      if (ml.IsActive() && dynamic_cast<const EXVKalDetector &>(ml.GetParent(kFALSE)).IsPowerOn()) {
         ml.ProcessHit(xx, *fHitBufPtr); // create hit point
      }
      if (lyr == nlayers - 1) break;
   }
}
