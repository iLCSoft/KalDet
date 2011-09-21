#ifndef __ILDParallelPlanarMeasLayer__
#define __ILDParallelPlanarMeasLayer__
//
//  ILDParallelPlanarMeasLayer.h
//  KalDet
//
//  Created by Steve Aplin on 9/20/11.
//  DESY
//

#include "ILDPlanarMeasLayer.h"

class ILDParallelPlanarMeasLayer : public ILDPlanarMeasLayer {

	Double_t _r;
	Double_t _phi;
	Double_t _cos_phi;
	Double_t _sin_phi;

		
public:
  
	
  ILDParallelPlanarMeasLayer(TMaterial &min,
									TMaterial &mout,
									Double_t   r,
									Double_t   phi,
									Double_t   Bz,
									Double_t   SortingPolicy,
									Double_t   xiwidth,
									Double_t   zetawidth,
									Double_t   xioffset,
									Bool_t     is_active,
									Int_t      layerID = -1,
									const Char_t    *name = "ILDParallelPlanarMeasLayer");

	
	
	virtual ~ILDParallelPlanarMeasLayer() {}
	
	virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
																		 TVector3 &xx,
																		 Double_t &phi,
																		 Int_t     mode,
																		 Double_t  eps = 1.e-8) const;
	
	virtual Int_t    CalcXingPointWith(const TVTrack  &hel,
																		 TVector3 &xx,
																		 Double_t &phi,
																		 Double_t  eps = 1.e-8) const{

		return CalcXingPointWith(hel,xx,phi,0,eps);

	}
	
};

#endif

