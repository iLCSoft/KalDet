#ifndef __ILDVMEASLAYER__
#define __ILDVMEASLAYER__
//*************************************************************************
//* ===================
//*  ILDVMeasLayer Class
//* ===================
//*
//* (Description)
//*   Virtual measurement layer class used by ILD[X]MeasLayer Classes.
//* (Requires)
//* (Provides)
//*     class ILDVMeasLayer
//*
//*************************************************************************
//
#include "TVector3.h"
#include "kaltest/TKalMatrix.h"
#include "kaltest/TCylinder.h"
#include "kaltest/TVMeasLayer.h"
#include "kaltest/TAttDrawable.h"
#include "kaltest/KalTrackDim.h"
#include "TString.h"

class TVTrackHit;
class TNode;

class ILDVMeasLayer : public TVMeasLayer {
public:

  static Bool_t kActive;
  static Bool_t kDummy;

  // Ctors and Dtor
 
  virtual ~ILDVMeasLayer();
  
  inline int getLayerID() const { return _layerID ; } 
  
  inline TString GetMLName() const { return _name;    }
  
  inline Double_t GetBz() const { return _Bz; } 

 protected:
  
  ILDVMeasLayer(TMaterial &min,
		TMaterial &mout,
		Double_t  Bz,
		Bool_t    is_active = ILDVMeasLayer::kActive,
		int layerID = -1 , 
		const Char_t    *name = "ILDMeasL");

  Double_t _Bz ;       // Magnitude of B-Field in Z
  int _layerID ;
  TString  _name;      // layer name
  
 private:


};

#endif
