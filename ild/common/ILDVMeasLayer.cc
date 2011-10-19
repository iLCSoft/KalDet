
#include "ILDVMeasLayer.h"

Bool_t   ILDVMeasLayer::kActive = kTRUE;
Bool_t   ILDVMeasLayer::kDummy = kFALSE;

//ClassImp(ILDVMeasLayer)

ILDVMeasLayer::ILDVMeasLayer(TMaterial &min,
                             TMaterial &mout,
                             Double_t   Bz,
                             Bool_t     isactive,
                             int        layerID ,
                             const Char_t    *name)  
: TVMeasLayer(min, mout, isactive),
_Bz(Bz),
_layerID( layerID ), 
_name(name)
{
}

ILDVMeasLayer::~ILDVMeasLayer()
{
}
