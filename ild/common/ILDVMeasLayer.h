#ifndef __ILDVMEASLAYER__
#define __ILDVMEASLAYER__

/** ILDVMeasLayer: Virtual measurement layer class used by ILD[X]MeasLayer Classes.
 *
 * @author S.Aplin DESY
 */

#include "TVector3.h"
#include "kaltest/TKalMatrix.h"
#include "kaltest/TCylinder.h"
#include "kaltest/TVMeasLayer.h"
#include "kaltest/TAttDrawable.h"
#include "kaltest/KalTrackDim.h"
#include "TString.h"

#include <vector>

class TVTrackHit;
class TNode;
class ILDVTrackHit;

namespace EVENT{
  class TrackerHit;
}

class ILDVMeasLayer : public TVMeasLayer {
public:
  
  static Bool_t kActive;
  static Bool_t kDummy;
  
  /** Get the layer ID */
  inline int getLayerID() const { return _layerID ; } 
  
  /** Get the Cell ID associated with this measurement layer */
  inline const std::vector<int>& getCellIDs() const { return _cellIDs ; }
  
  /** Get the number of Cell ID associated with this measurement layer */
  inline unsigned int getNCellIDs() const { return _cellIDs.size() ; }
  
  /** Get the Measurement Layers Name */
  inline TString GetMLName() const { return _name;    }
  
  /** Get the Magnetic field at the measurement surface */
  inline Double_t GetBz() const { return _Bz; } 
  
  /** Convert LCIO Tracker Hit to an ILDPLanarTrackHit  */
  virtual ILDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const = 0 ;
  
  /** Check whether the measurement layer represents a series of detector elements */
  bool isMultilayer() { return _isMultiLayer; } 
  
protected:
  
  ILDVMeasLayer(TMaterial &min,
                TMaterial &mout,
                Double_t  Bz,
                Bool_t    is_active = ILDVMeasLayer::kActive,
                int CellID = -1 , 
                const Char_t    *name = "ILDMeasL");
  
  ILDVMeasLayer(TMaterial &min,
                TMaterial &mout,
                Double_t  Bz,
                const std::vector<int>& cellIDs,
                Bool_t    is_active = ILDVMeasLayer::kActive,
                const Char_t    *name = "ILDMeasL");
  
  
  
  Double_t _Bz ;       // Magnitude of B-Field in Z
  int _layerID ;
  std::vector<int> _cellIDs ;

  TString  _name;      // layer name
  bool _isMultiLayer;
  
private:
  
  
};

#endif
