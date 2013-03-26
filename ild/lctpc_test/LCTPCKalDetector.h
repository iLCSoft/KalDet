#ifndef LCTPCKALDETECTOR_H
#define LCTPCKALDETECTOR_H

#include "TVKalDetector.h"

#include "ILDVMeasLayer.h"

#include <map>

namespace gear{
  class GearMgr ;
}

namespace kaldet{

  /**
   * The LCTPC implementation for a TPC which is completely instantiated from GEAR.
   * 
   */
class LCTPCKalDetector : public TVKalDetector {

public:

    LCTPCKalDetector() {};

    /** 
     * The constructor. All information to initialise the TPC is taken from GEAR.
     *
     * The class has been copied from GearTPCKalDetector class and adopted for the use of MarlinTrk
     * You can find comments and necessary information in the original class
     * 
     */
    LCTPCKalDetector(const gear::GearMgr& gearMgr);

    /// The destructor.
    virtual ~LCTPCKalDetector();

    /**
     * Get access to the measurement layers using moduleID and row.
     * Do not directly access the measurement layers using At() 
     * because the order depends on the order in the gear file.
     * Throws a gear::Exception if the row on the module is not defined.
     */
    virtual ILDVMeasLayer const * GetMeasLayer(int moduleID, int row) const;

protected:
    /// Map which contains the information which measurement layer is stored
    /// at which position in the array.
    std::map< std::pair<int, int >, Int_t > moduleRowToMeasurementLayerMap;
};

}// namespace kaldet
#endif //LCTPCKALDETECTOR_H
