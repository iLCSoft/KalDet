# v01-14-01

* 2018-01-31 Frank Gaede ([PR#2](https://github.com/iLCSoft/KalDet/pull/2))
  - silence all compiler warnings by using:
      - `"-Wno-effc++ -Wno-unused-parameter"`
      - and SYSTEM includes
  - future users/maintainers should consider to actually **fix these warnings**

# v01-14

# v01-13-03
F. Gaede
- made compatible with c++11 and ROOT6
- removed -ansi -pedantic -Wno-long-long
- added current source dir to include path for dictionary buil

# 01-13-02
- patch release with some additional debug output for diagnostics

# v01-13-01
- fix for ROOT 5.34.18 (TVector3 c'tor) (F.Gaede)

# v01-13
- modified ILDVXDKalDetector:  
- modified to optionally change the relative position of the measurement surface within the sensitive layer => needed for example in the FPCCD
- optionally users can add `<parameter name="relative_position_of_measurement_surface" type="double" value="0.15"/>` to the VXD parameter section in the GEAR file

# v01-12
- added detector geometry for lctpc which is now based on ild classes. Notes:
- it is compatible with Clupatra/MarlinTrk
- the correspondence between row number and module is represented in a different waycompared to the original class GearTPCKalDetector. In the current implementation, row numbers are global, i.e. rows of different modules but with the same center are now the same rows! This is needed to unify the measurement layers that are combinations of several modules.

# v01-11
- ILD
  - final change to FTD and Beam-pipe material

# v01-10
- ILD
  - Reduced the tolerance for the is on surface to 1micron in ILDCylinderMeasLayer.
  - Corrected sign of phi in CalcXingPointWith. ILDDiscMeasLayer, ILDSegmentedDiscMeasLayer, ILDParallelPlanarMeasLayer.
  - Added necessary mode check for determining intersection correctly based on the direction of travel, i.e. fwd or bwd mode =+1 or -1 ILDConeMeasLayer, ILDDiscMeasLayer, ILDPolygonBarrelMeasLayer, ILDSegmentedDiscMeasLayer, ILDParallelPlanarMeasLayer.
  - Tolerance fix for CalcXingPointWith to guard against unphysical values of phi when considering intersection forward or backwards with mod=+/-1.
  - Final changes to the material description to match ILD_oX_XX tracking geometry.

# v01-09 
 - ILD
   - Re-worked MaterialDataBase registerForService will need to be called at least once, triggering initialisation from the gear manager, subsequent calls will check that the gear managers are the same as those which was used to initialise. 
  - VXD support material added. 
  - Optimisation of Debug verbosity.
  - Allow SIT to be pixel based instead of strip. Switch acts on reading strip_angle_deg from GEAR File.

# v01-08
- ILD
  - Added SET detector


# v01-07
- ILD 
  - Added Strip Measurements for Segmented Disc Measurement Layers. 
  - Added proper local coordinates for segments in ILDSegmentedDiscMeasLayer, so that it is no longer restricted to global X and Y. 
- Corrected getIntersectionAndCellID so that it returns the correct module ID, previously it was returning the Nth segment in the segmented disk, so in the case of odd and even disks, it was returning moduleID/2.
- In ILDMeasurementSurfaceStoreFiller FTD now reads from GEAR parameters. Fixed counting problem for sensor. std::vector access out of range.
- Renamed ILDParallelStripPlanarMeasLayer as ILDParallelPlanarStripMeasLayer.
- Cleaned up use of cout and removed unnecessary includes.

# v01-06
- ILD 
  - ILDMeasurementSurfaceStoreFiller created to fill GEAR MeasurementSurfaceStore.
  - Variable number of sensors to the FTD, SIT and SET.
  - Added Rotated Strip measurements for ILDParallelStripPlanarMeasLayer.
- Allow offset in z for planes. 
- Allow to set the origin of the u coordinate which lies in the plane, necessary for a plane composed of two pieces.
- ILDFTDKalDetector adapted to deal with new gear output.
- Added ILDConeMeasLayer using TCutCone. Still needs correcting in ILDSupportKalDetector

# v01-05
- ILD
  - Simple cylinder based SIT design added for loi data
- ILDCylinderMeasLayer, added overriding IsOnSurface method to permit control of the tolerance for the surface comparison.
	    Added use of TrackerHitZCylinder from LCIO.
- ILDDiscMeasLayer, added dedicated CalcXingPointWith method to override the one from TPlane which is based on newtonian method.

# v01-04
-  added pointer to EVENT::TrackerHit in ILDVTrackHit to for better navigation between kaldet and lcio
-  corrected ecal barrel length in support detector


# v01-03
- ILD: 
  - Multilayer measurement layers added, which allow a single kaltest layer to represent numbers detector layers. 
  - Added names to support measurement layers
- ILDParallelStripPlanarMeasLayer and ILDPlanarStripHit First shot at 1-D strip hits. Only measurements purely in the transverse plane are implemeted so far.
- ILDVMeasLayer Name attribute moved up to TVMeasLayer in KalTest, added new method getIntersectionAndCellID. This enables multilayers to efficiently return the crossing point and cellID of the intersected sub layer.
- ILDSupportKalDetector added calo face to support detector, added beryllium and corrected rad length of aluminium
- ILDVXDKalDetector corrected calculation of overlap region
- ILDSITKalDetector added support for double layer SIT with layers facing away from the IP
- ILDTPCKalDetector corrected TPC field cage for new TPC10.cc from Mokka
- ILDFTDKalDetector corrected sign issue and corrected sortpolicy, fixed bug in material, carbon vs air, corrected measurement layer order
- ILDSegmentedDiscMeasLayer Added segmented disc needed for FTD
- Removed unnecessary casting


# v01-02
- Updated for ILD detector reconstruction with KalTest and MarlinTrk
- ILD detectors added in new ild directory. 
- VXD, TPC, FTD(both simple disk and petal based)
- SIT is currently only a placeholder using simple planar design
- common hit and plane classes added to ild/common subdirectory
- uses UTIL/ILDConf.h from lcio for detector element encoding
- Added ILDParallelPlanarMeasLayer for VXD and SIT, which has closed solution for calculating the crossing point of a helix. Constructor ensures that it is parallel by supplying only r and phi, which are used to calculate the normal and centre for planar surface.


# v01-01-01
 - bug fixed in cmake files for mac osx compatibility- no code changes.

# v01-01
 - improvements in cmake files- no code changes.

# v01-00
- first release as part of iLCSoft- no code changes.(except header in subdirectory 'include/kaldet', when installed with cmake)


