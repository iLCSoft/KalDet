#ifndef __ILDSEGMENTEDDISCMEASLAYER_H__
#define __ILDSEGMENTEDDISCMEASLAYER_H__

//
//  ILDSegmentedDiscMeasLayer.h
//  KalDet
//
//  Created by Steve Aplin on 11/5/11.
//*************************************************************************
//* ===================
//*  ILDSegmentedDiscMeasLayer Class
//* ===================
//*
//* (Description)
//*   Segemented Disk Planar measurement layer class used with ILDPLanarTrackHit.
//*   Segments are isosolese trapezoids whose axis of symmetry points to the origin
//* (Requires)
//*   ILDVMeasLayer
//* (Provides)
//*     class ILDSegmentedDiscMeasLayer
//*
//*************************************************************************


#include "TVector3.h"
#include "TKalMatrix.h"
#include "TPlane.h"
#include "ILDVMeasLayer.h"
#include "KalTrackDim.h"
#include "TMath.h"
#include <sstream>

class TVTrackHit;

class ILDSegmentedDiscMeasLayer : public ILDVMeasLayer, public TPlane {
public:
  // Ctors and Dtor
  
  ILDSegmentedDiscMeasLayer(TMaterial &min,
                            TMaterial &mout,
                            double   Bz,
                            double   SortingPolicy,
                            int      nsegments,
                            double   zpos,
                            double   phi0, // defined by the axis of symmerty of the first petal
                            double   trap_rmin,
                            double   trap_height,
                            double   trap_innerBaseLength,
                            double   trap_outerBaseLength,
                            bool     is_active,
                            int      layerID = -1,
                            const Char_t    *name = "ILDDiscMeasL");
  
  
  // Parrent's pure virtuals that must be implemented
  
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                const TVector3   &xv) const;
  
  virtual TKalMatrix XvToMv    (const TVector3   &xv) const;
  
  virtual TVector3   HitToXv   (const TVTrackHit &ht) const;
  
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)  const;
  
  virtual ILDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const ;
  
  inline virtual Bool_t   IsOnSurface (const TVector3 &xx) const;
  
  double GetSortingPolicy() const { return _sortingPolicy; }
  
private:
  
  double angular_range_2PI( double phi ) const;
  
  double _sortingPolicy;
  int    _nsegments;
  double _trap_rmin;
  double _trap_height;
  double _trap_inner_base_length;
  double _trap_outer_base_length;
  double _trap_tan_beta; // tan of the openning angle of the petal

  double _rmax;
  double _start_phi; // trailing edge of the first sector
  double _segment_dphi;
  
};



#endif