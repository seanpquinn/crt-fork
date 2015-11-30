//////////////////////////////////////////////////////////
// Filename: detector.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the detector singleton class
//                that provides detector parameters and functions
//
//////////////////////////////////////////////////////////

#ifndef _DETECTOR_H
#define _DETECTOR_H

#include "globals.h"
#include "gvector.h"
#include <map>
#include <string>
#include <sstream>

namespace AugerSouth {
  const double Dec = -35.2*DEG2RAD;
  const double maxZenithEff = 60.*DEG2RAD;
  const double maxExpo = 1.8109155297707715;
};

class Detector
{
 private:
  static Detector* pinstance;
  gvector pD0;
  double detR;
  double tmin;
  double tmax;
  double tstep;
  bool exposure; //0=isotropic, 1=auger
  double detDec;  // Declincation of detector site
  double detFOV;  // Angular field of view of detector
  double detMaxExpo;  // Angular field of view of detector
  std::string dtype;
  std::map<long,long > srccnts;
  short npnts;
  void *rng_ptr;         // pointer to rng
  double zeta(double dec);
  double alpham(double z);
  double relw(double dec);
  
 public:
  Detector(void):pD0(gvector(0.,0.,0.)),detR(0.1),tmin(0.),tmax(1.e9),tstep(10.),exposure(0),
  detDec(AugerSouth::Dec),detFOV(AugerSouth::maxZenithEff),
  detMaxExpo(AugerSouth::maxExpo),
  dtype("default"),npnts(20),rng_ptr(0)
    {};
  virtual ~Detector(){};
  static Detector* Instance();
  void SetP0(double fx,double fy, double fz);
  void SetP0(gvector& fgv);
  gvector GetP0(void);
  double GetPosMag(void);
  void SetDetR(double);
  double GetDetR(void);
  void SetTmax(double);
  double GetTmax(void);
  void SetTmin(double);
  double GetTmin(void);
  void SetTstep(double);
  double GetTstep(void);
  void SetDetDec(double);
  double GetDetDec(void);
  void SetDetFOV(double);
  double GetDetFOV(void);
  void SetDetMaxExpo(double);
  double GetDetMaxExpo(void);
  void SetRndPtr(void * rng);

  void SetType(std::string);
  std::string GetType(void);
  bool GetExposure(void);

  gvector getshift(gvector &lx, gvector &ly, 
                   long srcnum, double frac=0., double angle=0.);

  double exposureValue(const double &ra, const double &dec);
  double exposureValueGal(const double &glong, const double &glat);
  void SampleExposure(double &glong,double &glat);

  std::string GetDescription(void);

  double MaxExpo(void) {
    return detMaxExpo;
  }
};
#endif //_DETECTOR_H

