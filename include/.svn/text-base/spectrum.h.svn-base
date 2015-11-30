//////////////////////////////////////////////////////////
// Filename: spectrum.h
// Authors:  Michael Sutherland
//           Brian Baughman
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation base
//                of the energy spectrum class
//
//              Energy spectra follow the form:
//                dN/dE ~   E ^ (a)
//
//////////////////////////////////////////////////////////

#ifndef _CRSPECTRUM_H
#define _CRSPECTRUM_H

#include <string>

class CRSpectrum
{

 protected:
  double  minE,maxE;     // max/min energy limits for spectrum
 public:
  CRSpectrum():minE(0.),maxE(0.) {
    
  };
  CRSpectrum(double fminE, double fmaxE):minE(fminE),maxE(fmaxE) {
    if (minE>maxE) {
      minE=fmaxE;
      maxE=fminE;
    }
  };
  virtual ~CRSpectrum() {
    
  };

  void SetMaxE(double max)
  {
    maxE = max;
  }
  
  double GetMaxE(void)
  {
    return maxE;
  }
  
  void SetMinE(double min)
  {
    minE = min;
  }
  
  double GetMinE(void)
  {
    return minE;
  }
  
  // Virtual functions all Spectra need to define
  virtual double GetRandomE(double)=0;
  virtual double GetRandomProb(double)=0;
  virtual std::string GetDescription(void)=0;
  virtual std::string GetTag(void)=0;
};

#endif
