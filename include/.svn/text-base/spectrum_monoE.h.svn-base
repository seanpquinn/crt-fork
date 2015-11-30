//////////////////////////////////////////////////////////
// Filename: spectrum_monoE.h
// Authors:  Brian Baughman
//           Michael Sutherland
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


#ifndef _MONOSPECTRUM_H
#define _MONOSPECTRUM_H

#include "spectrum.h"
#include "globals.h"

class MonoSpectrum: public CRSpectrum
  {
  public:
    MonoSpectrum():CRSpectrum(){
      
    };
    MonoSpectrum(double fE):CRSpectrum(fE,fE){
    };
    virtual ~MonoSpectrum() {
      
    };
    
    double GetRandomE(double r) {
      return minE;
    };
    
    double GetRandomProb(double e) {
      if (e==minE)
        return 1.;
      else
        return 0.;
    };

    std::string GetDescription(void) {
      std::ostringstream ff;
      ff << minE/EeV ;
      return ff.str();
    };
    std::string GetTag(void) {
      return "";
    }
  };

#endif
