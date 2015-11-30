//////////////////////////////////////////////////////////
// Filename: spectrum_monoplaw.h
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

#ifndef _PLAWSPECTRUM_H
#define _PLAWSPECTRUM_H

#include "spectrum.h"
#include "globals.h"

class PLawSpectrum: public CRSpectrum
  {
    
  private:
    double index;          // spectral index
    double norm;           // normalization
    
  public:
    PLawSpectrum():CRSpectrum(),index(0.),norm(0.){
      
    };
    PLawSpectrum(double fminE, double fmaxE, double findex)
    :CRSpectrum(fminE,fmaxE),index(findex){
      SetNorm();
    };
    virtual ~PLawSpectrum() {
      
    };
    
    void SetIndex(double k)
    {
      index = k;
    }
    
    double GetIndex(void)
    {
      return index;
    }
    
    void SetNorm(void)
    {
      norm = (1.+index)/( pow(maxE,(1.+index)) - pow(minE,(1.+index)) );
    }
    
    double GetNorm(void)
    {
      return norm;
    }
    
    double GetRandomE(double ran)
    {
      if ( minE == maxE )
        return minE;
      
      double pterm = ran*( pow(maxE,(1.+index)) - pow(minE,(1.+index)) );
      pterm += pow(minE, (1.+index));
      return pow(pterm, (1./(1.+index)));
    }
    
    double GetRandomProb(double en)
    {
      double num = pow(en, 1.+index) - pow(minE, 1.+index);
      double denom = pow(maxE, 1.+index) - pow(minE, 1.+index);
      return (num/denom);
    }
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << index << ' ' <<   minE/EeV << ' ' << maxE/EeV;
      return ff.str();
    }
    std::string GetTag(void) {
      return "";
    }
    
  };

#endif
