//////////////////////////////////////////////////////////
// Filename: spectrum_brokenplaw.h
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

#ifndef _BPLAWSPECTRUM_H
#define _BPLAWSPECTRUM_H

#include "spectrum.h"
#include "globals.h"

class BPLawSpectrum: public CRSpectrum
  {
    
  private:
    double index0;          // first spectral index [minE,ebreak]
    double norm0;           // first normalization
    double prob0;           // probablity for first E range
    double index1;          // second spectral index  (ebreak,maxE]
    double norm1;           // second normalization
    double prob1;           // probablity for second E range
    double ebreak;          // Energy at break
    double ebreaknorm;      // break normalization
    double probnorm;        // probablity normalization
    
  public:
    BPLawSpectrum():CRSpectrum(),index0(0.),norm0(0.),index1(0.),norm1(0.){
      
    };
    BPLawSpectrum(double fminE, double fmaxE, double findex0, double findex1,
    double febreak)
    :CRSpectrum(fminE,fmaxE),index0(findex0),norm0(0),index1(findex1),norm1(0),
    ebreak(febreak){
      if (ebreak<minE) {
        ebreak=fminE;
        minE=febreak;
      } else if (ebreak>maxE) {
        ebreak=fmaxE;
        maxE=febreak;
      }
      ebreaknorm = pow(ebreak,index0-index1);
      norm0 = GetNorm(index0,minE,ebreak);
      prob0 = 1./norm0;
      norm1 = GetNorm(index1,ebreak,maxE);
      prob1 = ebreaknorm/(norm1);
      probnorm = prob0+prob1;
      prob0/=probnorm;
      prob1/=probnorm;
    };
    virtual ~BPLawSpectrum() {
      
    };
    

    
    double GetIndex0(void)
    {
      return index0;
    }
    double GetIndex1(void)
    {
      return index1;
    }    
    double GetNorm(double dex, double low, double high)
    {
      return (1.+dex)/( pow(high,(1.+dex)) - pow(low,(1.+dex)) );
    }
    double GetProb(double n, double dex, double low, double high)
    {
      return n*(pow(high,1.+dex)-pow(low,1.+dex))/(dex+1);
    }

    
    double GetRandomE(double ran)
    {
      if ( minE == maxE )
        return minE;
      if (ran<prob0) 
        return GetProbEne(ran/(prob0),index0,minE,ebreak);
      else 
        return GetProbEne((ran-prob0)/(prob1),index1,ebreak,maxE);
    }
    double GetProbEne(double ran,double dex, double low, double high) {
      double pterm = ran*( pow(high,(1.+dex)) - pow(low,(1.+dex)) );
      pterm += pow(low, (1.+dex));
      return pow(pterm, (1./(1.+dex)));
    }
    double GetRandomProb(double en)
    {
      double num = 0.;
      double denom = 1.;
      if (en<=ebreak) {
        num = pow(en, 1.+index0) - pow(minE, 1.+index0);
        denom = pow(ebreak, 1.+index0) - pow(minE, 1.+index0);
      } else {
        num = ebreaknorm*pow(en, 1.+index1) - pow(minE, 1.+index1);
        denom = pow(maxE, 1.+index1) - pow(ebreak, 1.+index1); 
      }
      return (num/denom);
    }
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << index0 << ' ' << index1 << ' ' << ebreak << ' ' <<   minE/EeV << ' ' << maxE/EeV;
      return ff.str();
    }
    std::string GetTag(void) {
      return "BP";
    }
    
  };

#endif
