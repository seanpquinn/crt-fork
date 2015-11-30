//////////////////////////////////////////////////////////
// Filename: source_isotropic.h
// Authors:  Brian Baughman
//           Michael Sutherland
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the isotropic source class
//
//////////////////////////////////////////////////////////

#ifndef _ISOTROPICSRC_H
#define _ISOTROPICSRC_H

#include "source.h"

class IsotropicSrc: public Source
  {
  public:
    IsotropicSrc(unsigned long np, double shift, gsl_rng * (&rn), 
                 CRSpectrum *spec, bool bktrk, double beam,
                 double pM, double pQ, bool exp):
    Source(np, shift, rn, spec, bktrk, pM, pQ)
    {
      beamAngle = beam;
      exposure = exp;
    };
    virtual ~IsotropicSrc(){};
    
    void getParticle(Particle &particle){
      //generate particle E, lorentz gamma, and Velocity
      //determine starting position
      double tL, tB;
      if(exposure) { //use Detector exposure function
        det->SampleExposure(tL,tB);
      } else { // Isotropic
        tL = gsl_ran_flat(rng_ptr,-pi,pi);
        tB = asin(gsl_ran_flat(rng_ptr,-1.,1.));
      }
      srcLoc.setcoords(tL,tB);
      initilized=false;
      Source::initParticle(particle);
      Source::crshift(particle);
    }
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << spectrum->GetTag();
      if(exposure==0)
        ff << "I ";
      else
        ff << "X ";
      
      ff << CRperSrc << ' ' <<    crShift << ' ' 
      <<  spectrum->GetDescription() << ' ' 
      << crMass/PROTONMASS << ' ' << crCharge/PROTONCHARGE;
      
      return ff.str();
      
    }
    
  };


#endif //_ISOTROPICSRC_H
