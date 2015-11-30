//////////////////////////////////////////////////////////
// Filename: source_longitudinal.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the longitudinal source class
//
//////////////////////////////////////////////////////////

#ifndef _LONGSRC_H
#define _LONGSRC_H

#include "source.h"

class LongSrc: public Source
  {
    
  private:
    long ngenerated;
  public:
    LongSrc(unsigned long np, double shift, gsl_rng * (&rn), CRSpectrum *spec, 
            double Long, bool bktrk, double beam, double pM, double pQ, bool exp):
    Source(np, shift, rn, spec, bktrk, pM, pQ),ngenerated(0)
    {
      srcLoc = gvector(Long,0.);
      if(bktrk)
        srcLoc.norm(1);
      else
        srcLoc.norm(SOURCEDISTANCE);
      beamAngle = beam;
      SetSrcLoc(srcLoc);
      exposure = exp;
    };
    virtual ~LongSrc(){};
    
    void getParticle(Particle &particle){
      //generate particle E, lorentz gamma, and velocity
      //determine starting position
      double tL, tB;
      tL = srcLoc.l();
      tB = asin(gsl_ran_flat(rng_ptr,-1.,1.));
      if(exposure) { //use Detector exposure function
        double randN = gsl_ran_flat(rng_ptr,0,det->MaxExpo());
        while(randN > det->exposureValueGal(tL,tB)) {
          randN = gsl_ran_flat(rng_ptr,0,det->MaxExpo());
          tB = asin(gsl_ran_flat(rng_ptr,-1.,1.));
        }
      }
      srcLoc.setcoords(tL,tB);
      initilized=false;
      Source::initParticle(particle);
      ngenerated++;
    }
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << spectrum->GetTag() ;
      if(exposure==0)
        ff << "L ";
      else
        ff << "Y ";
      
      ff <<  CRperSrc << ' ' <<    crShift << ' ' 
      << spectrum->GetDescription() << ' ' 
      << srcLoc.l()/DEG2RAD << ' ' 
      << crMass/PROTONMASS << ' ' << crCharge/PROTONCHARGE;
      
      return ff.str();
      
    }
    
  };


#endif //_LONGSRC_H
