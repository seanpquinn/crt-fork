//////////////////////////////////////////////////////////
// Filename: source_latitudinal.h
// Authors:  Michael Sutherland
//           Brian Baughman
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the latitudinal source class
//
//////////////////////////////////////////////////////////

#ifndef _LATSRC_H
#define _LATSRC_H

#include "source.h"


class LatSrc: public Source
  {
    
  private:
    long ngenerated;
  public:
    LatSrc(unsigned long np, double shift, gsl_rng * (&rn), CRSpectrum *spec, 
           double Lat, bool bktrk, double beam, double pM, double pQ, bool exp):
    Source(np, shift, rn, spec, bktrk, pM, pQ),ngenerated(0)
    {
      srcLoc = gvector(0.,Lat);
      if(bktrk)
        srcLoc.norm(1);
      else
        srcLoc.norm(SOURCEDISTANCE);
      beamAngle = beam;
      SetSrcLoc(srcLoc);
      exposure = exp;
    };
    virtual ~LatSrc(){};
    
    void getParticle(Particle &particle){
      //generate particle E, lorentz gamma, and Velocity
      //determine starting position
      double tL, tB;
      tL = gsl_ran_flat(rng_ptr,-pi,pi);
      tB = srcLoc.b();
      if(exposure) { //use Detector exposure function
        double randN = gsl_ran_flat(rng_ptr,0,det->MaxExpo());
        while(randN > det->exposureValueGal(tL,tB)) {
          randN = gsl_ran_flat(rng_ptr,0,det->MaxExpo());
          tL = gsl_ran_flat(rng_ptr,-pi,pi);
        }
      }
      srcLoc.setcoords(tL,tB);
      initilized=false;
      Source::initParticle(particle);
      ngenerated++;
    }
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << spectrum->GetTag();
      if(exposure==0)
        ff << "B ";
      else
        ff << "Z ";
      
      ff <<  CRperSrc << ' ' <<    crShift << ' ' 
      <<  spectrum->GetDescription() << ' ' 
      << srcLoc.b()/DEG2RAD << ' '
      << crMass/PROTONMASS << ' ' << crCharge/PROTONCHARGE;
      
      return ff.str();
      
    }
    
  };


#endif //_LONGSRC_H
