//////////////////////////////////////////////////////////
// Filename: source_arbtraj.h
// Authors:  Brian Baughman
//           Michael Sutherland
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the arbitrary trajectory source class
//
//////////////////////////////////////////////////////////

#ifndef _ARBTRAJSRC_H
#define _ARBTRAJSRC_H

#include "source.h"
class ArbTrajSrc: public Source
  {
  protected:
    gvector srcVel;
  public:
    ArbTrajSrc(gsl_rng *rn, CRSpectrum *spec, 
          bool bktrk, 
          double px, double py, double pz,
          double vx, double vy, double vz,
          double pM, double pQ):
    Source(1, 0., rn, spec, bktrk, pM, pQ)
    {
      
      distance = sqrt(px*px+py*py+pz*pz);
      srcLoc = gvector(px,py,pz);
      srcVel = gvector(vx,vy,vz);
      partVel = srcVel;
      SetSrcLoc(srcLoc);
    };

    virtual ~ArbTrajSrc(){};
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << spectrum->GetTag() << "A " << spectrum->GetDescription() << ' '
      << srcLoc.x() << ' ' << srcLoc.y() << ' ' << srcLoc.z() << ' '
      << srcVel.x() << ' ' << srcVel.y() << ' ' << srcVel.z() << ' '
      << crMass/PROTONMASS << ' ' << crCharge/PROTONCHARGE;
      return ff.str();
    }
    
    void getParticle(Particle &particle){
      //generate particle E, lorentz gamma, and Velocity
      //determine starting position
      Source::initParticle(particle);
      srcVel.norm(particle.GetVelMag());
      particle.SetV0(srcVel);
      partVel=srcVel;
    }
    
  };
#endif //_ARBTRAJSRC_H
