//////////////////////////////////////////////////////////
// Filename: source_point.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the point source class
//
//////////////////////////////////////////////////////////

#ifndef _POINTSRC_H
#define _POINTSRC_H

#include <istream>

#include "source.h"

class PointSrc: public Source
  {
  public:
    PointSrc(unsigned long np, double shift, gsl_rng *rn, CRSpectrum *spec, 
             bool bktrk, double d, double Long, double Lat, double beam, 
             double pM, double pQ):
    Source(np, shift, rn, spec, bktrk, pM, pQ)
    {
      beamAngle = beam;
      distance = d;
      srcLoc = gvector(Long,Lat);
      srcLoc.norm(distance);
      SetSrcLoc(srcLoc);
    };
    virtual ~PointSrc(){};
    
    //  this will create a particle, process it
    void getParticle(Particle &particle){ 
      //generate particle E, lorentz gamma, and Velocity
      //determine starting position
      Source::initParticle(particle);
      Source::crshift(particle);
      Source::uniform_smear(particle);
    }
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << spectrum->GetTag() << "P " <<  CRperSrc << ' ' <<    crShift << ' ' 
      <<  spectrum->GetDescription() << ' ' 
      << distance << ' ' 
      << srcLoc.l()/DEG2RAD << ' ' << srcLoc.b()/DEG2RAD << ' ' 
      << beamAngle/DEG2RAD << ' ' 
      << crMass/PROTONMASS << ' ' << crCharge/PROTONCHARGE;
      return ff.str();
    }

  };
#endif //_POINTSRC_H
