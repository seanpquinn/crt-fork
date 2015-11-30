//////////////////////////////////////////////////////////
// Filename: source_cr.h
// Authors:  Brian Baughman
//           Michael Sutherland
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the cosmic ray source class
//
//////////////////////////////////////////////////////////

#ifndef _CRSRC_H
#define _CRSRC_H

#include "source.h"
#include <fstream>
#include "spectrum_monoE.h"
class CRSrc: public Source
  {
  public:
    CRSrc(double shift, gsl_rng *rn, CRSpectrum *spec, 
          bool bktrk, double Long, double Lat,
          double pM, double pQ):
    Source(1, shift, rn, spec, bktrk, pM, pQ)
    {
      
      distance = SOURCEDISTANCE;
      srcLoc = gvector(Long,Lat);
      srcLoc.norm(distance);
    };
    CRSrc(std::ifstream &ins, gsl_rng *rn, bool bktrk):
    Source()
    {
      rng_ptr=rn;
      backtrack=bktrk;
      double tE, tL, tB;
      ins >> tE >> tL >> tB >> crMass >> crCharge;
      tL*=DEG2RAD;
      tB*=DEG2RAD;
      crMass *= PROTONMASS;
      crCharge *= PROTONCHARGE;
      partVel=gvector(0,0,0);
      CRperSrc=1; 
      crShift=0.;
      spectrum = new MonoSpectrum(tE*EeV);
      distance = SOURCEDISTANCE;
      srcLoc = gvector(tL,tB);
      srcLoc.norm(distance);
    };
    virtual ~CRSrc(){};
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << spectrum->GetTag() << "C " << spectrum->GetDescription() << ' ' 
      << srcLoc.l()/DEG2RAD << ' ' << srcLoc.b()/DEG2RAD << ' ' 
      << crMass/PROTONMASS << ' ' << crCharge/PROTONCHARGE;
      return ff.str();
    }
    
    void getParticle(Particle &particle){
      //generate particle E, lorentz gamma, and Velocity
      //determine starting position
      Source::initParticle(particle);
    }
    
  };
#endif //_CRSRC_H
