//////////////////////////////////////////////////////////
// Filename: source_eg.h
// Authors:  Brian Baughman
//           Michael Sutherland
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the extragalactic point source class
//
//////////////////////////////////////////////////////////

#ifndef _EGSRC_H
#define _EGSRC_H

#include "source.h"
#include <istream>

class EGSrc: public Source
  {
  private:
    double Brand;           // Magnitude of random field
    double Lcorr;           // Correlation length of magnetic field in Kpc
    double Leg;             // Distance away from Galaxy in Kpc
    double deltaprefactor;  // Prefactor to save computation time
    int approxtype;         // Variable holding type of approximation 
    // 0==min(R_larmor)>>Lcorr 1==min(R_larmor)~Lcorr
    // Returns the approxmiate average deflection angle 
    double deltarms(double Z, double E) {
      return deltaprefactor*fabs(Z)/E;
    }
    // Returns the approximate Larmor radius
    double rlarmor(double Z, double E) {
      return 1e-6*E/(Z*Brand);
    }
  public:
    EGSrc(unsigned long np, double shift, gsl_rng *rn, CRSpectrum *spec, 
          bool bktrk, double d, double Long, double Lat,
          double B, double Lc, double Le,
          double pM, double pQ):
    Source(np, shift, rn, spec, bktrk, pM, pQ),Brand(B),Lcorr(Lc),Leg(Le),
    deltaprefactor(0.),approxtype(1)
    {
      if (Brand==0) {
        approxtype=-1;
        deltaprefactor=0.;
      } else {
        double rlar = rlarmor(spectrum->GetMaxE(),crCharge);
        if((rlar-Lcorr)>1e3) {
          approxtype=0;
          deltaprefactor=0.533*DEG2RAD*Brand*1e2*Leg;
        } else {
          approxtype=1;
          deltaprefactor=0.0122*DEG2RAD*Brand*1e2*sqrt(Lcorr*Leg);
        }
      }
      beamAngle = 0.;
      distance = d;
      srcLoc = gvector(Long,Lat);
      
    };
    virtual ~EGSrc(){};
    
    //  this will create a particle, process it
    void getParticle(Particle &particle){ 
      //generate particle E, lorentz gamma, and Velocity
      //determine starting position
      Source::initParticle(particle);
      Source::crshift(particle);
      Source::gaussian_smear(particle,deltarms(particle.GetCharge(),
                                               particle.GetEnergy()));
      
    }
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << spectrum->GetTag() << "E " <<  CRperSrc << ' ' <<    crShift << ' ' 
      <<  spectrum->GetDescription() << ' ' 
      << distance << ' ' 
      << srcLoc.l()/DEG2RAD << ' ' << srcLoc.b()/DEG2RAD << ' ' 
      << Brand << ' ' << Lcorr << ' ' << Leg << ' '
      << crMass/PROTONMASS << ' ' << crCharge/PROTONCHARGE;
      return ff.str();
    }
    
    
  };
#endif //_EGSRC_H

