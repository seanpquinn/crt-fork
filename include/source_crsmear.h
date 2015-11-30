//////////////////////////////////////////////////////////
// Filename: source_crsmear.h
// Authors:  Michael Sutherland
//           Brian Baughman
//
// Copyright 2011 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the cosmic ray source class using energy
//                and angular resolution smearing
//
//////////////////////////////////////////////////////////

#ifndef _CRSRCSMEAR_H
#define _CRSRCSMEAR_H

#include "source.h"
#include <fstream>
#include "spectrum_monoE.h"
class CRSrcSmear: public Source
  {
  protected:
    double Eres;            // Energy resolution
    double Ares;            // Angular resolution

  public:
    CRSrcSmear(unsigned long np, double shift, gsl_rng *rn, CRSpectrum *spec, 
          bool bktrk, double Long, double Lat,
          double pM, double pQ, double eres, double ares):
    Source(np, shift, rn, spec, bktrk, pM, pQ)
    {
      
      distance = SOURCEDISTANCE;
      srcLoc = gvector(Long,Lat);
      srcLoc.norm(distance);
      SetSrcLoc(srcLoc);
      CRperSrc = np;
      Eres = eres;
      Ares = ares;
    };
    CRSrcSmear(std::ifstream &ins, gsl_rng *rn, bool bktrk):
    Source()
    {
      rng_ptr=rn;
      backtrack=bktrk;
      double tE, tL, tB;
      ins >> tE >> tL >> tB >> crMass >> crCharge >> CRperSrc >> Eres >> Ares;
      tL*=DEG2RAD;
      tB*=DEG2RAD;
      crMass *= PROTONMASS;
      crCharge *= PROTONCHARGE;
      partVel=gvector(0,0,0);
      //CRperSrc=1; 
      crShift=0.;
      spectrum = new MonoSpectrum(tE*EeV);
      distance = SOURCEDISTANCE;
      srcLoc = gvector(tL,tB);
      srcLoc.norm(distance);
    };
    virtual ~CRSrcSmear(){};
    
    std::string GetDescription(void){
      std::ostringstream ff;
      ff << spectrum->GetTag() << "CS " << spectrum->GetDescription() << ' ' 
      << srcLoc.l()/DEG2RAD << ' ' << srcLoc.b()/DEG2RAD << ' ' 
      << crMass/PROTONMASS << ' ' << crCharge/PROTONCHARGE << ' '
      << CRperSrc << ' ' << Eres << ' ' << Ares*RAD2DEG;
      return ff.str();
    }
    
    void getParticle(Particle &particle){
      //generate particle E, lorentz gamma, and Velocity
      //determine starting position
      Source::initParticle(particle);
      Source::gaussian_smear(particle, Ares);
      Source::gaussian_smear_energy(particle, Eres);
    }
    
  };
#endif //_CRSRCSMEAR_H
