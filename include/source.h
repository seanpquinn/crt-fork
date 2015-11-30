//////////////////////////////////////////////////////////
// Filename: source.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation base
//                of the different source classes
//
//////////////////////////////////////////////////////////

#ifndef _SOURCE_H
#define _SOURCE_H

#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include "globals.h"
#include "spectrum.h"
#include "particle.h"
#include "gvector.h"
#include "detector.h"
class Source
  {
    
    
  protected:
    int src_number;           // Source number
    double distance;          // Distance to source in kpc
    gvector srcLoc;           // vector of source location
    gvector partVel;          // vector to store velocity
    double crShift;           // Maximum shift of position from srcLoc
    unsigned long CRperSrc;   // Number of CRs injected by each source
    CRSpectrum* spectrum;     // Spectrum used to generate particles
    double beamAngle;         // [rad], used only for galactic 'pointsrc'
    double relFlux;           // Relative flux of source (not used)
    gsl_rng *rng_ptr;         // pointer to rng
    double y_rk4[N];	        // vector of y functions 
    bool backtrack;           // Flag for back tracking
    double crMass;            // Mass of particles emitted
    double crCharge;          // charge of particles emitted
    Particle p;               // base particle characteristics
    gvector lx;               // "x" unit vector around fixedsrcLoc
    gvector ly;               // "y" unit vector around fixedsrcLoc
    bool initilized;          // boolean which keeps track of initilzation state
    Detector* det;
    bool exposure;
    
    void initParticle(Particle &particle){
      initSource();  // Make sure our source is ready
      particle=p;
      particle.SetEnergy( spectrum->GetRandomE( gsl_ran_flat(rng_ptr,0,1) ) );
      gvector tv0 = partVel;
      tv0.norm(particle.GetVelMag());
      particle.SetV0(tv0);
    }
    void initSource(void) {
      if (initilized) 
        return;      
      if(backtrack){
        srcLoc.norm(1.);
        partVel = srcLoc;
        p.SetCharge(-1.*crCharge);
        p.SetHit(true);
        p.SetDOCA(50.);
        crShift=0;
        p.SetP0( det->GetP0().x(),det->GetP0().y(),det->GetP0().z() );
      }
      else{
        srcLoc.norm(distance);
        partVel = -1.*srcLoc;
        p.SetCharge(crCharge);
        p.SetHit(false);
        p.SetDOCA(9999999.);
        p.SetP0(srcLoc);
        
      }
      for(int i=0;i<N;++i)
        y_rk4[i]=0.;
      
      partVel.norm(1.);
      p.SetMass(crMass);
      SetOrtho();  // Make sure orthogonal vectors set properly.
      initilized=true;
    }
    
    void crshift(Particle &particle){
      //if backtracking, account for injection disk size
      static gvector pstart;
      pstart = p.GetP0();
      static double discFrac, injDiscAngle;
      discFrac = sqrt(gsl_ran_flat(rng_ptr,0.,1.));
      injDiscAngle = gsl_ran_flat(rng_ptr,-pi,pi);
      if(backtrack && det->GetDetR()>0 ){
        pstart += det->getshift(lx, ly, src_number, discFrac, injDiscAngle);
      } 
      //deal with forward-tracking crshift
      else if(crShift > 0){
        discFrac *= crShift;
        pstart += (lx*sin(injDiscAngle)+ly*cos(injDiscAngle))*discFrac;
      }
      particle.SetP0(pstart);
    }
    
    void uniform_smear(Particle &particle){
      if (beamAngle>0) {
        double phi = gsl_ran_flat(rng_ptr,-pi,pi);
        double theta = gsl_ran_flat(rng_ptr,0,beamAngle);
        gvector vstart = (lx*sin(phi)+ly*cos(phi))*sin(theta)+partVel*cos(theta);
        vstart.norm(particle.GetVelMag());
        particle.SetV0(vstart);
      }
    }
    
    void gaussian_smear(Particle &particle, double drms) {
      if (drms>0) {
        double phi = gsl_ran_flat(rng_ptr,-pi,pi);
        double theta = gsl_ran_gaussian_tail(rng_ptr,0,drms);
        while (theta>pi*.5) {  // This must be done to limit possible loop arounds
          theta = gsl_ran_gaussian_tail(rng_ptr,0,drms);
        }
        gvector vstart = (lx*sin(phi)+ly*cos(phi))*sin(theta)+partVel*cos(theta);
        vstart.norm(particle.GetVelMag());
        particle.SetV0(vstart);
      }
    }

    void gaussian_smear_energy(Particle &particle, double drms) {
      if (drms>0) {
        double sign = gsl_ran_flat(rng_ptr,-1.,1.);
        sign /= fabs(sign);
        particle.SetEnergy( particle.GetEnergy()* (1. + sign*gsl_ran_gaussian_tail(rng_ptr,0,drms)) );
        while( particle.GetEnergy() < 0 ){
          sign = gsl_ran_flat(rng_ptr,-1.,1.);
          sign /= fabs(sign);
          particle.SetEnergy( particle.GetEnergy() + particle.GetEnergy()*sign*gsl_ran_gaussian_tail(rng_ptr,0,drms) );
        }
      }
    }
    
  public:
    Source(void):
    src_number(-1),distance(SOURCEDISTANCE),srcLoc(gvector(SOURCEDISTANCE,0,0)),
    partVel(gvector(0,0,0)),crShift(0),CRperSrc(0), 
    spectrum(0),beamAngle(0),relFlux(1),rng_ptr(0),
    backtrack(false),crMass(0),crCharge(0),
    p(Particle()),lx(gvector()),ly(gvector()),
    initilized(false),det(0),exposure(0)
    {
    };
    Source(unsigned long np,double shift,gsl_rng *rn,CRSpectrum* spec,
           bool bktrk, double pM, double pQ):
    src_number(-1),distance(SOURCEDISTANCE),srcLoc(gvector(SOURCEDISTANCE,0,0)),
    partVel(gvector(0,0,0)),crShift(shift),CRperSrc(np), 
    spectrum(spec),beamAngle(0),relFlux(1),rng_ptr(rn),
    backtrack(bktrk),crMass(pM),crCharge(pQ),
    p(Particle()),lx(gvector()),ly(gvector()),
    initilized(false),exposure(0)
    {
      
    };
    virtual ~Source() {
    };
    
    void SetSrcNumber(unsigned int n) {
      src_number=n;
    }
    unsigned int GetSrcNumber(void) {
      return src_number;
    }
    // General functions commmon to all sources
    unsigned long GetNumCRperSrc(void) const
    {
      return CRperSrc;
    }
    
    gvector GetSrcLoc(void) const
    {
      return srcLoc;
    }
    void SetSrcLoc(gvector& fsrcLoc) {
      if (srcLoc!=fsrcLoc) {
        srcLoc=fsrcLoc;
      }
    }
    void SetOrtho(void) {
      lx = gvector(-sin(partVel.l()),
                   cos(partVel.l()),
                   0.);
      lx.norm(1.);
      ly = gvector(-cos(partVel.l())*sin(partVel.b()),
                   -sin(partVel.l())*sin(partVel.b()),
                   cos(partVel.b()));
      ly.norm(1.);
    }
    gvector GetPartVel(void) const
    {
      return partVel;
    }
    void SetPartVel(gvector& fpartVel) {
      partVel=fpartVel;
    }
    
    
    void SetNumCRperSrc(unsigned long n)
    {
      CRperSrc = n;
    }
    
    
    void SetCRshift(double d)
    {
      crShift = d;
    }
    
    double GetCRshift(void) const
    {
      return crShift;
    }
    
    double GetPmass(void){
      return crMass;
    }
    
    double GetPcharge(void){
      return crCharge;
    }
    
    bool GetBktrkStatus(void){
      return backtrack;
    }
    
    void SetDetector(Detector* d){
      det = d;
    }
    
    // Virtual Functions each source needs to define    
    virtual void getParticle(Particle &particle) =0;
    virtual std::string GetDescription(void)=0;
    
  };


#endif
