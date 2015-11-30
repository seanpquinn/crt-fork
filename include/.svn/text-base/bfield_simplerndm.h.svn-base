//////////////////////////////////////////////////////////
// Filename: bfield_simplerndm.h
// Authors:  Brian Baughman
//           Michael Sutherland
//
// Copyright 2009 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the simple random magnetic field
//
//////////////////////////////////////////////////////////

#ifndef _SimpleRndmField_H
#define _SimpleRndmField_H

#include "bfield.h"
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

class SimpleRndmField: public BFIELD
  {
    
  private:
    double corrlen;    // Correlation lenth to use
    double sigmalen;   // Sigma for cell size variation
    double sigmanorm;  // Sigma for normalization variation
    double curlen;     // Cell "size" for current location
    gvector curfield;  // Field for current location
    gvector curedge;   // Entering edge of current cell
    gvector curloc;    // Current location
    gvector curvel;    // Current velocity
    gsl_rng *rng_ptr;  // pointer to rng

  public:
    SimpleRndmField():BFIELD(),corrlen(0.),sigmalen(0.),sigmanorm(0.),
    curlen(0.),curfield(gvector(0,0,0)),curedge(gvector(0,0,0)),
    curloc(gvector(0,0,0)),curvel(gvector(0,0,0)),rng_ptr(0) {};
    
    SimpleRndmField(gsl_rng *rn):BFIELD(),corrlen(0.),sigmalen(0.),sigmanorm(0.),
    curlen(0.),curfield(gvector(0,0,0)),curedge(gvector(0,0,0)),
    curloc(gvector(0,0,0)),curvel(gvector(0,0,0)),rng_ptr(rn) {};

    virtual ~SimpleRndmField() {};
    
    void Init(void) {
      if (rng_ptr==0)
        SetRand();
      do {
        curlen = corrlen - gsl_ran_gaussian(rng_ptr,sigmalen);
      } while (curlen<0);
      curedge=curloc+curvel*.5*curlen;
      curfield.setcoords(gsl_ran_flat(rng_ptr,-pi,pi),
                         asin(gsl_ran_flat(rng_ptr,-1.,1.)));
      double curnorm = norm - gsl_ran_gaussian(rng_ptr,sigmanorm);
      do {
        curnorm = norm - gsl_ran_gaussian(rng_ptr,sigmanorm);
      } while (curnorm<0);
      curfield.norm(curnorm);
    }
    
    void SetRand(void) {
      rng_ptr = gsl_rng_alloc(gsl_rng_ranlux389);
    }
    
    void SetCorrLength(const double& cl)
    {
      corrlen = cl;
    }
    
    double GetCorrLength(void) const
    {
      return corrlen;
    }
    void SetSigmaLen(const double& sl)
    {
      sigmalen = sl;
    }
    
    double GetSigmaLen(void) const 
    {
      return sigmalen;
    }   
    void SetSigmaNorm(const double& sn)
    {
      sigmanorm = sn;
    }
    double GetSigmaNorm(void) const 
    {
      return sigmanorm;
    }  
    double Bx(const double& fx, const double& fy, const double& fz)
    {
      SetField(fx,fy,fz);
      return curfield.x();
    }
    
    double By(const double& fx, const double& fy, const double& fz)
    {
      SetField(fx,fy,fz);
      return curfield.y();
    }
    
    double Bz(const double& fx, const double& fy, const double& fz)
    {
      SetField(fx,fy,fz);
      return curfield.z();
    }
    void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
      SetField(fx,fy,fz);
      bvec[0]=curfield.x();
      bvec[1]=curfield.y();
      bvec[2]=curfield.z();
    }
    
    std::string GetType(void) const 
    {
      return "simprand";
    }
    
    std::string GetDescription(void) const {
      std::ostringstream ff;
      ff << "F " <<  GetType()  << ' ' << norm << ' '
      << corrlen << ' ' << sigmalen << ' ' << sigmanorm;
      return ff.str();
    }
    
    unsigned int checkStep(const gvector& P,const gvector& V, double &step)  {
      gvector tpoint = P;
      tpoint=tpoint+step*V;
      tpoint=tpoint-curedge;
      if (tpoint.r()<step)
        step=tpoint.r();
      curvel.setcoords(V.x(),V.y(),V.z());
      return 0;
    }
    
  private:
    void SetField(const double& x, const double& y, const double& z) {
      curloc.setcoords(x,y,z);
      if ((curedge-curloc).r()>=curlen) {
        Init();
      }
    }

  };

#endif
