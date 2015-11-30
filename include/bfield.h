//////////////////////////////////////////////////////////
// Filename: bfield.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation base
//                of all magnetic field models
//
//////////////////////////////////////////////////////////

#ifndef _BFIELD_H
#define _BFIELD_H


#include <string>
#include <sstream>

#include "particle.h"
#include <cmath>
#include "globals.h"
#include "gvector.h"

class BFIELD
{
protected:
  // Temp vars to save on allocation time.
  // Be EXTREMELY careful using these as they can be set within functions
  // DO NOT USE THEM unless you check all called fucntions within the function
  // you are working on.
  double tx; 
  double ty;
  double tz;
  double tr;
  double tprefactor;
  double tphi;
  double ttheta;
  double txcomp;
  double tycomp;
  double tzcomp;
  double tstheta;
  double tctheta;
  double tsphi;
  double tcphi;
  double cvec[3];
  double tmag;
  double tHmag;
  double tcap;
  double tdenom;
  double tA;
  double tsp;
  double tcp;
  double retval;
  double ttempscale;
  double tbeta;
protected:
  double norm;   // [microG]
  double maxR;   // [kpc]
  int Nstps; // max number of steps
  std::string descript;
public:
  BFIELD():norm(0.),maxR(20.),Nstps(1000000000),descript(""),
  gc(gvector(solarDistance,0.,0.)) {
    
  };
  virtual ~BFIELD(){
    
  };
  
  // This function should be overloaded if you want to make your field check
  // the next step size.
  // P is the current position of the particle
  // V is the current velocity
  // step is the next proposed step
  // set step to what you want and return 1 if you changed it
  virtual unsigned int checkStep(const gvector P,const gvector V, double &step) const {
    return 0;
  }
  
  int GetNstps(void) const 
  {
    return Nstps;
  }
  void SetNstps(const int& n)
  {
    Nstps = n;
  }
  
  void SetNorm(const double& a)
  {
    norm = a;
  }
  
  double GetNorm(void) const 
  {
    return norm;
  }
  
  void SetMaxR(const double& r){
    maxR = r;
  }
  
  double GetMaxR(void) const {
    return maxR;
  }
  
  double fieldR(const double& x, const double& y, const double& z)
  {
    return sqrt(x*x + y*y + z*z);
  }
  
  double fieldRho(const double& x, const double& y)
  {
    return sqrt(x*x + y*y);
  }
  
  double fieldPhi(const double& x, const double& y)
  {
    return -1.*atan2(y,x);
  }
  double fieldTheta(const double& x, const double& y, const double& z)
  {
    return asin(z/fieldR(x,y,z));
  }
  
  void gal2gc(double &x, double &y, double &z)
  {
    x = -1.*(x-solarDistance);
    y = -1.*y;
  }
  void gc2gal(double &x, double &y, double &z)
  {
    x = -1.*(x+solarDistance);
    y = -1.*y;
  }  
  bool withinField(const double& x, const double& y, const double& z)
  {
    if(fieldR(x,y,z) < maxR)
      return true;
    else
      return false;
  }
  
  // Returns a string of B magnitudes 
  std::string Bvecstr(const double& x, const double& y, const double& z) {
    Bvec(x,y,z,cvec);
    std::ostringstream ff;
    for(int i=0;i<3;++i) {
      ff << cvec[i];
      if (i!=2)
        ff << ',';
    } 
    return ff.str();
  }
  
  // Returns a string of B magnitudes 
  gvector Bgvec(const double& x, const double& y, const double& z) {
    Bvec(x,y,z,cvec);
    return gvector(cvec[0],cvec[1],cvec[2]);
  }
  
  virtual double Bx(const double& fx, const double& fy, const double& fz)=0;
  virtual double By(const double& fx, const double& fy, const double& fz)=0;
  virtual double Bz(const double& fx, const double& fy, const double& fz)=0;
  virtual void Bvec(const double& fx, const double& fy, const double& fz,
                    double *bvec)=0;
  
  virtual std::string GetType(void) const =0;
  virtual std::string GetDescription(void) const =0;
  
  // Galactic center
  gvector gc;
  
};

#endif
