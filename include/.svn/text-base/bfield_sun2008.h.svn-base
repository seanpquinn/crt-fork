//////////////////////////////////////////////////////////
// Filename: bfield_sun2008.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2011 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the Sun2008 magnetic field
//
//              If the field vector points inward toward
//                galactic center (clockwise tangent to
//                the spiral viewed from North Galactic Pole
//                is inside the circle of radius R), pitch
//                angle is NEGATIVE
//
//////////////////////////////////////////////////////////

#ifndef _SUN2008_H
#define _SUN2008_H

#include "globals.h"
#include "bfield.h"
#include <string>

#include <iostream>
#include <iomanip>
#include <fstream>
class SUN2008: public BFIELD
{
  
protected:
  double pitchAngle;  // this grouping is for
  double spitchAngle; //   ASS+RING
  double cpitchAngle;
  double Bc;
  double Rc;
  // norm == B0   in the paper's notation
  double R0;
  double z0;
  double R1, R2, R3;
  
  double pitchAngle1;  // this grouping is
  double spitchAngle1; //  for BSS
  double cpitchAngle1;
  double pitchAngle2;
  double spitchAngle2;
  double cpitchAngle2;   
  double Rb, Rb1, Rb2;
  
  double BH0;               // this grouping is
  double zH0, zH1, zH1a;     //  for Halo component
  double RH0;
  
  bool   model;
  double haloprenorm;
  
public:
  SUN2008():BFIELD(),
  pitchAngle(0.),spitchAngle(0.),cpitchAngle(1.),
  Bc(0.),Rc(0.),R0(0.),z0(0.),R1(0.),R2(0.),R3(0.),
  pitchAngle1(0.),spitchAngle1(0.),cpitchAngle1(1.),
  pitchAngle2(0.),spitchAngle2(0.),cpitchAngle2(1.),
  Rb(0.),Rb1(0.),Rb2(0.),
  BH0(0.),zH0(0.),zH1(0.),RH0(0.),
  model(1), haloprenorm(0.)
  {
    
  };
  virtual ~SUN2008(){
    
  };    
  
  void SetModel(const bool& m)
  {
    model=m;
  }
  
  std::string GetType(void) const
  {
    if(model)
      return "ass+ring";
    else
      return "bss";
  }
  
  void SetNorm(const double& n)
  {
    // set DISK norm for (r>Rc) ONLY;
    // for (r < Rc) the field is whatever Bc is set
    norm = n;
  }
  void SetHaloNorm(const double& n)
  {
    // set halo norm only
    BH0 = n;
  }
  
  void SetPitchAngle(const double& p)
  {
    pitchAngle = p;
    spitchAngle = sin(p);
    cpitchAngle = cos(p);
    
  }
  double GetPitchAngle(void) const
  {
    return pitchAngle;
  }
  void SetBc(const double& a)
  {
    Bc = a;
  }
  void SetRc(const double& a)
  {
    Rc = a;
  }
  void SetR0(const double& a)
  {
    R0 = a;
  }
  void Setz0(const double& a)
  {
    z0 = a;
  }
  void SetR1(const double& a)
  {
    R1 = a;
  }
  void SetR2(const double& a)
  {
    R2 = a;
  }
  void SetR3(const double& a)
  {
    R3 = a;
  }
  double GetBc(void) const
  {
    return Bc;
  }
  double GetRc(void) const
  {
    return Rc;
  }
  double GetR0(void) const
  {
    return R0;
  }
  double Getz0(void) const
  {
    return z0;
  }
  double GetR1(void) const
  {
    return R1;
  }
  double GetR2(void) const
  {
    return R2;
  }
  double GetR3(void) const
  {
    return R3;
  }
  
  
  void SetPitchAngle1(const double& p)
  {
    pitchAngle1 = p;
    spitchAngle1 = sin(p);
    cpitchAngle1 = cos(p);
    
  }
  double GetPitchAngle1(void) const
  {
    return pitchAngle1;
  }
  void SetPitchAngle2(const double& p)
  {
    pitchAngle2 = p;
    spitchAngle2 = sin(p);
    cpitchAngle2 = cos(p);
    
  }
  double GetPitchAngle2(void) const
  {
    return pitchAngle2;
  }
  void SetRb(const double& a)
  {
    Rb = a;
  }
  double GetRb(void) const
  {
    return Rb;
  }
  void SetRb1(const double& a)
  {
    Rb1 = a;
  }
  double GetRb1(void) const
  {
    return Rb1;
  }
  void SetRb2(const double& a)
  {
    Rb2 = a;
  }
  double GetRb2(void) const
  {
    return Rb2;
  }
  
  
  void SetzH0(const double& a)
  {
    zH0 = a;
  }
  double GetzH0(void) const
  {
    return zH0;
  }
  void SetzH1(const double& a)
  {
    zH1 = a;
  }
  double GetzH1(void) const
  {
    return zH1;
  }
  void SetzH1a(const double& a)
  {
    zH1a = a;
  }
  double GetzH1a(void) const
  {
    return zH1a;
  }
  void SetRH0(const double& a)
  {
    RH0 = a;
  }
  double GetRH0(void) const
  {
    return RH0;
  }
  
  std::string GetDescription(void) const {
    std::ostringstream ff;
    ff << "F sun2008 ";
    std::string tp = GetType();
    if (tp=="ass+ring") {
      ff << tp << ' '
      << pitchAngle*RAD2DEG << ' '
      << norm << ' '
      << R0 << ' '
      << Bc << ' '
      << Rc << ' '
      << z0 << ' '
      << R1 << ' '
      << R2 << ' '
      << R3 << ' '
      << BH0 << ' '
      << zH0 << ' '
      << zH1 << ' '
      << zH1a << ' '
      << RH0;
    } else if (tp=="bss") {
      ff << tp << ' '
      << norm << ' '
      << R0 << ' '
      << Bc << ' '
      << Rc << ' '
      << z0 << ' '
      << Rb << ' '
      << pitchAngle1*RAD2DEG << ' '
      << Rb1 << ' '
      << pitchAngle2*RAD2DEG << ' '
      << Rb2 << ' '
      << BH0 << ' '
      << zH0 << ' '
      << zH1 << ' '
      << zH1a << ' '
      << RH0;
    }
    return ff.str();
  }
  
  
  double D1(const double& r, const double& z) {
    retval = 1.;
    if(r > Rc)
      retval *= norm * exp(-1.*((r-solarDistance)/R0) - (fabs(z)/z0));
    else
      retval *= Bc;
    
    return retval;
  }
  
  double D2(const double& r, const double& phi) {
    retval = -1.;
    if(model){
      if(r > R1)
        return (retval);
      else if( (R2 < r) && (r <= R1) )
        return (-1.*retval);
      else if( (R3 < r) && (r <= R2) )
        return (retval);
      else  // r <= R3
        return (-1.*retval);
    } else {
      retval *= -1.;
      ttempscale = 1.;
      tbeta = 1.;
      if(r > Rb){
        ttempscale *= Rb1;
        tbeta *= cpitchAngle1 / spitchAngle1;
      } else {
        ttempscale *= Rb2;
        tbeta *= cpitchAngle2 / spitchAngle2;
      }
      retval *= -1. * cos(phi + tbeta*log(r/ttempscale) );
      return retval;
    }
  }
  
  double halomag(const double& r, const double& z) {
    retval = 0;
    if( fabs(z) < zH0 )
      tA = zH1;
    else
      tA = zH1a;
    tcap = (fabs(z) - zH0) / tA;
    tdenom = 1. + (tcap*tcap);
    retval = -1.*(BH0 * (1./tdenom) * (r/RH0) * exp(-1.*((r-RH0)/RH0)));
    if(z < 0.)
      retval *= -1.;
    return retval;
  }
  
  // Returns a vector of B magnitudes faster
  void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
    tx = fx;
    ty = fy;
    tz = fz;
    bvec[2] = 0.;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      tmag = 1.;
      tHmag = 1.;
      tr = BFIELD::fieldRho(tx,ty);
      tphi = BFIELD::fieldPhi(tx,ty);
      tcphi = cos(tphi);
      tsphi = sin(tphi);
      tsp = 1.;
      tcp = 1.;
      
      tmag *= D1(tr,tz);
      tmag *= D2(tr,tphi);
      tHmag *= halomag(tr,tz);
      if(model){
        tsp *= spitchAngle;
        tcp *= cpitchAngle;
      } else {
        if(tr > Rb1){
          tsp *= spitchAngle1;
          tcp *= cpitchAngle1;
        } else {
          tsp *= spitchAngle2;
          tcp *= cpitchAngle2;
        }
      }
      bvec[0] = tmag*(tsp*tcphi - tcp*tsphi) + (tHmag*tsphi);
      bvec[1] = tmag*(-1.*tsp*tsphi - tcp*tcphi) + (tHmag*tcphi);
    } else {
      bvec[0]=0;
      bvec[1]=0;
    }
  }
  
  
  double Bx(const double& fx, const double& fy, const double& fz)
  {
    Bvec(fx,fy,fz,cvec);
    return cvec[0];
  }
  
  double By(const double& fx, const double& fy, const double& fz)
  {
    Bvec(fx,fy,fz,cvec);
    return cvec[1];
  }
  
  double Bz(const double& fx, const double& fy, const double& fz)
  {
    return 0;
  }
  
  
};


#endif
