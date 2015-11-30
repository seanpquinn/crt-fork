//////////////////////////////////////////////////////////
// Filename: bfield_ssfield.h
// Original Authors:  Brian Baughman
//           Michael Sutherland 
// Added Pshirkov's Model: Azadeh Keivani (Feb 2012)
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the spiral magnetic field models
//
//              If the field vector points inward toward
//                galactic center (clockwise tangent to
//                the spiral viewed from North Galactic Pole
//                is inside the circle of radius R), pitch
//                angle is NEGATIVE
//
//////////////////////////////////////////////////////////

#ifndef _PSHIRKOV_DISK_H
#define _PSHIRKOV_DISK_H

#include "globals.h"
#include "bfield.h"
#include <string>

#include <iostream>
#include <iomanip>
#include <fstream>

class PSHIRKOVDISK: public BFIELD
  {
    
  protected:
    double pitchAngle;  // [rad]
    double spitchAngle;
    double cpitchAngle;
    double Rc;      // the radius of the central region where the disk field is assumed to be constant
    double d; //the distance to the first field reversal
    double z0;      // z attenuation length
 //   bool zparity; //?
    bool axisymetric;
    bool model;
    double prenorm;
    
  public:
    PSHIRKOVDISK():BFIELD(),
    pitchAngle(0.),spitchAngle(0.),cpitchAngle(1.),
    Rc(0.),d(0.),z0(0.),
	model(1),prenorm(0.)
    {

    };
    virtual ~PSHIRKOVDISK(){
      
    };    
    
	  void SetModel(const bool& mo)
	  {
		  model=mo;
	  }
	  
	  std::string GetType(void) const
	  {
	  if(model)
		  return "bss_s";
	  else 
		  return "ass_s";
	  }
    
    void SetNorm(const double& n)
    {
      norm = n;
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
    
    void SetRc(const double& s1)
    {
      Rc = s1;
    }
    
    double GetRc(void) const
    {
      return Rc;
    }
    
	  void Setd(const double& s2)
	  {
		  d = s2;
	  }
	  
	  double Getd(void) const
	  {
		  return d;
	  }
    
    void Setz0(const double& s3)
    {
      z0 = s3;
    }
    
    double Getz0(void) const
    {
      return z0;
    }
    
    double Brmag(const double& r, const double& z)
    {
     retval = 1.;
        
		if(r>Rc)
          retval *= solarDistance/(r*cos((cpitchAngle/spitchAngle)*log(1+d/solarDistance)-pi/2)); 
        else
          retval *= solarDistance/(Rc*cos((cpitchAngle/spitchAngle)*log(1+d/solarDistance)-pi/2));
		
      return retval;
    }
	  
    double xfrac(const double& theta)
    {
      return spitchAngle*cos(theta) - cpitchAngle*sin(theta);
    }
    double yfrac(const double& theta)
    {
      return -1.*spitchAngle*sin(theta) - cpitchAngle*cos(theta);
    }
	  
    double prefactor(const double& r,const double& theta)
    {
      tprefactor=1.;
		
          tprefactor *= cos(theta - (cpitchAngle/spitchAngle)*(log(r/solarDistance)-log(1+d/solarDistance))-pi/2 );
		
        if (!model)
		tprefactor= fabs(tprefactor);

      return tprefactor;
    }    
    
    double zAtten(const double& Z)
    { 
        retval = fabs(Z);
		return ( exp(-1.*retval/z0) );
    }
    // Returns a vector of B magnitudes faster
    void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
      tx = fx;
      ty = fy;
      tz = fz;
      BFIELD::gal2gc(tx,ty,tz);
      if(withinField(tx,ty,tz)){
        tr = BFIELD::fieldRho(tx,ty);
        ttheta = BFIELD::fieldPhi(tx,ty);
        tmag = norm*prefactor(tr,ttheta)*Brmag(tr,tz)*zAtten(tz);
        tctheta = tx/tr;//cos(theta);
        tstheta = ty/tr;//sin(theta);
        bvec[0] = tmag*(spitchAngle*tctheta - cpitchAngle*tstheta);
        bvec[1] = tmag*(-1.*spitchAngle*tstheta - cpitchAngle*tctheta);
        bvec[2] = 0.;
      }  else {
        bvec[0]=0;
        bvec[1]=0;
        bvec[2]=0;
      }
    }
    
    
    double Bx(const double& fx, const double& fy, const double& fz)
    {
      tx = fx;
      ty = fy;
      tz = fz;
      BFIELD::gal2gc(tx,ty,tz);
      if(withinField(tx,ty,tz)){
        tr = BFIELD::fieldRho(tx,ty);
        ttheta = BFIELD::fieldPhi(tx,ty);
        txcomp = norm*prefactor(tr,ttheta)*Brmag(tr,tz)*xfrac(ttheta)*zAtten(tz);
        return txcomp;
      } else
        return 0.;
    }
    
    double By(const double& fx, const double& fy, const double& fz)
    {
      tx = fx;
      ty = fy;
      tz = fz;
      BFIELD::gal2gc(tx,ty,tz);
      if(withinField(tx,ty,tz)){
        tr = BFIELD::fieldRho(tx,ty);
        ttheta = BFIELD::fieldPhi(tx,ty);
        tycomp = norm*prefactor(tr,ttheta)*Brmag(tr,tz)*yfrac(ttheta)*zAtten(tz);
        return tycomp;
      } else
        return 0.;
    }
    
    double Bz(const double& fx, const double& fy, const double& fz)
    {
      return 0.;
    }
    
    std::string GetDescription(void) const  {
      std::ostringstream ff;
      ff << "F pshirkov-disk "
         << GetType() << ' '
         << norm << ' '
         << pitchAngle/DEG2RAD << ' ' 
         << Rc << ' ' 
         << d << ' ' 
	     << z0;
      return ff.str();
    }
  };


#endif
