//////////////////////////////////////////////////////////
// Filename: bfield_ssfield.h
// Authors:  Brian Baughman
//           Michael Sutherland 
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

#ifndef _SSFIELD_H
#define _SSFIELD_H

#include "globals.h"
#include "bfield.h"
#include <string>

#include <iostream>
#include <iomanip>
#include <fstream>
class SSFIELD: public BFIELD
  {
    
  protected:
    double pitchAngle;  // [rad]
    double spitchAngle;
    double cpitchAngle;
    double scale1;      // galactocentric distance along solar axis
    //  at which field is maximum strength
    double scale2;      // tanh cutoff
    double scale3;      // first z attenuation length
    double scale4;      // second z attenuation length
    double zcut;        // height in z where z-atten length goes
    //  from scale3 -> scale4
    int zparity;
    bool axisymetric;
    int    model;         //0=stanev, 1=hmr
    double prenorm;
    
  public:
    SSFIELD(int fzparity, bool faxisymetric):BFIELD(),
    pitchAngle(0.),spitchAngle(0.),cpitchAngle(1.),
    scale1(0.),scale2(0.),scale3(0.),scale4(0.),
    zcut(0.),zparity(fzparity),axisymetric(faxisymetric),
    model(1), prenorm(0.)
    {

    };
    virtual ~SSFIELD(){
      
    };    
    
    void SetModel(const bool& m)
    {
      model=m;
    }
    
    std::string GetModel(void) const
    {
      if(model)
        return "hmr";
      else
        return "stanev";
    }
    
    void SetNorm(const double& n)
    {
      // calculate raw strength at solar position, then rescale
      prenorm = fabs(Brmag(solarDistance, 0)*prefactor(solarDistance, 0));
      norm = n/prenorm;
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
    
    void SetScale1(const double& s1)
    {
      scale1 = s1;
    }
    
    double GetScale1(void) const
    {
      return scale1;
    }
    
    void SetScale2(const double& s2)
    {
      scale2 = s2;
    }
    
    double GetScale2(void) const
    {
      return scale2;
    }
    
    void SetScale3(const double& s3)
    {
      scale3 = s3;
    }
    
    double GetScale3(void) const
    {
      return scale3;
    }
    
    void SetScale4(const double& s4)
    {
      scale4 = s4;
    }
    
    double GetScale4(void) const
    {
      return scale4;
    }
    
    void SetZcut(const double& zc){
      zcut = zc;
    }
    
    double GetZcut(void) const {
      return zcut;
    }
    
    double Brmag(const double& r, const double& z)
    {
      retval = 1.;
      if(z<0)
        retval *= static_cast<double>(zparity);
      if(!model){ //stanev
        if(r<4.)
          retval *= -1.*(3.*solarDistance/4.);
        else
          retval *= -1.*(3.*solarDistance/r);
        return retval;
      }
      else{ //hmr
        retval *= -1.*(3.*solarDistance/r)*pow(tanh(r/scale2), 3.);
        return retval;
      }
      
    }
    double xfrac(const double& phi)
    {
      return spitchAngle*cos(phi) - cpitchAngle*sin(phi);
    }
    double yfrac(const double& phi)
    {
      return -1.*spitchAngle*sin(phi) - cpitchAngle*cos(phi);
    }
    
    double prefactor(const double& r,const double& phi)
    {
      tprefactor=1.;
      if(model==0){ //stanev
        if(r<4.)
          tprefactor *= cos(phi - (cpitchAngle/spitchAngle)*log(4./scale1) );
        else
          tprefactor *= cos(phi - (cpitchAngle/spitchAngle)*log(r/scale1) );
        if (axisymetric==true)
          tprefactor= fabs(tprefactor);
      }
      
      else if( model==1){ //hmr
        tprefactor *= cos(phi - (cpitchAngle/spitchAngle)*log(r/scale1) );
        if(axisymetric==true)
          tprefactor *= tprefactor;
      } 
      return tprefactor;
    }    
    
    double zAtten(const double& Z)
    { 
      if(!model){ //stanev
        retval = fabs(Z);
        if( retval <= zcut)
          return ( exp(-1.*retval/scale3) );
        else{
          return ( exp(-1.*retval/scale4) );
        }	
      }
      else{ //hmr
        retval = 1.;
        if(zparity == -1)
          retval *= fabs(tanh(Z/0.02)); // This line should not have 0.02 hard coded in.
        retval *= 0.5*( (1./cosh(Z/scale3) ) + (1./cosh(Z/scale4) ) );
        return retval;
      }
    }
    // Returns a vector of B magnitudes faster
    void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
      tx = fx;
      ty = fy;
      tz = fz;
      BFIELD::gal2gc(tx,ty,tz);
      if(withinField(tx,ty,tz)){
        tr = BFIELD::fieldRho(tx,ty);
        tphi = BFIELD::fieldPhi(tx,ty);
        tmag = norm*prefactor(tr,tphi)*Brmag(tr,tz)*zAtten(tz);
        tcphi = tx/tr;//cos(phi);
        tsphi = ty/tr;//sin(phi);
        bvec[0] = tmag*(spitchAngle*tcphi - cpitchAngle*tsphi);
        bvec[1] = tmag*(-1.*spitchAngle*tsphi - cpitchAngle*tcphi);
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
        tphi = BFIELD::fieldPhi(tx,ty);
        txcomp = norm*prefactor(tr,tphi)*Brmag(tr,tz)*xfrac(tphi)*zAtten(tz);
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
        tphi = BFIELD::fieldPhi(tx,ty);
        tycomp = norm*prefactor(tr,tphi)*Brmag(tr,tz)*yfrac(tphi)*zAtten(tz);
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
      ff << "F " <<  GetType() << ' '
      << GetModel() << ' '
      << norm*prenorm << ' '
      << pitchAngle/DEG2RAD << ' ' 
      << scale1 << ' ' 
      << scale2 << ' ' 
      << scale3 << ' '
      << scale4 << ' ' 
      << zcut;
      return ff.str();
    }
  };


#endif
