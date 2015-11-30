//////////////////////////////////////////////////////////
// Filename: bfield_toroid.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the toroidal magnetic field
//
//////////////////////////////////////////////////////////

#ifndef _ToroidalBField_H
#define _ToroidalBField_H
#include "bfield.h"
#include <string.h>

class ToroidalBField: public BFIELD
{
  
protected:
  double toroidalRadius;
  double scale1;
  double lorentzHalfWidth;
  double Bmax;
  
public:
  
  ToroidalBField():BFIELD()
  {
  }
  
  virtual ~ToroidalBField()
  {
  }
  
  void SetToroidalRadius(const double& a)
  {
    toroidalRadius = a;
  }
  
  double GetToroidalRadius(void) const
  {
    return toroidalRadius;
  }
  
  void SetHeight(const double& a)
  {
    scale1 = a;
  }
  
  double GetHeight(void) const
  {
    return scale1;
  }
  
  void SetHalfWidth(const double& a)
  {
    lorentzHalfWidth = a;
  }
  
  double GetHalfWidth(void) const
  {
    return lorentzHalfWidth;
  }
  
  void SetBmax(const double& a)
  {
    Bmax = a;
  }
  
  double GetBmax(void) const
  {
    return Bmax;
  }
  
  double lorentzFactor(const double& z)
  {
    return ( ((fabs(z)-scale1)/lorentzHalfWidth)*((fabs(z)-scale1)/lorentzHalfWidth) );
  }
  
  double Bmag(const double& x, const double& y, const double& z)
  {
    if( fieldRho(x,y) < toroidalRadius )
    {
      return (Bmax / (1.+lorentzFactor(z)));
    }
    else
    {
      return (Bmax * exp(-fieldRho(x,y)/toroidalRadius) / (1.+lorentzFactor(z)));
    }
  }
  
  double Bx(const double& fx, const double& fy, const double& fz)
  {
    tx = fx;
    ty = fy;
    tz = fz;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      return (-Bmag(tx,ty,tz)*sin(-1.*fieldPhi(tx,ty)));
    }
    else{return 0;}
  }
  double By(const double& fx, const double& fy, const double& fz)
  {
    tx = fx;
    ty = fy;
    tz = fz;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      return (Bmag(tx,ty,tz)*cos(-1.*fieldPhi(tx,ty)));
    }
    else{return 0;}
  }
  
  double Bz(const double& fx, const double& fy, const double& fz)
  {
    return 0;
  }
  void Bvec(const double& fx, const double& fy, const double& fz, double *bvec){
    tx = fx;
    ty = fy;
    tz = fz;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      tmag = Bmag(tx,ty,tz);
      tphi = fieldPhi(tx,ty);
      bvec[0] = (-tmag*sin(-1.*tphi));
      bvec[1] = (tmag*cos(-1.*tphi));
      bvec[2] = 0;
    } else {
      bvec[0]=0;
      bvec[1]=0;
      bvec[2]=0;
    }
  }  
  std::string GetType(void) const
  {
    return "toroidal";
  }
  
  std::string GetDescription(void) const {
    std::ostringstream ff;
    ff << "F " <<  GetType() << ' ' 
    << GetToroidalRadius() << ' '
    << GetHeight() << ' '
    << GetHalfWidth() << ' ' 
    << GetBmax();
    return ff.str();
  }
  
};


#endif
