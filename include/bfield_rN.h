//////////////////////////////////////////////////////////
// Filename: bfield_rN.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2011 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the radial power law magnetic field
//
//////////////////////////////////////////////////////////

#ifndef _rNBfield_H
#define _rNBfield_H

#include "globals.h"
#include "bfield.h"
#include <string>

#include <iostream>
#include <iomanip>
#include <fstream>
class R2BField: public BFIELD
{
  
protected:
  
  double index;
  double innerEdge; //[kpc]
  double pitchAngle;  // [rad]
  double spitchAngle;
  double cpitchAngle;
  int    model;         //0=stanev, 1=hmr
  double scaleHeight;
  
public:
  
  R2BField():BFIELD()
  {
  }
  
  virtual ~R2BField()
  {
  }
  
  void SetIndex(const double& s)
  {
    index = s;
  }
  
  double GetIndex(void) const
  {
    return index;
  }
  
  void SetInnerEdge(const double& s)
  {
    innerEdge = s;
  }
  
  double GetInnerEdge(void) const
  {
    return innerEdge;
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
  
  void SetModel(const int& m)
  {
    model=m;
  }
  
  std::string GetModel(void) const 
  {
    if (model==1)
      return "expo";
    else if (model==2)
      return "smooth";
    else
      return "none";
  }
  
  void SetScaleHeight(const double& z)
  {
    scaleHeight = z;
  }
  
  double GetScaleHeight(void) const
  {
    return scaleHeight;
  }
  
  
  
  double zAtten(const double& Z)
  { 
    if(model==1)
      return ( exp(-fabs(Z)/scaleHeight) );
    else if(model==2)
      return ( 1./(1.+(fabs(Z)/scaleHeight))  );
    else
      return 1.;
  }
  
  void Bvec(const double& fx, const double& fy, const double& fz, double *bvec){
    bvec[2]=0.;
    tx = fx;
    ty = fy;
    tz = fz;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      tr = fieldR(tx,ty,tz);      
      if(tr < innerEdge){
        bvec[0]=0;
        bvec[1]=0;
      } else {
        tmag = 1.;
        tphi = BFIELD::fieldPhi(tx,ty);
        tmag *= (norm * zAtten(tz) * pow((tr/solarDistance),index) );
        tcphi = cos(tphi);
        tsphi = sin(tphi);
        bvec[0] = -1.*tmag*(spitchAngle*tcphi - cpitchAngle*tsphi);
        bvec[1] = -1.*tmag*(-1.*spitchAngle*tsphi - cpitchAngle*tcphi);
      }
    }  else {
      bvec[0]=0.;
      bvec[1]=0.;
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
    return 0.;
  }
  
  std::string GetType(void) const 
  {
    return "r^N";
  }
  
  std::string GetDescription(void) const {
    std::ostringstream ff;
    ff << "F " <<  GetType() << ' ' 
    << norm << ' '
    << index << ' '
    << innerEdge << ' '
    << pitchAngle*RAD2DEG << ' '
    << GetModel() << ' '
    << scaleHeight;
    return ff.str();
  }
  
};


#endif
