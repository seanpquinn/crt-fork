//////////////////////////////////////////////////////////
// Filename: bfield_dipole.h
// Authors:  Brian Baughman
//           Michael Sutherland
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the dipole magnetic field
//
//////////////////////////////////////////////////////////

#ifndef _DipoleBField_H
#define _DipoleBField_H
#include "bfield.h"
#include <string.h>

class DipoleBField: public BFIELD
{

 public:

  DipoleBField():BFIELD()
  {
  }

  virtual ~DipoleBField()
  {
  }

  double Bx(const double& fx, const double& fy, const double& fz)
  {
    tx = fx;
    ty = fy;
    tz = fz;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      tr = fieldR(tx,ty,tz);
      tprefactor = norm/(tr*tr*tr);
      tphi = fieldPhi(tx,ty);
      ttheta = fieldTheta(tx,ty,tz);
      txcomp = 3.*sin(ttheta)*cos(ttheta)*cos(tphi);
      return tprefactor*txcomp;
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
      tr = fieldR(tx,ty,tz);
      tprefactor = norm/(tr*tr*tr);
      tphi = fieldPhi(tx,ty);
      ttheta = fieldTheta(tx,ty,tz);
      tycomp = 3.*sin(ttheta)*cos(ttheta)*sin(tphi);
      return -1.*tprefactor*tycomp;
    }
    else{return 0;}
  }
  
  double Bz(const double& fx, const double& fy, const double& fz)
  {
    tx = fx;
    ty = fy;
    tz = fz;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      tr = fieldR(tx,ty,tz);
      tprefactor = norm/(tr*tr*tr);
      ttheta = fieldTheta(tx,ty,tz);
      tctheta = cos(ttheta);
      tzcomp = 1. - 3.*tctheta*tctheta;
      return -1.*tprefactor*tzcomp;
    }
    else{return 0;}
  }
  void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
    tx = fx;
    ty = fy;
    tz = fz;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      tr = fieldR(tx,ty,tz);
      tprefactor = norm/(tr*tr*tr);
      tphi = fieldPhi(tx,ty);
      ttheta = fieldTheta(tx,ty,tz);
      tstheta = sin(ttheta);
      tctheta = cos(ttheta);
      bvec[0] = 3.*tprefactor*tstheta*tctheta*cos(tphi);
      bvec[1] = -3.*tprefactor*tstheta*tctheta*sin(tphi);
      bvec[2] = -1.*tprefactor*(1. - 3.*tctheta*tctheta);
    } else {
      bvec[0]=0;
      bvec[1]=0;
      bvec[2]=0;
    }
  }  
  std::string GetType(void) const 
  {
    return "dipole";
  }
  
  std::string GetDescription(void) const {
    std::ostringstream ff;
    ff << "F " <<  GetType() << ' ' 
    << norm/baseDipoleNorm;
    return ff.str();
  }

};


#endif
