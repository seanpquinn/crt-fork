//////////////////////////////////////////////////////////
// Filename: bfield_uniform.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the uniform magnetic field
//
//////////////////////////////////////////////////////////

#ifndef _UniformBField_H
#define _UnformBField_H

#include "bfield.h"
#include <string.h>

class UniformBField: public BFIELD
{

 private:
  double magBx;
  double magBy;
  double magBz;


 public:
  UniformBField():BFIELD(),magBx(0.),magBy(0.),magBz(0.) {};
  virtual ~UniformBField() {};

  void SetMagBx(const double& bx)
  {
    magBx = bx;
  }
  
  double GetMagBx(void) const
  {
    return magBx;
  }
  
  void SetMagBy(const double& by)
  {
    magBy = by;
  }
  
  double GetMagBy(void) const
  {
    return magBy;
  }
  
  void SetMagBz(const double& bz)
  {
    magBz = bz;
  }
  
  double GetMagBz(void) const
  {
    return magBz;
  }
  
  double Bx(const double& fx, const double& fy, const double& fz)
  {
    return magBx;
  }
  
  double By(const double& fx, const double& fy, const double& fz)
  {
    return magBy;
  }
  
  double Bz(const double& fx, const double& fy, const double& fz)
  {
    return magBz;
  }
  void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
    bvec[0]=magBx;
    bvec[1]=magBy;
    bvec[2]=magBz;
  }  
  std::string GetType(void) const
  {
    return "uniform";
  }
  
  
  std::string GetDescription(void) const {
    std::ostringstream ff;
    ff << "F " <<  GetType() << ' '
    << magBx << ' ' 
    << magBy << ' ' 
    << magBz;
    return ff.str();
  }
};

#endif
