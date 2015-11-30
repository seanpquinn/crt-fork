//////////////////////////////////////////////////////////
// Filename: bfield_pshirkov_halo.h
// Authors:  Azadeh Keivani
//
// Copyright 2012 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the halo component of the magnetic field
//					based on Pshirkov's model
//
//////////////////////////////////////////////////////////

#ifndef _PSHIRKOV_HALO_H
#define _PSHIRKOV_HALO_H

#include "globals.h"
#include "bfield.h"
#include <string>

#include <iostream>
#include <iomanip>
#include <fstream>

class PSHIRKOVHALO: public BFIELD
{
  
protected:
  
  double BH0;               // this grouping is
  double zH0, zH1, zH1a;     //  for Halo component
  double RH0;
  
  bool   model;
  double haloprenorm;
  
public:
  PSHIRKOVHALO():BFIELD(),
  BH0(0.),zH0(0.),zH1(0.),zH1a(0.),RH0(0.),
  model(1), haloprenorm(0.)
  {
    
  };
	
  virtual ~PSHIRKOVHALO(){
    
  };    
  
  void SetModel(const bool& m)
  {
    model=m;
  }
  
  std::string GetType(void) const
  {
    if(model)
      return "sq";
    else
      return "abs";
  }
  

  void SetHaloNorm(const double& n)
  {
    // set halo norm only
    BH0 = n;
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
    ff << "F pshirkov-halo ";
    std::string tp = GetType();
      ff << tp << ' '
	  << BH0 << ' '
      << zH0 << ' '
      << zH1 << ' '
      << zH1a << ' '
      << RH0;
      return ff.str();
  }
  
  
  /*double D1(const double& r, const double& z) {
    retval = 1.;
    if(r > Rc)
      retval *= norm * exp(-1.*((r-solarDistance)/R0) - (fabs(z)/z0));
    else
      retval *= Bc;
    
    return retval;
  }*/
  
  /*double D2(const double& r, const double& phi) {
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
  }*/
  
  double halomag(const double& r, const double& z) {
    retval = 0;
    if( fabs(z) < zH0 )
      tA = zH1;
    else
      tA = zH1a;
    tcap = (fabs(z) - zH0) / tA;
    tdenom = 1. + (tcap*tcap);
	  if(model){
    retval = (BH0 * (1./tdenom) * exp(-1.*((r-RH0)/RH0)*((r-RH0)/RH0)));
    if(z < 0.)
      retval *= -1.;
    return retval;
	  }
	  else {
		  retval = (BH0 * (1./tdenom) * exp(-1.*fabs(((r-RH0)/RH0))));
		  if(z < 0.)
			  retval *= -1.;
		  return retval;
	  }
  }
  
  // Returns a vector of B magnitudes faster
  void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
    tx = fx;
    ty = fy;
    tz = fz;
    bvec[2] = 0.;
    BFIELD::gal2gc(tx,ty,tz);
    if(withinField(tx,ty,tz)){
      tHmag = 1.;
      tr = BFIELD::fieldRho(tx,ty);
      tphi = BFIELD::fieldPhi(tx,ty);
      tcphi = cos(tphi);
      tsphi = sin(tphi);
  
      tHmag *= halomag(tr,tz);

      bvec[0] = tHmag*tsphi;
      bvec[1] = tHmag*tcphi;
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
