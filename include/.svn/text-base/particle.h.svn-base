//////////////////////////////////////////////////////////
// Filename: particle.h
// Authors:  Michael Sutherland
//           Brian Baughman
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the particle class
//
//////////////////////////////////////////////////////////

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "globals.h"
#include "gvector.h"

class Particle
{

 private:
  double pMass;            //[kg]
  double pCharge;          //[C]
  double pEnergy;          //[eV]
  double pGamma;
  double pPrefactor;
  double pVelMag;          //[kpc/yr]
  gvector pP0;             //[kpc]
  gvector pV0;             //[kpc/yr]
  gvector pP;             //[kpc]
  gvector pV;             //[kpc/yr]  
  bool hitSt;
  double DOCA;              //[kpc]

 public:
  Particle():pMass(0.),pCharge(0.),
  pEnergy(0.),pGamma(0.),pPrefactor(1.),pVelMag(0.),
  pP0(gvector(0.,0.,0.)),pV0(gvector(0.,0.,0.)),
  pP(gvector(0.,0.,0.)),pV(gvector(0.,0.,0.)),hitSt(false),DOCA(0.)

    {};
  virtual ~Particle() {
  
  };

  void SetHit(bool h)
  {
    hitSt = h;
  }
  
  bool GetHit(void) const
  {
    return hitSt;
  }
  
  void SetMass(double m){
    pMass=m;
  }
  
  double GetMass(void)  const
  {  
    return pMass;
  }
  
  void SetCharge(double q)
  {
    pCharge = q;
  }
  
  double GetCharge(void)  const
  {  
    return pCharge;
  }
  
  void SetEnergy(double en)
  { 
    pEnergy = en;
    pGamma = pEnergy/pMass;
    pVelMag = SoLgal * sqrt(1.-(1./(pGamma*pGamma) ) );
    pPrefactor = crt2gaussian*pCharge/pEnergy;
  }
  
  double GetEnergy(void) const{
    return pEnergy;
  }
  
  
  double GetGamma(void) const
  {
    return pGamma;
  }
  
  double GetPrefactor(void) const
  {
    return pPrefactor;
  }
  
  double GetVelMag(void) const
  {
    return pVelMag;
  }
  gvector GetP0(void) const
  {
    return pP0;
  }
  gvector GetP(void) const
  {
    return pP;
  }
  
  gvector GetV0(void) const
  {
    return pV0;
  }
  gvector GetV(void) const
  {
    return pV;
  }
  
  void SetP0(double fx,double fy, double fz) {
    pP0.setcoords(fx,fy,fz);
    pP.setcoords(fx,fy,fz);
  }
  void SetP(double fx,double fy, double fz) {
    pP.setcoords(fx,fy,fz);
  }
  void SetV0(double fx,double fy, double fz) {
    pV0.setcoords(fx,fy,fz);
    pV.setcoords(fx,fy,fz);
  }
  void SetV(double fx,double fy, double fz) {
    pV.setcoords(fx,fy,fz);
  }
  void SetP0(gvector& fgv) {
    pP0=fgv;
    pP=fgv;
  }
  void SetP(gvector& fgv) {
    pP=fgv;
  }
  void SetV0(gvector& fgv) {
    pV0=fgv;
    pV=fgv;
  }
  void SetV(gvector& fgv) {
    pV=fgv;
  }
  void SetDOCA(double doca)
  {
    DOCA=doca;
  }
  
  double GetDOCA(void) const
  {
    return DOCA;
  }
  
  double GetPosMag(void) const
  {
    return pP.r();
  }
  
  //double GetPdotV(void) const
  //{
  // return pP.dot(pV);
  //}

  double GetTimetoDOCA(gvector fgv)
  {
    // This algorithm assumes that pV is constant over the next time step
    // That is pP(t-now)=pP(now)+pV*(t-now)
    double time = -1./(pV.dot(pV));
    time *= pV.dot( (pP - fgv) );
    return time;
  }
};


#endif
