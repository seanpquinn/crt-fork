//////////////////////////////////////////////////////////
// Filename: detector.cpp
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the detector
//
//////////////////////////////////////////////////////////

#include "detector.h"
#include "transforms.h"
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h> 
#include <gsl/gsl_randist.h>

Detector* Detector::pinstance = 0;
Detector* Detector::Instance () 
{
  if (pinstance == 0)  // is it the first call?
  {  
    pinstance = new Detector; // create sole instance
  }
  return pinstance; // address of sole instance
}


void Detector::SetP0(double fx,double fy, double fz) {
  pD0.setcoords(fx,fy,fz);
}

void Detector::SetP0(gvector& fgv) {
  pD0=fgv;
}

gvector Detector::GetP0(void)
{
  return pD0;
}

double Detector::GetPosMag(void)
{
  return pD0.r();
}

void Detector::SetDetR(double d){
  detR = d;
}

double Detector::GetDetR(void){
  return detR;
}

void Detector::SetTmax(double a){
  tmax = a;
}

double Detector::GetTmax(void){
  return tmax;
}

void Detector::SetTmin(double b){
  tmin = b;
}

double Detector::GetTmin(void){
  return tmin;
}

void Detector::SetTstep(double s){
  tstep = s;
}

double Detector::GetTstep(void){
  return tstep;
}
void Detector::SetDetDec(double d){
  detDec = d;
}

double Detector::GetDetDec(void){
  return detDec;
}
void Detector::SetDetFOV(double d){
  detFOV = d;
}

double Detector::GetDetFOV(void){
  return detFOV;
}

std::string Detector::GetType(void){
  return dtype;
}

void Detector::SetType(std::string s){
  dtype = s;
}
void Detector::SetRndPtr(void * rng) {
  rng_ptr = rng;
};

bool Detector::GetExposure(void){
  return exposure;
}

gvector Detector::getshift(gvector &lx, gvector &ly,
                           long srcnum, double frac, double angle) {
  if (dtype=="area") {
    std::map<long,long>::iterator srciter = srccnts.find(srcnum);
    long nevt=0;
    if(srciter != srccnts.end() ) {
      srciter->second++;
      nevt=srciter->second%npnts;
    } else {
      srccnts.insert(std::pair<long,long>(srcnum,0));
    }
    angle=0;
    if (nevt==0) {
      frac=0.;
    }
    else {
      frac=1.;
      angle=2.*pi*float(nevt-1)/float(npnts-1);
    }
  }
  return (lx*sin(angle)+ly*cos(angle))*detR*frac;
}


double Detector::zeta(double dec) {
  double z = cos(detFOV)-sin(detDec)*sin(dec);
  z /= cos(detDec)*cos(dec);
  return z;
}

double Detector::alpham(double z) {
  if (z>1)
    return 0;
  else if (z<-1)
    return pi;
  else
    return acos(z);
}

double Detector::relw(double dec) {
  double z = zeta(dec);
  double am = alpham(z);
  double rw = cos(detDec)*cos(dec)*sin(am);
  rw += am*sin(detDec)*sin(dec);
  return rw;
}

double Detector::exposureValueGal(const double &glong, const double &glat){
  double ra,dec;
  transforms::gal2eq(glong,glat,ra,dec);
  return relw(dec);
}

double Detector::exposureValue(const double &ra, const double &dec){
  return relw(dec);
}

void Detector::SampleExposure(double &glong, double &glat) {
  if (rng_ptr==0)
    return;
  static double tra, tdec,minbnd,maxbnd,randN;
  gsl_rng* cast_rng_ptr = (gsl_rng*) rng_ptr;
  tra = gsl_ran_flat(cast_rng_ptr,-pi,pi);
  minbnd = detDec-detFOV;
  if (minbnd<-pi*.5) 
    minbnd=-pi*.5;
  minbnd = sin(minbnd);
  maxbnd = detDec+detFOV;
  if (maxbnd>pi*.5) 
    maxbnd=pi*.5;
  maxbnd = sin(maxbnd);
  tdec = asin(gsl_ran_flat(cast_rng_ptr,minbnd,maxbnd));
  randN = gsl_ran_flat(cast_rng_ptr,0,detMaxExpo);
  while(randN > relw(tdec)) {
    randN = gsl_ran_flat(cast_rng_ptr,0,detMaxExpo);
    tdec = asin(gsl_ran_flat(cast_rng_ptr,minbnd,maxbnd));
  }
  transforms::eq2gal(tra,tdec,glong,glat);
}

std::string Detector::GetDescription(void) {
  std::ostringstream ff;
  ff << "D " <<  GetP0().x() << ' ' << GetP0().y() << ' ' << GetP0().z() << ' ' << GetDetR();
  return ff.str();
}


