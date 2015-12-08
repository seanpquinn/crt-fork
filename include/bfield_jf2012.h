//////////////////////////////////////////////////////////
// Filename: bfield_jf2012.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2012 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the Jansson/Farrar magnetic field
//
//
//////////////////////////////////////////////////////////

#ifndef _JF2012_H
#define _JF2012_H

#include "globals.h"
#include "bfield.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdexcept>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

class JF2012: public BFIELD
{
  
protected:
  gsl_rng *rng_ptr;      // pointer to rng
  bool bl_use_stri;
  bool bl_generate;      // true == generate cells, otherwise read from file
  std::string rndmfile;

  // disk field
  double bdisk[8];  // field strengths of arms at r=5 kpc
  double fxsec[8];  // relative cross-sectional areas of spirals
                                          //   b8 is determined from other 7
  double rc_B[8];                        // radii where each arm crosses the
                                          //   negative x-axis
  double pitch;
  double bring;         // ring field strength 3<r<5 kpc
  double hdisk, wdisk;  // disk/halo transistion and width

  // toroidal halo
  double Bn, Bs;        // northern, southern halo field strength
  double rn, rs;        // northern, southern transistion radius
  double wh, z0;        // transistion width and vertical scale height

  // X halo
  double BX;            // field strength at origin
  double thetaX0;       // elevation angle at z=0, r > rXc
  double rXc;           // radius where elevation angle = thetaX0
  double rX;            // exponential scale height

  double beta, fsigma;
  double sqrtbeta;

  double tR, tDmag, tXmag;
  double spitch, cpitch;

  static const size_t Nsign = 64000000;
  char mastergrid[Nsign];
  
public:
  // CONSTRUCTOR FOR REGULAR FIELD ONLY (NO STRIATED)
  JF2012(gsl_rng *r, 
    double sigma, 
    bool stri,
    double JF12b1,
    double JF12b2,
    double JF12b3,
    double JF12b4,
    double JF12b5,
    double JF12b6,
    double JF12b7,
    double JF12bring,
    double JF12hdisk,
    double JF12wdisk,
    double JF12Bn,
    double JF12Bs,
    double JF12rn,
    double JF12rs,
    double JF12wh,
    double JF12z0,
    double JF12Bx,
    double JF12ThX,
    double JF12rXc,
    double JF12rX
    ):BFIELD(),
  rng_ptr(r),bl_use_stri(stri),bl_generate(false),rndmfile("")  
  {
  // these parameters are 'constant' within the model
  //  (no sigma values provided)
  pitch = 11.5;
  fsigma = sigma;

  fxsec[0] = 0.130;
  fxsec[1] = 0.165;
  fxsec[2] = 0.094;
  fxsec[3] = 0.122;
  fxsec[4] = 0.13;
  fxsec[5] = 0.118;
  fxsec[6] = 0.084;
  fxsec[7] = 0.156;
  
  // these parameters are subject to sigma shifting
  
  bdisk[0] = JF12b1;  rc_B[0] =  5.1;
  bdisk[1] = JF12b2;  rc_B[1] =  6.3;
  bdisk[2] = JF12b3;  rc_B[2] =  7.1;
  bdisk[3] = JF12b4;  rc_B[3] =  8.3;
  bdisk[4] = JF12b5;  rc_B[4] =  9.8;
  bdisk[5] = JF12b6;  rc_B[5] = 11.4;
  bdisk[6] = JF12b7;  rc_B[6] = 12.7;
  for(size_t i = 0; i < 7; i++)
  {
    bdisk[7] += -fxsec[i]*bdisk[i]/fxsec[7];
  }
  rc_B[7] = 15.5;

  bring = JF12bring;
  hdisk = JF12hdisk;
  wdisk = JF12wdisk;

  Bn = JF12Bn;
  Bs = JF12Bs;
  rn = JF12rn;
  rs = JF12rs;                    // rs > 16.7 kpc, no sigma given
  wh = JF12wh;
  z0 = JF12z0;

  BX = JF12Bx;
  thetaX0 = JF12ThX;
  rXc = JF12rXc;
  rX = JF12rX;

  beta = 1.38 + fsigma*0.41;
  sqrtbeta = sqrt(beta);

  // convert angular units
  spitch=sin(pitch*DEG2RAD);
  cpitch=cos(pitch*DEG2RAD);
  thetaX0 *= DEG2RAD;
  }; //END CONSTRUCTOR


  // CONSTRUCTOR FOR REGULAR+STRIATED FIELD
  JF2012(gsl_rng *r, double sigma, bool stri, bool f, std::string s):BFIELD(),
  rng_ptr(r),bl_use_stri(stri),bl_generate(f),rndmfile(s.c_str())  
  {
  // these parameters are 'constant' within the model
  //  (no sigma values provided)
  pitch = 11.5;
  fsigma = sigma;

  // these parameters are subject to sigma shifting
  bdisk[0] =  0.1 + fsigma*1.8;  rc_B[0] =  5.1;
  bdisk[1] =  3.0 + fsigma*0.6;  rc_B[1] =  6.3;
  bdisk[2] = -0.9 + fsigma*0.8;  rc_B[2] =  7.1;
  bdisk[3] = -0.8 + fsigma*0.3;  rc_B[3] =  8.3;
  bdisk[4] = -2.0 + fsigma*0.1;  rc_B[4] =  9.8;
  bdisk[5] = -4.2 + fsigma*0.5;  rc_B[5] = 11.4;
  bdisk[6] =  0.0 + fsigma*1.8;  rc_B[6] = 12.7;
  bdisk[7] =  2.7 + fsigma*1.8;  rc_B[7] = 15.5;

  bring = 0.1 + fsigma*0.1;
  hdisk = 0.40 + fsigma*0.03;
  wdisk = 0.27 + fsigma*0.08;

  Bn =  1.4 + fsigma*0.1;
  Bs = -1.1 + fsigma*0.1;
  rn = 9.22 + fsigma*0.08;
  rs = 17;                    // rs > 16.7 kpc, no sigma given
  wh = 0.20 + fsigma*0.12;
  z0 = 5.3 + fsigma*1.6;

  BX = 4.6 + fsigma*0.3;
  thetaX0 = 49.0 + fsigma*1;
  rXc = 4.8 + fsigma*0.2;
  rX = 2.9 + fsigma*0.1;

  beta = 1.38 + fsigma*0.41;
  sqrtbeta = sqrt(beta);

  // convert angular units
  spitch=sin(pitch*DEG2RAD);
  cpitch=cos(pitch*DEG2RAD);
  thetaX0 *= DEG2RAD;

  if (bl_use_stri) {
    if (bl_generate) {
      for(size_t i = 0; i < Nsign; i++ )
        mastergrid[i] = GetStriatedSign();
      std::cout << "Dumping striated component to " << rndmfile << std::endl;
      std::ofstream fout;
      fout.open(rndmfile.c_str(), std::ofstream::binary);
      fout.write(mastergrid, Nsign * sizeof(char));
      fout.close();
      std::cout << "Done" << std::endl;
    } else {
      std::cout << "Reading striated component from " << rndmfile << std::endl;
      std::ifstream fin;
      fin.open(rndmfile.c_str(), std::ifstream::binary);
      if (!fin.is_open())
        throw std::runtime_error("Error reading file");
      fin.read (mastergrid, Nsign * sizeof(char));
      fin.close();
      std::cout << "Done" << std::endl;
    } //end striation setup
  }

  };  // END constructor

  virtual ~JF2012(){};    
  
    int GetMasterGridIndex(const double &a, const double &b, const double &c){
      int bin = (int)(10.*(a+20))
                + 400*((int)(10.*(b+20)))
                + 160000*((int)(10.*(c+20)));
      if( (bin<0)||(bin>63999999) ){
        return -1;
      }
      else
        return bin;
    }

  char GetStriatedSign() {
    return gsl_ran_flat(rng_ptr,0,1) > .5 ? 1: -1;
  }

  std::string GetType(void) const
  {
    return "jf2012";
  }

  std::string GetDescription(void) const {
    std::ostringstream ff;
    ff << "F " << GetType() << ' '
       << GetFSigma() << ' ';

    if(bl_use_stri==0){
      ff << 0 << " " << bdisk[0] << " " << bdisk[1] << " " << bdisk[2] << " " 
      << bdisk[3] << " " << bdisk[4] << " " << bdisk[5] << " " << bdisk[6]
      << " " << bring << " " << hdisk << " " << wdisk << " " << Bn
      << " " << Bs << " " << rn << " " << rs << " " << wh
      << " " << z0 << " " << BX << " " << thetaX0 << " " << rXc
      << " " << rX << "\n";
    }
    else{
       ff << bl_use_stri << ' '
          << bl_generate << ' '
          << rndmfile.c_str();
    }
    return ff.str();
  }

  double GetFSigma(void) const {
    return fsigma;
  }

  void gal2fieldgc(double &x, double &y, double &z)
  {
    x = x-solarDistance;
    // y and z remain unchanged
    // this coordinate system is not the same as for other fields
  }
  double JFfieldPhi(const double& x, const double& y)
  {
    return atan2(y,x);
  }
  double L(const double& z, const double& h, const double& w){
    double expo = exp( -2.*(fabs(z)-h) / w );
    return ( 1./(1.+expo) );
  }

  double toroidalhalo(const double& r, const double& z){
    double mag = -1.;
    mag *= exp( -1.*fabs(z) / z0 );
    mag *= L(z, hdisk, wdisk);
    if(z >= 0.)
      mag *= (Bn*(1.-L(r,rn,wh)));
    else
      mag *= (Bs*(1.-L(r,rs,wh)));
    return mag;
  }
  
  // Returns a vector of B magnitudes faster
  void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
    tx = fx;
    ty = fy;
    tz = fz;
    bvec[0] = 0.;
    bvec[1] = 0.;
    bvec[2] = 0.;
    gal2fieldgc(tx,ty,tz);
    tR = BFIELD::fieldR(tx,ty,tz);
    tr = BFIELD::fieldRho(tx,ty);
    if(tR < 1.){                // zero field in 1 kpc GC SPHERE
      return;
    }
    if(withinField(tx,ty,tz)){  // this only checks for r < maxR (=20 kpc)
      tDmag = 0.;
      tHmag = 0.;
      tXmag = 0.;

      tphi = JFfieldPhi(tx,ty);
      ttheta = BFIELD::fieldTheta(tx,ty,tz);
      tcphi = cos(tphi);
      tsphi = sin(tphi);
      tctheta = cos(ttheta);
      tstheta = sin(ttheta);


      // compute toroidal halo magnitude
      tHmag = toroidalhalo(tr,tz);


      // compute X halo mag
      double zsign=1.;
      double rc_x = rXc + fabs(tz)/tan(thetaX0);
      double rp, thetaX;
      if( tr < rc_x ){
        rp = tr * rXc / (rc_x);
        tXmag = BX * pow(rXc/rc_x, 2.) * exp(-1.*rp/rX);
        thetaX = atan2( fabs(tz), (tr-rp) );
        if(tz==0){
          thetaX=pi/2.;
        }
      }
      else{
        rp = tr - fabs(tz)/tan(thetaX0);
        tXmag = BX * (rp/tr) * exp(-1.*rp/rX);
        thetaX = thetaX0;
      }
      if(tz<0){
        zsign *= -1.;
      }


      // assign the halo components
      bvec[0] =     tHmag*tsphi  +  zsign * tXmag*cos(thetaX)*tcphi;
      bvec[1] = -1.*tHmag*tcphi  +  zsign * tXmag*cos(thetaX)*tsphi;
      bvec[2] =                     tXmag*sin(thetaX);


      // compute disk component magnitudes
      if( tr > 3. ){
        if (tr < 5.){   // "molecular ring"
          tDmag = -1.*bring * (5./tr) * (1.-L(tz, hdisk, wdisk));
          bvec[0] =     tDmag*tsphi  + bvec[0];
          bvec[1] = -1.*tDmag*tcphi  + bvec[1];
          // do nothing with bvec[2]
        }
        else{           // spiral region
          double r_negx = tr * exp( -1./tan( (90.-pitch)*DEG2RAD ) * (tphi-pi) );
          if( r_negx > rc_B[7] ){
            r_negx = tr * exp( -1./tan( (90.-pitch)*DEG2RAD ) * (tphi+pi) );
          }
          if( r_negx > rc_B[7] ){
            r_negx = tr * exp( -1./tan( (90.-pitch)*DEG2RAD ) * (tphi+3.*pi) );
          }
          for(int i=7; i>=0; i--){
            if(r_negx < rc_B[i]){
              tDmag = bdisk[i] * (5./tr) * (1.-L(tz, hdisk, wdisk));
            }
          }
          bvec[0] = -1.*(-1.*tDmag*spitch*tcphi + tDmag*cpitch*tsphi)  + bvec[0];
          bvec[1] = (tDmag*spitch*tsphi + tDmag*cpitch*tcphi)  + bvec[1];
          // do nothing with bvec[2]
        } //end else spiral
      } // end r>3

    // incorporate striated component
    if(bl_use_stri){
      double scale = 1. + sqrtbeta * mastergrid[GetMasterGridIndex(tx, ty, tz)];
      bvec[0] = bvec[0] * scale;
      bvec[1] = bvec[1] * scale;
      bvec[2] = bvec[2] * scale;
      } // end if withinFields
    }

  } // end Bvec function
  
  
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
    Bvec(fx,fy,fz,cvec);
    return cvec[2];
  }
  
  
};


#endif
