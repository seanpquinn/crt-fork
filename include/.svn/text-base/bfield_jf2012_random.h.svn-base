//////////////////////////////////////////////////////////
// Filename: bfield_jf2012_random.h
// Authors:  Michael Sutherland
//
// Copyright 2012 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the random component of the
//                Jansson/Farrar magnetic field
//              The methods here do not include the striated
//                component (include/bfield_jf2012.h is a
//                more appropriate location)
//
//////////////////////////////////////////////////////////

#ifndef _JF2012Rndm_H
#define _JF2012Rndm_H

#include "bfield.h"
#include "gvector.h"

#include <vector>
#include <iostream>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

class JF2012Rndm: public BFIELD
  {
    
  private:
    gsl_rng *rng_ptr;      // pointer to rng
    bool bl_generate;      // true == generate cells, otherwise read from file
    std::string rndmfile;

    double b0;             // Sigma for normalization variation
    double b5;
    double bdisk[8];       // field strengths of arms at r=5 kpc
                           //   b8 is determined from other 7
    double rc_B[8];        // radii where each arm crosses the
                           //   negative x-axis
    double pitch; 
    double scale_spiral;
    double scale_radial;   // Scale lengths for exponential radial and 
    double scale_vertical; //  vertical (away from Gal. plane) strength attenuation

    unsigned long Ncells;
    double corrlen;        // Correlation length to use
    double sigmalen;       // Sigma for cell size variation
    class cell {           // An individual cell
      private:
        double curlen;     // Cell "size"
        gvector cellpos;   // Cell position
        gvector curfield;  // Cell field vector
      public:
        cell(){
          curlen = 1.;
          cellpos.setcoords(10., 10., 10.);
          curfield.setcoords(10., 10., 10.);
        };
        virtual ~cell() {};
        void SetCurLen(double a){curlen = a;}
        double GetCurLen(void){return curlen;}
        gvector GetCellPosition(void){return cellpos;}
        void SetCellPosition(gvector a){cellpos = a;}
        gvector GetCellField(void){return curfield;}
        void SetCellField(gvector a){curfield = a;}
    } ;

    cell* cell_ptr;
    std::vector<cell*>::iterator celliter;

    std::vector<cell*> mastergrid[64000];

    gvector curloc;        // particle current location

  public:
    JF2012Rndm(gsl_rng *rn, bool f, unsigned long n, std::string s):BFIELD(),
    rng_ptr(rn),bl_generate(f),rndmfile(s.c_str()),Ncells(n),
    corrlen(0.060),sigmalen(0.025),
    curloc(gvector(0,0,0))
    {
      pitch = 11.5;
      b5 = 7.691;
      b0 = 4.774;
      scale_spiral = 0.607;
      scale_radial = 10.603;
      scale_vertical = 2.511;
      bdisk[0] = 10.637;  rc_B[0] =  5.1;
      bdisk[1] = 6.677;   rc_B[1] =  6.3;
      bdisk[2] = 9.65;    rc_B[2] =  7.1;
      bdisk[3] = 6.911;   rc_B[3] =  8.3;
      bdisk[4] = 1.691;   rc_B[4] =  9.8;
      bdisk[5] = 16.797;  rc_B[5] = 11.4;
      bdisk[6] = 37.057;  rc_B[6] = 12.7;
      bdisk[7] = 10.867;  rc_B[7] = 15.5;

      Init();
    };

    virtual ~JF2012Rndm() {};
    
    void Init(void) {
      if( bl_generate ){          // generate new cells
        std::ofstream out;
        out.open(rndmfile.c_str());
        out << "#Lcorr[kpc]\tXpos[kpc]\tYpos[kpc]\tZpos[kpc]\tBx[muG]\tBy[muG]\tBz[muG]" << std::endl;
        for(unsigned long n=0; n<Ncells; n++){
          cell_ptr = new cell();
          gvector P, B;
          cell_ptr->SetCurLen( sampleCorrelationLength() );
          P = samplePosition();
          cell_ptr->SetCellPosition(P);
          B = sampleField(P);
          cell_ptr->SetCellField(B);
          mastergrid[GetMasterGridIndex(P.x(), P.y(), P.z())].push_back(cell_ptr);
          out << cell_ptr->GetCurLen() << "\t"
              << cell_ptr->GetCellPosition().x() << "\t"
              << cell_ptr->GetCellPosition().y() << "\t"
              << cell_ptr->GetCellPosition().z() << "\t"
              << cell_ptr->GetCellField().x() << "\t"
              << cell_ptr->GetCellField().y() << "\t"
              << cell_ptr->GetCellField().z() << std::endl;
        }
        out.close();
      }
      else{                      // read from file
        std::ifstream in;
        std::string line = "";
        in.open(rndmfile.c_str());
        getline(in, line);
        gvector P;
        for(unsigned long n=0; n<Ncells; n++){
          double x,y,z;
          gvector P;
          cell_ptr = new cell();
          in >> x;
          cell_ptr->SetCurLen(x);
          in >> x >> y >> z;
          P.setcoords(x,y,z);
          cell_ptr->SetCellPosition(P);
          in >> x >> y >> z;
          P.setcoords(x,y,z);
          cell_ptr->SetCellField(P);
          P = cell_ptr->GetCellPosition();
          mastergrid[GetMasterGridIndex(P.x(), P.y(), P.z())].push_back(cell_ptr);

	  // In case the the filesize is less than Ncells,
	  //    don't generate cells with the default values
	  if (in.eof() && line==""){
	    break;
	  }
        }
        in.close();
      }
    }//end Init()


    int GetMasterGridIndex(const double &a, const double &b, const double &c){
      int bin = ((int)(a+20)) + 40*((int)(b+20)) + 1600*((int)(c+20));
      if( (bin<0)||(bin>63999) ){
        return -1;
      }
      else
        return bin;
    }

    double sampleCorrelationLength(void){
      double x;
      do {
        x = corrlen - gsl_ran_gaussian(rng_ptr,sigmalen);
      } while (x<0);
      return x;
    }

    gvector samplePosition(void){
      gvector x;
      x.setcoords(gsl_ran_flat(rng_ptr,-pi,pi),
                         asin(gsl_ran_flat(rng_ptr,-1.,1.)));
      x.norm( maxR*pow(gsl_ran_flat(rng_ptr,0,1),(1./3.)) );
      return x;
    }

    gvector sampleField(gvector p){
      gvector x;
      x.setcoords(gsl_ran_flat(rng_ptr,-pi,pi),
                  asin(gsl_ran_flat(rng_ptr,-1.,1.)));
      // sample spiral magnitude
      double spiralnorm = 0;
      double tr = BFIELD::fieldRho(p.x(),p.y());
      double tz = p.z();
      double tphi = JFfieldPhi(p.x(),p.y());
      if (tr < 5.){   // "molecular ring"
        spiralnorm = norm - gsl_ran_gaussian(rng_ptr,b5);
        do {
          spiralnorm = norm - gsl_ran_gaussian(rng_ptr,b5);
        } while (spiralnorm<0);
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
            spiralnorm = bdisk[i];
          }
        }
        spiralnorm = (norm - gsl_ran_gaussian(rng_ptr,spiralnorm)) * (5./tr);
      }
      spiralnorm = spiralnorm * exp(-0.5*(tz/scale_spiral)*(tz/scale_spiral));
      // sample smooth magnitude
      double smoothnorm = norm - gsl_ran_gaussian(rng_ptr,b0);
      do {
        smoothnorm = norm - gsl_ran_gaussian(rng_ptr,b0);
      } while (smoothnorm<0);
      smoothnorm = smoothnorm * exp(-1.*(tr/scale_radial)) * exp(-0.5*(tz/scale_vertical)*(tz/scale_vertical)) ;
      // add the magnitudes in quadrature and rescale vector
      x.norm( sqrt((smoothnorm*smoothnorm) + (spiralnorm*spiralnorm)) );
      return x;
    }
 
    void gal2fieldgc(double &x, double &y, double &z)
    {
      x = x-solarDistance;
      // y and z remain unchanged
    }

    double JFfieldPhi(const double& x, const double& y)
    {
      return atan2(y,x);
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
    Bvec(fx,fy,fz,cvec);
    return cvec[2];
    }

    void Bvec(const double& fx, const double& fy, const double& fz, double *bvec) {
      tx = fx;
      ty = fy;
      tz = fz;
      bvec[0]=0.;
      bvec[1]=0.;
      bvec[2]=0.;
      gal2fieldgc(tx,ty,tz);
      curloc.setcoords(tx,ty,tz);
      if(withinField(tx,ty,tz)){
        gvector P;
        int bin;
        for(double ix=tx-1; ix<(tx+1.5); ix++){
          for(double iy=ty-1; iy<(ty+1.5); iy++){
            for(double iz=tz-1; iz<(tz+1.5); iz++){
              bin = GetMasterGridIndex(ix, iy, iz);
              if(bin==-1)
                continue;
              std::vector<cell*> cells = mastergrid[bin];
              celliter = cells.begin();
              for(; celliter!=cells.end(); ++celliter){
                cell_ptr = (*celliter);
                P = cell_ptr->GetCellPosition();
                if( (P-curloc).r() < cell_ptr->GetCurLen() ){
                  bvec[0] = bvec[0] + cell_ptr->GetCellField().x();
                  bvec[1] = bvec[1] + cell_ptr->GetCellField().y();
                  bvec[2] = bvec[2] + cell_ptr->GetCellField().z();
                }
              } //end for celliter
            }
          }
        } //end for ix
      }
    } //end Bvec
    


    std::string GetType(void) const 
    {
      return "jf2012_random";
    }
    
    std::string GetDescription(void) const {
      std::ostringstream ff;
      ff << "F " <<  GetType() << " "
         << bl_generate << " "
         << Ncells << " "
         << rndmfile.c_str() << " ";
      
      return ff.str();
    }

  };

#endif
