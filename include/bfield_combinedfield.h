//////////////////////////////////////////////////////////
// Filename: bfield_combinedfield.h
// Authors:  Brian Baughman
//           Michael Sutherland
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the total magnetic field, comprised of the
//                individual field components listed in the
//                input configuration file
//
//////////////////////////////////////////////////////////

#ifndef _combinedfield_H
#define _combinedfield_H

#include "bfield.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

class combinedfield: public BFIELD
  {
    
  protected:
    BFIELD** fields;
    unsigned int nfields;
    
  public:
    combinedfield():BFIELD(),fields(0),nfields(0)
    {};
    virtual ~combinedfield()
    {};
    
    void addfield(BFIELD *newfield) {
      // Iterate number of fields
      ++nfields;
      // Allocate space for another field
      fields = (BFIELD**) realloc (fields, nfields * sizeof(BFIELD*));
      if (fields == NULL) {
        // Fields are messed up!!
        // This should do something important
      }
      fields[nfields-1]=newfield;
      if(newfield->GetMaxR()>maxR)
        maxR=newfield->GetMaxR();
      if(newfield->GetNstps()>Nstps)
        Nstps=newfield->GetNstps();
    }
    unsigned int checkStep(const gvector& P,const gvector& V, double &step) {
      unsigned int retval=0;
      for(unsigned int i=0;i<nfields;++i) {
        retval += fields[i]->checkStep(P,V,step);
      }
      return retval;
    }
    double Bx(const double& fx, const double& fy, const double& fz)
    {
      tx=0.;
      for(unsigned int i=0;i<nfields;++i) {
        tx += fields[i]->Bx(fx,fy,fz);
      }
      return tx;
    }
    
    double By(const double& fx, const double& fy, const double& fz)
    {
      ty=0.;
      for(unsigned int i=0;i<nfields;++i) {
        ty += fields[i]->By(fx,fy,fz);
      }
      return ty;
    }
    
    double Bz(const double& fx, const double& fy, const double& fz)
    {
      tz=0.;
      for(unsigned int i=0;i<nfields;++i) {
        tz += fields[i]->Bz(fx,fy,fz);
      }
      return tz;
    }
    void Bvec(const double& fx, const double& fy, const double& fz, double *bvec)
    {
      bvec[0]=0;
      bvec[1]=0;
      bvec[2]=0;
      for(unsigned int i=0;i<nfields;++i) {
        fields[i]->Bvec(fx,fy,fz,cvec);
        bvec[0]+=cvec[0];
        bvec[1]+=cvec[1];
        bvec[2]+=cvec[2];
      }
    }
    std::string GetType(void) const 
    {
      std::string a = "";
      for(unsigned int i=0;i<nfields;++i) {
        a += fields[i]->GetType();
        a += ", ";
      }
      return a;
    }
    
    std::string GetDescription(void) const {
      std::string a = "";
      for(unsigned int i=0;i<nfields;++i) {
        a += fields[i]->GetDescription();
        a += "\n";
      }
      return a;
    }
    std::string GetOutDescription(void) const {
      std::string a = "";
      for(unsigned int i=0;i<nfields;++i) {
        a += "# ";
        a += fields[i]->GetDescription();
        if (i!=nfields-1)
          a += "\n";
      }
      return a;
    }
  };


#endif
