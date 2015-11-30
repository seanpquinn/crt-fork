//////////////////////////////////////////////////////////
// Filename: bfield_ring.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2009 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the ring (or annulus) magnetic field
//
//////////////////////////////////////////////////////////

#ifndef _ringBfield_H
#define _ringBfield_H
#include "bfield.h"
#include <string.h>

class RingBField: public BFIELD
  {
    
  protected:
    
    double innerEdge; //[kpc]
    double outerEdge;
    double scaleHeight;
    
  public:
    
    RingBField():BFIELD()
    {
    }
    
    virtual ~RingBField()
    {
    }
    
    void SetInnerEdge(const double& s)
    {
      innerEdge = s;
    }
    
    double GetInnerEdge(void) const
    {
      return innerEdge;
    }
    
    void SetOuterEdge(const double& s)
    {
      outerEdge = s;
    }
    
    double GetOuterEdge(void) const
    {
      return outerEdge;
    }
    
    void SetScaleHeight(const double& z)
    {
      scaleHeight = z;
    }
    
    double GetScaleHeight(void) const
    {
      return scaleHeight;
    }
    
    double Bx(const double& fx, const double& fy, const double& fz)
    {
      tx = fx;
      ty = fy;
      tz = fz;
      BFIELD::gal2gc(tx,ty,tz);
      if(withinField(tx,ty,tz)){
        tr = fieldR(tx,ty,tz);
        if( (tr<innerEdge) || (tr>outerEdge) )
          return 0.;
        else
          return ( norm*sin(fieldPhi(tx,ty))*exp(-fabs(tz)/scaleHeight) );
      } else 
        return 0.;
      
    }
    
    double By(const double& fx, const double& fy, const double& fz)
    {
      tx = fx;
      ty = fy;
      tz = fz;
      BFIELD::gal2gc(tx,ty,tz);
      if(withinField(tx,ty,tz)){
        tr = fieldR(tx,ty,tz);
        if( (tr<innerEdge) || (tr>outerEdge) )
          return 0;
        else
          return ( norm*cos(fieldPhi(tx,ty))*exp(-fabs(tz)/scaleHeight) );
      } else 
        return 0.;
    }
    
    double Bz(const double& fx, const double& fy, const double& fz)
    {
      return 0;
    }
    
    void Bvec(const double& fx, const double& fy, const double& fz, double *bvec)
    {
      bvec[2] = 0;
      tx = fx;
      ty = fy;
      tz = fz;
      BFIELD::gal2gc(tx,ty,tz);
      if(withinField(tx,ty,tz)){
        tr = fieldR(tx,ty,tz);
        if( (tr<innerEdge) || (tr>outerEdge) ) {
          bvec[0] = norm*sin(fieldPhi(tx,ty))*exp(-fabs(tz)/scaleHeight);
          bvec[1] = norm*cos(fieldPhi(tx,ty))*exp(-fabs(tz)/scaleHeight);
        }
        else{
          bvec[0] = 0;
          bvec[1] = 0;
        }
      }
      else{
        bvec[0]=0;
        bvec[1]=0;
      }
    }
    
    std::string GetType(void) const 
    {
      return "ring";
    }
    
    std::string GetDescription(void) const {
      std::ostringstream ff;
      ff << "F " <<  GetType() << ' ' 
      << norm << ' '
      << innerEdge << ' '
      << outerEdge << ' '
      << scaleHeight;
      return ff.str();
    }
    
  };


#endif
