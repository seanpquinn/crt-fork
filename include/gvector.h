//////////////////////////////////////////////////////////
// Filename: gvector.cc
// Authors:  Brian Baughman 
//           Michael Sutherland
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the gvector class and
//                its functions
//
//////////////////////////////////////////////////////////

#ifndef _GVECTOR_H
#define _GVECTOR_H

#include <math.h>

class gvector
  {

  protected:
    double mx;
    double my;
    double mz;
    double ml;
    double mb;
    double mr;
  private:
    bool set;
  public:
    gvector(void):mx(0.),my(0.),mz(0.),ml(0.),mb(0.),mr(0.),set(false){
    };
    gvector(const double& fx,const double& fy,const double& fz)
    :mx(fx),my(fy),mz(fz),ml(0),mb(0),mr(sqrt(fx*fx+fy*fy+fz*fz)),set(false){
    };
    gvector(const double& fl,const double& fb):ml(fl),mb(fb),mr(1.){
      setcart();
    };
    
    virtual ~gvector(void)
    { };
    // Set coords
    void setcoords(const double& fx,const double& fy,const double& fz) {
      mx = fx;
      my = fy;
      mz = fz;
      mr = sqrt(mx*mx+my*my+mz*mz);
      set = false;
    };
    void setangles(void) {
      if(mr>0){
        ml = atan2(my,mx);
        mb = asin(mz/mr);
      } else {
        ml=0.;
        mb=0.;
      }
      set=true;
    }
    void setcart(void) {
      if(mr==0)
        mr=1.;
      mx = mr*cos(mb)*cos(ml);
      my = mr*cos(mb)*sin(ml);
      mz = mr*sin(mb);
      set=true;
    }
    void setcoords(const double& fl,const double& fb) {
      ml = fl;
      mb = fb;
      setcart();
    };
    
    // Returns length of vector
    double r(void) const{
      return mr;
    }
    // Get coords
    double x(void) const{
      return mx;
    }
    double y(void) const{
      return my;
    }
    double z(void) const{
      return mz;
    }
    double l(void) {
      if (set==false)
        setangles();
      return ml;
    }
    double b(void)  {
      if (set==false)
        setangles();
      return mb;
    }
    // vector operators
    void norm(const double& d) {
      if (set==false)
        setangles();        
      if (d!=mr){
        mx = mx*d/mr;
        my = my*d/mr;
        mz = mz*d/mr;
        mr = d;
      }
    };
    double dot(const gvector& ov) const {
      return (mx*ov.mx+my*ov.my+mz*ov.mz);
    };
    double arclen(const gvector& ov) const {
      double x = dot(ov)/(mr*ov.mr);
      if(x>1.)
        x=1.;
      return acos(x);
    };
    
    // arithmetic operators
    gvector operator +(const gvector& ov)
    {
      return gvector(mx+ov.mx,my+ov.my,mz+ov.mz);
    };
    gvector & operator +=(const gvector& ov)
    {
      this->setcoords(mx+ov.mx,my+ov.my,mz+ov.mz);
      return *this;
    };
    gvector operator -(const gvector& ov)
    {
      return gvector(mx-ov.mx,my-ov.my,mz-ov.mz);
    };
    gvector & operator -=(const gvector& ov)
    {
      this->setcoords(mx-ov.mx,my-ov.my,mz-ov.mz);
      return *this;
    };
    bool operator ==(const gvector& ov)
    {
      if (mx==ov.mx && my==ov.my && mz==ov.mz 
        && ml==ov.ml && mb==ov.mb && mr==ov.mr)
        return true;
      else
        return false;
    };
    bool operator !=(const gvector& ov)
    {
      if (mx==ov.mx && my==ov.my && mz==ov.mz 
          && ml==ov.ml && mb==ov.mb && mr==ov.mr)
        return false;
      else
        return true;
    };    
    friend gvector operator*(const double& x,const gvector& ov){
      return gvector(ov.mx*x,ov.my*x,ov.mz*x);
    };
    friend gvector operator*(const gvector& ov,const double& x){
      return gvector(ov.mx*x,ov.my*x,ov.mz*x);
    };
    friend double dot(const gvector& v1, const gvector& v2){
      return v1.dot(v2);
    };
    
  };

#endif /* _GVECTOR_H */



