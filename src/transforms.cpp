//////////////////////////////////////////////////////////
// Filename: transforms.cpp
// Authors:  Brian Baughman
//           Michael Sutherland
//
// Copyright 2009 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementations for
//                transforming between galactic and equatorial
//                coordinate systems
//
//////////////////////////////////////////////////////////

#include "transforms.h"
#include <math.h>

void transforms::gal2eq(const double &l, const double &b,
            double &ra, double &dec,bool supergalactic) {
  double dNGP=GdNGP;
  double a0=Ga0;
  double L0=GL0;
  if (supergalactic==true) {
    dNGP=SGdNGP;
    a0=SGa0;
    L0=SGL0;
  }
  dec=asin(cos(b)*cos(dNGP)*sin(l-L0)+sin(b)*sin(dNGP));
  double x=cos(b)*cos(l-L0);
  double y=(-sin(b)*cos(dNGP)+cos(b)*sin(dNGP)*sin(l-L0));
  ra=atan2(y,x)+a0;
  if (ra>pi)
    ra-=pi*2;
  
}

void transforms::eq2gal(const double &ra, const double &dec,
            double &l, double &b,bool supergalactic) {
  double dNGP=GdNGP;
  double a0=Ga0;
  double L0=GL0;
  if (supergalactic==true) {
    dNGP=SGdNGP;
    a0=SGa0;
    L0=SGL0;
  }
  b=asin(-cos(dec)*cos(dNGP)*sin(ra-a0)
         +sin(dec)*sin(dNGP));
  double y=(sin(dec)*cos(dNGP)+cos(dec)*sin(dNGP)*sin(ra-a0));
  double x=cos(dec)*cos(ra-a0);
  l=(atan2(y,x)+L0);
  if (l>pi)
    l-=pi*2;
}
