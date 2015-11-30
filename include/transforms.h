//////////////////////////////////////////////////////////
// Filename: transforms.h
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

#ifndef _TRANSFORMS_H
#define _TRANSFORMS_H

#include "globals.h"

namespace transforms {
  // Define transformation constants
  const double GdNGP=27.128310056550401*DEG2RAD;
  const double Ga0=(192.8595212503885+90.)*DEG2RAD;
  const double GL0=(122.93193411101866-90.)*DEG2RAD;
  const double SGdNGP=15.64407736*DEG2RAD;
  const double SGa0=(283.18940711+90.)*DEG2RAD;
  const double SGL0=(26.73153707-90.)*DEG2RAD;
  
  void gal2eq(const double &l, const double &b,
              double &ra, double &dec,bool supergalactic=false);  
  
  void eq2gal(const double &ra, const double &dec,
              double &l, double &b,bool supergalactic=false);  
};

#endif

