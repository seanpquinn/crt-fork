//////////////////////////////////////////////////////////
// Filename: infileparser.h
// Authors:  Brian Baughman 
//           Michael Sutherland
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the infileparser class, which handles the
//                reading and interpretation of the input
//                configuration file as well as creating the
//                user-requested source and magnetic fields
//
//////////////////////////////////////////////////////////

#ifndef _INFILEPARSER_H
#define _INFILEPARSER_H

#include <fstream>
#include <vector>
#include <gsl/gsl_rng.h>

#include "source.h"
#include "bfield_combinedfield.h"
#include "bfield_ssfield.h"


// Fucntion to parse input
void parseinput(Detector* det, gsl_rng *r,bool backtrk, std::ifstream &in, 
                combinedfield &totalField, std::vector<Source*> &srcs);
// Function to trim leading and trailing whitespace from a string
std::string trim(std::string str);
// Function to remove all whitespace from string
std::string rmwhite(std::string str);
// Function to remove double spaces and replace tabs with single spaces
std::string clean(std::string istr);
// Functions to check and set parameters for spiral magnetic fields
bool CheckSpiralParameters(double p, double l1, double l2, double l3,
			 double l4, double zcut);
void SetSpiralParameters(SSFIELD *bf, double p, double l1, double l2,
			 double l3, double l4, double zcut, double N);

#endif
