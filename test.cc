//////////////////////////////////////////////////////////
// Filename: test.cpp
// Authors:  Brian Baughman        bbaugh@mps.ohio-state.edu
//           Michael Sutherland    msutherland@phys.lsu.edu
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the source code for generating
//                the 'test' executable, which propagates 2 cosmic
//                rays through a uniform field for determining the
//                integration routine accuracy and other features
//                of CRT
//
//////////////////////////////////////////////////////////

// Define TEST 
#define TEST 1
// include files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cmath>
#include <math.h>
#include <time.h>
#include <typeinfo>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include "particle.h"
#include "spectrum_monoE.h"
#include "source.h"
#include "source_arbtraj.h"
#include "bfield.h"
#include "bfield_uniform.h"

#include "globals.h"
#include "commandparser.h"
#include "infileparser.h"

#include "nr3.h"
#include "odeint.h"
#include "stepperdopr5.h"
#include "stepperdopr853.h"

#include "OutputSrvc.h"
#include "detector.h"

#include "converter.h"
//*************************** main program ***************************
int main (int argc, char * const argv[]) {

  // parse arguments
  CommandParser arginfo(argc,argv);
  
  // Declare some variables
  int nTot=0, evtcount=0; //total number of observed CR events
  double obsL, obsB, sep, srcL, srcB, x, y, z, impactAngle;
  
  //##################//  
  // Setup test field
  //##################//
  combinedfield totalField;   //total magnetic field
  UniformBField* field =  new UniformBField();
  field->SetMagBx(0.);
  field->SetMagBy(0.);
  field->SetMagBz(1.);
  totalField.addfield(field);
  
  //#####################################//
  // Initialize the random number generator
  //#####################################//
  gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlux389);
  // Set seed for random number generator
  std::string seedstr = arginfo.getArg("seed");
  unsigned long seed=time(NULL);
  from_string<unsigned long>(seed,seedstr);
  gsl_rng_set(r, seed );
  
  
  //#####################################//
  // Choose the integration routine
  //#####################################//  
  bool rk5=false;
  std::string rk5str = arginfo.getArg("rk5");
  if(rk5str == "t")
    rk5 = true;
  
  bool rk853=false;
  std::string rk853str = arginfo.getArg("rk853");
  if(rk853str == "t")
    rk853 = true;
  
  if( !rk5 && !rk853) //default to rk853
    rk853 = true;
  
  
  //#####################################//
  // Propagation parameters
  //#####################################//
  Detector* detector = Detector::Instance();

  double tmax = 1e6;
  std::string tmaxstr = arginfo.getArg("tmax"); // In years
  from_string<double>(tmax,tmaxstr);
  detector->SetTmax(tmax);
  
  double tstep = 10.;
  std::string tstepstr = arginfo.getArg("tstep"); // In years
  from_string<double>(tstep,tstepstr);
  detector->SetTstep(tstep);
  
  double err = 1e-8;
  std::string errstr = arginfo.getArg("err"); // In years
  from_string<double>(err,errstr);
  
  bool backtrk=false;
  
  //#####################################//
  // Input/Output options
  //#####################################//
  std::string infile = arginfo.getArg("infile");
  std::string ofile = arginfo.getArg("ofile");

  //#####################################//
  // Setup output stream
  //#####################################//
  OutputSrvc* outsrvc = OutputSrvc::Instance();
  if (ofile!="") {
    outsrvc->setOut(ofile);
  }
  std::ostream out(outsrvc->getOut());
  out.precision(10);
  std::ostream Err(outsrvc->getErr());
  //##################//  
  // Setup useage
  //##################//
  // ALL arguments should be requested before this point!
  // This ensures a proper usage output
  std::string usagestr = 
	"USEAGE:\n$> "+std::string(argv[0])+" [OPTIONS]\nAvailable options:\n"
	+arginfo.getHelp()+"\n";
  
  //##################//  
  // Setup test sources
  //##################//
  std::vector<Source*> srcs;  //source vector
  double m = 1.*PROTONMASS;
  double q = 1.*PROTONCHARGE;
  
  ArbTrajSrc* src1 = new  ArbTrajSrc(r, new MonoSpectrum(1*EeV),
                                     backtrk, 
                                     10., 0.,0.,
                                     0.,1.,0.,
                                     m, q);
  srcs.push_back(src1);
 
  ArbTrajSrc* src2 = new  ArbTrajSrc(r, new MonoSpectrum(1*EeV),
                                     backtrk, 
                                     10., 0.,0.,
                                     0.,0.,1.,
                                     m, q);

  srcs.push_back(src2);

  

  
  
  //##########################################/
  //##########################################/
  //##########################################/
  //
  //  Test integration routines
  //
  // Create headers inside datafiles
  // Add list of used parameters to header
  // Loop over all sources and Generate events
  Err << "\nSeed is: " << seed << std::endl;
  std::vector<Source*>::iterator srciter=srcs.begin();
  Source* src_ptr;
  for(; srciter!=srcs.end(); ++srciter){
    src_ptr = (*srciter);
    nTot += src_ptr->GetNumCRperSrc();
  }
  srciter=srcs.begin();

  out << "# Run command: " << std::endl << "# ";
  for (int i=0;i<argc;++i) {
    out << argv[i] << " ";
  }
  out << std::endl << "#" << std::endl
  << "# RNGseed: " <<seed<<std::endl
  << "# " <<std::distance(srcs.begin(),srcs.end())
  << " srcs, "<<nTot<<" total events, detector size="
  << detector->GetDetR()<<"kpc"<<std::endl;
  for(; srciter!=srcs.end(); ++srciter){
    src_ptr = (*srciter);
    out << "# " << src_ptr->GetDescription()<<std::endl;
  }
  srciter=srcs.begin();
  out << "#\n" << totalField.GetOutDescription() << "#"<<std::endl
  << "# (x,y,z) refers to:   backtracking-> final coordinates,   "
  << "forwardtracking-> injection site" << std::endl
  << "# Step \t"
  << "Time [yr] \t"
  << "Xcoor. [kpc] \t"
  << "Ycoor. [kpc] \t"
  << "Zcoor. [kpc] \t"
  << "Xexact. [kpc] \t"
  << "Yexact. [kpc] \t"
  << "Zexact. [kpc] \t"
  << "VX [kpc/yr] \t"
  << "VY [kpc/yr] \t"
  << "VZ [kpc/yr] \t"
  << "VXexact [kpc/yr] \t"
  << "VYexact [kpc/yr] \t"
  << "VZexact [kpc/yr] \t"
  << "Particle Mass \t"
  << "Particle Charge"
  <<std::endl;
  
  //##################//  
  // Declare some vars
  //##################//
  unsigned long g=1;
  Particle cr;
  rhs d;
  Output outpt(10);
  VecDoub_IO ystart(NMAX);
  //##################//  
  // Setup some constants
  //##################//
  const Doub atol=1e-50; //err;
  const Doub rtol=err; //atol;
  const Doub h1=detector->GetTstep();
  const Doub hmin=0.;
  const Doub x1=0.;
  const Doub x2=detector->GetTmax();
  
  // begin event looping
  for(; srciter!=srcs.end(); ++srciter){
    src_ptr = (*srciter);
    for(g=1;g<src_ptr->GetNumCRperSrc()+1;++g){   //event counter
      src_ptr->getParticle(cr);
      evtcount++;
      
      // Initialize output variables
      impactAngle = -999;
      obsL = -999;
      obsB = -999;
      srcL = -999;
      srcB = -999;
      sep = -999;
      x = -999;
      y = -999;
      z = -999;
      
      // Initialize differential equation array
      ystart[0]=cr.GetP0().x();
      ystart[1]=cr.GetP0().y();
      ystart[2]=cr.GetP0().z();
      ystart[3]=cr.GetV0().x();
      ystart[4]=cr.GetV0().y();
      ystart[5]=cr.GetV0().z();
      d=rhs(cr.GetPrefactor(),&totalField);

      // Perform integration routine
      if(rk5){
        Odeint < StepperDopr5 < rhs > > ode5(ystart, x1, x2, atol, rtol, h1, hmin,
                                        outpt,d);
        ode5.integrate(&cr, &totalField, src_ptr->GetBktrkStatus(), detector);
      }
      else if(rk853){
        Odeint < StepperDopr853 < rhs > > ode853(ystart, x1, x2, atol, rtol, h1,
                                            hmin, outpt,d);
        ode853.integrate(&cr, &totalField, src_ptr->GetBktrkStatus(), detector);
      }
      else{
        Err<<"Incorrect integration routine chosen. "
        << "Please set -rk5 or -rk853 !!"<<std::endl;
        return(1);
      }
      
    }//end Events loop
  }//end Sources loop
  
  
  // Close output service
  outsrvc->Close();
  Err <<"Wrote: "<< ofile <<std::endl
  <<"All finished"<<std::endl;

  return (0);			//successful completion! 
}
//#############################//
//##   End of main  routine ##//
//#############################//

