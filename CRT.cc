//////////////////////////////////////////////////////////
// Filename: main.cpp
// Authors:  Michael Sutherland    msutherland@phys.lsu.edu
//           Brian Baughman        bbaugh@mps.ohio-state.edu
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the source code for generating
//                the 'CRT' executable, which propagates cosmic rays
//                from source distributions through Galactic magnetic
//                fields
//
//              Information can be found in the accompanying
//                README-CRT file, Auger GAP-2008-099
//                and in Sutherland et al., Astroparticle Physics,
//                Volume 34, Issue 4, p. 198-204. (2010).
//
//////////////////////////////////////////////////////////

// define RESPECTFIELDS if one of your fields changes the step size
#define RESPECTFIELDS
// define TRACKS if you want to output the individual postions of particles at 
// each step
#undef TRACKS
//#define TRACKS 1

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
#include "source.h"
#include "bfield.h"

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

//######################### main program ###########################
int main (int argc, char * const argv[]) {
  // Define constants
  //static const double pi = acos(-1.);
  // parse arguments
  CommandParser arginfo(argc,argv);
  
  
  int nHits=0, nTot=0, evtcount=0; //total number of observed CR events
  double obsL, obsB, srcL, srcB, sep;
  double ox, oy, oz, sx, sy, sz, impactAngle;
  
  std::vector<Source*> srcs;  //source vector
  std::string allFields="";
  combinedfield totalField;   //total magnetic field
  
  
  std::string helpstr = "f";
  helpstr = arginfo.getArg("help,h,H",helpstr,"Flag to display help.");
  
  //#####################################//
  // Initialize the random number generator
  //#####################################//
  gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlux389);
  // Set seed for random number generator
  unsigned long seed=time(NULL);
  std::string seedstr = arginfo.getArg("seed","",
                                       "Option to set seed for random number "
                                       "generator.\nDefault: time.");
  from_string<unsigned long>(seed,seedstr);
  gsl_rng_set(r, seed );
  
  
  
  //#####################################//
  // Choose the integration routine
  //#####################################//  
  bool rk5=false;
  std::string rk5str = arginfo.getArg("rk5",stringify(rk5),
                                      "Flag request the use of stepperdopr5. "
                                      "Superseded by rk853.");
  if(rk5str == "t")
    rk5 = true;
  
  bool rk853=false;
  std::string rk853str = arginfo.getArg("rk853",stringify(rk853),
                                        "Flag to force the use of "
                                        "stepperdopr853."
                                        "Supersedes rk5.");
  if(rk853str == "t")
    rk853 = true;

  if (rk5 && rk853)
    rk5=false;

  if( !rk5 && !rk853) //default to rk853
    rk853 = true;
  
  double err = 1e-8;
  std::string errstr = arginfo.getArg("err",stringify(err),
                                      "Option to set maximum allowed error for "
                                      "dynamic stepping.");
  from_string<double>(err,errstr);
  
  
  
  //#####################################//
  // To backtrk or not to backtrk   (it's better in the original Klingon)
  //#####################################//   
  bool backtrk=false;
  std::string bk = arginfo.getArg("backtrack",stringify(backtrk),
                                  "Flag to forces events to be backtracked. ");
  if(bk == "t")
    backtrk = true;
  
  
  
  //#####################################//
  // Propagation parameters and detector model
  //#####################################//
  Detector* detector = Detector::Instance();
  detector->SetRndPtr(r);
  std::string detRstr = arginfo.getArg("detR","",
                 "Option to set the size of detector in kpc.\nDefault: 0 for "
                 "backtracking and 0.1 for forward tracking.\nDefining a specific detector "
                 "model in your input file overrides this flag"); //In kiloparsecs
  double detR=0.1;
  if(backtrk)
    detR = 0.;
  from_string<double>(detR,detRstr);
  detector->SetDetR(detR);

  std::string dettypestr = arginfo.getArg("dettype","default",
                                          "Option to set the detector type.\n"
                                          "'default' and 'area' available."
  );
  detector->SetType(dettypestr);

  
  double tmax = 1e6;  // In years
  std::string tmaxstr = arginfo.getArg("tmax",stringify(tmax),
                   "Option to set maximum time to track particle in years. "); 

  from_string<double>(tmax,tmaxstr);
  detector->SetTmax(tmax);
  
  double tstep = 10.;
  std::string tstepstr = arginfo.getArg("tstep",stringify(tstep),
                                        "Option to set minimum time step in "
                                        "years."); // In years
  from_string<double>(tstep,tstepstr);
  detector->SetTstep(tstep);

  bool detectall=false;
  std::string dall = arginfo.getArg("detectall",stringify(detectall),
                                    "Flag to attempt to detect number of "
                                    "particles set in a\nsource definition. ");
  if(dall == "t")
    detectall = true;
  
  int tpersrc=10000;
  std::string tpersrcstr = arginfo.getArg("tpersrc",stringify(tpersrc),
                                          "Option to set the number of failed "
                                          "attempts to detect each\nsource "
                                          "before moving on."
  );
  from_string<int>(tpersrc,tpersrcstr);
  
  
  
  //#####################################//
  // Input/Output options
  //#####################################//    
  std::string infile = arginfo.getArg("infile,input","","Option to set input "
                                      "configuration file.\nREQUIRED");
  
  std::string ofile = arginfo.getArg("ofile,output","","Option to set output file.\n"
                                     "Default: std::cout.");
#ifdef TRACKS
  std::string trackfile = arginfo.getArg("trackfile","",
                                         "Option to set file to write tracks "
                                         "to.\nDefault: std::cout.");
#endif
  // ALL arguments should be requested before this point!
  // This ensures a proper usage output
  std::string usagestr = 
	"USEAGE:\n$> "+std::string(argv[0])+" [OPTIONS]\nAvailable options:\n"
	+arginfo.getHelp()+"\n";
  
  
  
  //#####################################//
  // Setup output stream
  //#####################################//
  OutputSrvc* outsrvc = OutputSrvc::Instance();
  if (ofile!="") {
    outsrvc->setOut(ofile);
  }
#ifdef TRACKS
  if (trackfile!="") {
    outsrvc->setLog(trackfile);
  }
  std::ostream Log(outsrvc->getLog());
  Log.precision(10);
#endif
  
  std::ostream out(outsrvc->getOut());
  out.precision(10);
  std::ostream Err(outsrvc->getErr());
  
  //#####################################//
  // Setup intput stream
  //#####################################// 
  std::ifstream in;
  in.open(infile.c_str());
  if(!in.is_open()){
    Err<<"Unable to open input file! "<< infile <<std::endl;
    Err << usagestr;
    return (1);
  } 
  
  // Read input file
  parseinput(detector, r, backtrk, in, totalField, srcs);
  // Close the file
  in.close();
  
  
  const Doub atol=1e-20;
  const Doub rtol=err; 
  const Doub h1=detector->GetTstep();
  const Doub hmin=0.;
  const Doub x1=0.;
  const Doub x2=detector->GetTmax();
  
  gvector detp = detector->GetP0();
  
  
  
  //#####################################//
  //#####################################//
  //#####################################//
  //
  //  The main of main
  //
  // Create headers inside datafiles
  // Add list of used parameters to header
  // Loop over all sources and Generate events
  // Report 
  Err << "Seed is: " << seed << std::endl;
  
  std::vector<Source*>::iterator srciter=srcs.begin();
  Source* src_ptr;
  for(; srciter!=srcs.end(); ++srciter){
    src_ptr = (*srciter);
    nTot += src_ptr->GetNumCRperSrc();
    src_ptr->SetSrcNumber(std::distance(srcs.begin(),srciter));
  }
  srciter=srcs.begin();
	
  if( helpstr.compare("t") == 0 ) {
    Err << usagestr << std::endl;
    exit(1);
  }
	
  // Check if invalid (not requested) arguments have been given
  if (arginfo.checkArgs()==false) {
    std::vector<std::string> argvec =  arginfo.getUnfound();
    std::vector<std::string>::iterator argiter = argvec.begin();
    Err<<"Invalid option(s) given:"<<std::endl;
    for (;argiter!=argvec.end();++argiter) {
      Err<<*(argiter)+ " ";
    }
    Err << std::endl << usagestr << std::endl;
    return (1);
  }
  
  out << "# Run command: " << std::endl << "# ";
  for (int i=0;i<argc;++i) {
    out << argv[i] << " ";
  }
  out << std::endl << "#" << std::endl
  << "# RNGseed: " <<seed<<std::endl
  << "# " <<std::distance(srcs.begin(),srcs.end())
  << " srcs, "<<nTot<<" total events"<<std::endl;
  out << "#\n# " << detector->GetDescription()<<std::endl;
  out << "#"<<std::endl;
  for(; srciter!=srcs.end(); ++srciter){
    src_ptr = (*srciter);
    out << "# " << src_ptr->GetDescription()<<std::endl;
  }
  srciter=srcs.begin();
  out << "#\n" << totalField.GetOutDescription() << "\n#"<<std::endl
  << "# Src Num \t"
  << "CR Energy [EeV] \t"
  << "Src. Long. [deg] \t"
  << "Src. Lat. [deg] \t"
  << "EXG Xcoor. [kpc] \t"
  << "EXG Ycoor. [kpc] \t"
  << "EXG Zcoor. [kpc] \t"
  << "Hit (1=Y,0=N) \t"
  << "Obs. Long. [deg] \t"
  << "Obs. Lat. [deg] \t"
  << "Ang. Sep. [deg] \t"
  << "Imp. Angle [deg] \t"
  << "Imp. Param. [kpc] \t"
  << "Particle Mass \t"
  << "Particle Charge \t"
  << "Evt. Number \t"
  << "Detector Xcoor. [kpc] \t"
  << "Detector Ycoor. [kpc] \t"
  << "Detector Zcoor. [kpc] \t"
  << "Event Quality"
  <<std::endl;
  
#ifdef TRACKS
  Log << "# Run command: " << std::endl << "# ";
  for (int i=0;i<argc;++i) {
    Log << arginfo.chToString(argv[i]) << " ";
  }
  Log << std::endl << "#" << std::endl
  << "# RNGseed: " <<seed<<std::endl
  << "# " <<std::distance(srcs.begin(),srcs.end())
  << " srcs, "<<nTot<<" total events"<<std::endl;
  out << "#\n# " << detector->GetDescription()<<std::endl;
  out << "#"<<std::endl;
  srciter=srcs.begin();
  for(; srciter!=srcs.end(); ++srciter){
    src_ptr = (*srciter);
    Log << "# " << src_ptr->GetDescription()<<std::endl;
  }
  srciter=srcs.begin(); 
  Log << "#\n" << totalField.GetOutDescription() << "#"<<std::endl
  << "# (x,y,z) refers to:   backtracking-> final coordinates,   "
  << "forwardtracking-> injection site" << std::endl
  << "# Event \t"                     
  << "Time [yr] \t"                  
  << "Xcoor. [kpc] \t"
  << "Ycoor. [kpc] \t"
  << "Zcoor. [kpc] \t"
  << "VX [kpc/yr] \t"
  << "VY [kpc/yr] \t"
  << "VZ [kpc/yr] \t"
  << "Energy [EeV] \t"
  << "Event Quality"
  <<std::endl;
#endif
  unsigned long g=1;
  Particle cr;
  gvector crv0 = gvector();
  gvector crv = gvector();
  gvector srcv = gvector();
  VecDoub_IO ystart(NMAX);
  Output outpt(10);
  rhs d;
  int badevt;
  int ntries;
  
  // begin event looping
  for(; srciter!=srcs.end(); ++srciter){
    src_ptr = (*srciter);
    src_ptr->SetDetector(detector);
    ntries=0;
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
      ox = -999;
      oy = -999;
      oz = -999;
      sx = -999;
      sy = -999;
      sz = -999;
      
      // Initialize differential equation array
      ystart[0]=cr.GetP0().x();
      ystart[1]=cr.GetP0().y();
      ystart[2]=cr.GetP0().z();
      ystart[3]=cr.GetV0().x();
      ystart[4]=cr.GetV0().y();
      ystart[5]=cr.GetV0().z();
      d = rhs(cr.GetPrefactor(),&totalField);
      
      // Perform integration routine
      if( rk5){
        Odeint < StepperDopr5 < rhs > > ode5(ystart, x1, x2, atol, rtol, h1, 
                                             hmin,outpt,d);
        badevt = ode5.integrate(&cr, &totalField, src_ptr->GetBktrkStatus(), 
                                detector);
      }
      else if( rk853){
        Odeint < StepperDopr853 < rhs > > ode853(ystart, x1, x2, atol, rtol, h1,
                                                 hmin, outpt,d);
        badevt = ode853.integrate(&cr, &totalField, src_ptr->GetBktrkStatus(), 
                                  detector);
      }
      else{
        Err<<"Incorrect integration routine chosen. "
        << "Please set -rk5 or -rk853 !!"<<std::endl;
        return(1);
      }
      
      if( (detectall)      // enforce 'detectall' condition if true
         && (badevt == 0)
         && (!cr.GetHit()) //  -redo perfectly good events that didn't happen to hit
         && (ntries < tpersrc) ){
        ntries++;
        g--;
        evtcount--;
        continue;	
      }
      
      // Output results and deal with differences between
      //   backtracked and forwardtracked events
      crv0 = cr.GetV0();
      crv = cr.GetV();
      srcv = src_ptr->GetSrcLoc();
      if( src_ptr->GetBktrkStatus()) {//backtracked
        cr.SetCharge(-1.*cr.GetCharge());
        if( cr.GetHit()){
          impactAngle = 90;
          nHits++;  
          obsL = (crv0).l() * RAD2DEG;
          obsB = (crv0).b() * RAD2DEG;
          srcL = (crv).l() * RAD2DEG;
          srcB = (crv).b() * RAD2DEG;
          sep  = (crv0).arclen((crv)) * RAD2DEG;
          ox = cr.GetP0().x();
          oy = cr.GetP0().y();
          oz = cr.GetP0().z();
          sx = cr.GetP().x();
          sy = cr.GetP().y();
          sz = cr.GetP().z();
        } else 
          Err << "Something is wrong!." << std::endl;
      }
      else{
        crv = (-1.*crv);
        srcL = (srcv).l() * RAD2DEG;
        srcB = (srcv).b() * RAD2DEG;
        if( cr.GetHit()){
          impactAngle = 90;
          nHits++;  
          obsL = (crv).l() * RAD2DEG;
          obsB = (crv).b() * RAD2DEG;
          sep  = (srcv).arclen((crv)) * RAD2DEG;
        }
        sx = cr.GetP0().x();
        sy = cr.GetP0().y();
        sz = cr.GetP0().z();
        ox = cr.GetP().x();
        oy = cr.GetP().y();
        oz = cr.GetP().z();
      }
      out << std::distance(srcs.begin(),srciter) << "\t"
      << cr.GetEnergy()/EeV << "\t"
      << srcL << "\t"
      << srcB << "\t"
      << sx << "\t"
      << sy << "\t"
      << sz << "\t"
      << cr.GetHit() << "\t"
      << obsL << "\t"
      << obsB << "\t"
      << sep << "\t"
      << impactAngle << "\t"
      << cr.GetDOCA() << "\t"
      << cr.GetMass()/PROTONMASS << "\t"
      << cr.GetCharge()/PROTONCHARGE << "\t"
      << evtcount << "\t"
      << ox << "\t"
      << oy << "\t"
      << oz << "\t"
      << badevt << std::endl;
    }//end Events loop
  }//end Sources loop
  out << "# Completed Successfully.\n";
  // Close output service
  outsrvc->Close();
  // Free random number generator
  gsl_rng_free(r);
  
  Err<<"# Hits / Total:\t "
  <<nHits<<" / "
  <<nTot<<std::endl
  <<"Wrote: "<< ofile <<std::endl
  <<"All finished"<<std::endl;
  
  return (0);			//successful completion! 
}
//#############################//
//##   End of main  routine ##//
//#############################//

