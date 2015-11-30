//////////////////////////////////////////////////////////
// Filename: getfield.cc
// Authors:  Brian Baughman 
//           Michael Sutherland
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the source code for generating
//                the 'getfield' executable, which returns magnetic
//                field values on a grid
//
//////////////////////////////////////////////////////////


// Define VECCHECK if you want errors indicating differences between access 
// methods of BFIELD
#undef VECCHECK
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

#include "bfield_combinedfield.h"
#include "bfield.h"
#include "bfield_bss_a.h"
#include "bfield_bss_s.h"
#include "bfield_ass_a.h"
#include "bfield_ass_s.h"
#include "bfield_uniform.h"
#include "bfield_dipole.h"
#include "bfield_simplerndm.h"
#include "bfield_ring.h"
#include "bfield_toroid.h"
#include "bfield_rN.h"
#include "bfield_sun2008.h"
#include "bfield_pshirkov_disk.h"
#include "bfield_pshirkov_halo.h"
#include "bfield_jf2012.h"
#include "bfield_jf2012_random.h"
#include <gsl/gsl_rng.h>
#include "source.h"

#include "globals.h"
#include "commandparser.h"
#include "infileparser.h"

#include "OutputSrvc.h"
#include "detector.h"

#include "converter.h"
//*************************** main program ***************************
int main (int argc, char * const argv[]) {
  // Define constants
  //static const double pi = acos(-1.);
  // parse arguments
  CommandParser arginfo(argc,argv);
  
  
  combinedfield totalField;   //total magnetic field
  //#########################//
  // Setup scan grid
  //#########################//
  double maxx=60;
  double minx=-60;
  long nx=120*2+1;
  from_string<double>(minx,arginfo.getArg("minx"));
  from_string<double>(maxx,arginfo.getArg("maxx"));
  from_string<long>(nx,arginfo.getArg("nx"));
  double maxy=60;
  double miny=-60;
  long ny=120*2+1;
  from_string<double>(miny,arginfo.getArg("miny"));
  from_string<double>(maxy,arginfo.getArg("maxy"));
  from_string<long>(ny,arginfo.getArg("ny"));
  double maxz=60;
  double minz=-60;
  long nz=120*2+1;
  from_string<double>(minz,arginfo.getArg("minz"));
  from_string<double>(maxz,arginfo.getArg("maxz"));
  from_string<long>(nz,arginfo.getArg("nz"));

  //#####################################//
  // Input/Output options
  //#####################################//    
  std::string infile = arginfo.getArg("infile");
  std::string ofile = arginfo.getArg("ofile");
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
  std::ostream out(outsrvc->getOut());
  out.precision(10);
  std::ostream Err(outsrvc->getErr());
  std::ostream Log(outsrvc->getLog());

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
  Detector* detector = Detector::Instance();
  std::vector<Source*> srcs;  //source vector
  gsl_rng *r = gsl_rng_alloc(gsl_rng_ranlux389);
  unsigned long seed=time(NULL);
  gsl_rng_set(r, seed );
  parseinput(detector,r,false, in, totalField,srcs);
  // Close the file
  in.close();


  //#####################################//
  //#####################################//
  //
  //  The main of getfield
  //
  // Create headers inside datafiles
  // Add list of used parameters to header
  // Loop over all sources and Generate events
  out << "# Run command: " << std::endl << "# ";
  for (int i=0;i<argc;++i) {
    out << argv[i] << " ";
  }
  out << std::endl 
    << "#" << detector->GetDetR()<<"kpc"<<std::endl
    << totalField.GetOutDescription() <<std::endl
    << "#"
    << "X [kpc] \t"
    << "Y [kpc] \t"
    << "Z [kpc] \t"
    << "B_X [muG] \t"
    << "B_Y [muG] \t"
    << "B_Z [muG] \t"
    << "B_tot [muG] \t"
    <<std::endl;

  // begin event looping
  double x,y,z;
  double dx = (maxx-minx)/(nx-1.);
  if (nx==1)
    dx=0.;
  double dy = (maxy-miny)/(ny-1.);
  if (ny==1)
    dy=0.;  
  double dz = (maxz-minz)/(nz-1.);
  if (nz==1)
    dz=0.;
  
  double bvec[3];
  for(int i=0;i<nx;++i){
    x = minx+dx*double(i);
    for(int j=0;j<ny;++j){
      y = miny+dy*double(j);
      for(int k=0;k<nz;++k){
        z = minz+dz*double(k);
        totalField.Bvec(x,y,z,bvec);
        out
        << x << "\t"
        << y << "\t"
        << z << "\t"
        << bvec[0] << "\t"
        << bvec[1] << "\t"
        << bvec[2] << "\t"
        << sqrt(bvec[0]*bvec[0] + bvec[1]*bvec[1] + bvec[2]*bvec[2]) << "\t"
        << std::endl;
#ifdef VECCHECK
        // Check of fields
        double bvec[3];
        totalField.Bvec(x,y,z,bvec);
        if (fabs(totalField.Bx(x,y,z)-bvec[0])>1e-8)
          Err << "Bad x component!" << bvec[0] << '\t' << totalField.Bx(x,y,z)<< std::endl;
        if (fabs(totalField.By(x,y,z)-bvec[1])>1e-8)
          Err << "Bad y component!" << bvec[1] << '\t' << totalField.By(x,y,z)<< std::endl;
        if (fabs(totalField.Bz(x,y,z)-bvec[2])>1e-8)
          Err << "Bad z component!" << bvec[2] << '\t' << totalField.Bz(x,y,z)<< std::endl;
#endif
      }// end z loop
    }//end y loop
  }//end x loop
  // Close output service
  outsrvc->Close();
  // Free random number generator
  Log << "Wrote: "<< ofile <<std::endl
  <<"All finished"<<std::endl;
  return (0);			//successful completion! 
}
//#############################//
//##   End of main  routine ##//
//#############################//

