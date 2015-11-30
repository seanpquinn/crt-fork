//////////////////////////////////////////////////////////
// Filename: infileparser.cpp
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

#include "infileparser.h"

#include <iostream>
#include "spectrum.h"
#include "spectrum_monoplaw.h"
#include "spectrum_brokenplaw.h"
#include "spectrum_monoE.h"

#include "source_point.h"
#include "source_longitudinal.h"
#include "source_latitudinal.h"
#include "source_isotropic.h"
#include "source_cr.h"
#include "source_arbtraj.h"
#include "source_eg.h"
#include "source_crsmear.h"

#include "bfield.h"
#include "bfield_bss_a.h"
#include "bfield_bss_s.h"
#include "bfield_ass_a.h"
#include "bfield_ass_s.h"
#include "bfield_uniform.h"
#include "bfield_dipole.h"
#include "bfield_combinedfield.h"
#include "bfield_simplerndm.h"
#include "bfield_ring.h"
#include "bfield_toroid.h"
#include "bfield_rN.h"
#include "bfield_sun2008.h"
#include "bfield_pshirkov_disk.h"
#include "bfield_pshirkov_halo.h"
#include "bfield_jf2012.h"
#include "bfield_jf2012_random.h"

#include "globals.h"
#include "OutputSrvc.h"
#include "detector.h"

void parseinput(Detector* det, gsl_rng *r, bool backtrk, std::ifstream &in,
                combinedfield &totalField, std::vector<Source*> &srcs) 
{ 
  // Get output
  OutputSrvc* outsrvc = OutputSrvc::Instance();
  std::ostream err(outsrvc->getErr());
  
  std::string first = "";
  std::string stype = "";
  std::string fieldtype = "";
  std::string line = "";
  std::string model="";
  std::stringstream ssline;

  // Source parameters
  unsigned long Np;
  double C, index, emin, emax, d, tE, tL, tB, beamangle, m, q, px, py, pz, vx, vy, vz,
    Brand, Lcorr, Leg, index1, ebreak, Eres, Ares;

  // Generic field parameters
  double norm, pitch, s1, s2, s3, s4, zc, bx, by, bz;

  // Sun2008 field parameters
  double p, Bc, Rc, B0, R0, z0, R1, R2, R3;
  double p1, p2, Rb, Rb1, Rb2;
  double BH0, zH0, zH1, zH1a, RH0;
  
  // JF2012 field parameters and flags
  // Init flags
  bool flag, stri;
  // Disk parameters
  double JF12b1,JF12b2,JF12b3,JF12b4,JF12b5,JF12b6,JF12b7;
  double JF12bring,JF12hdisk,JF12wdisk;
  // Toroidal halo parameters
  double JF12Bn,JF12Bs,JF12rn,JF12rs,JF12wh,JF12z0;
  // X halo
  double JF12Bx,JF12ThX,JF12rXc,JF12rX;

  Source *src_ptr;
  std::string tstring;
  while( !in.eof() ){
    // Only deal with whole lines
    getline(in,line);
    if (in.eof() && line=="") 
      break;
    if (line=="")    // Skip blank lines
      continue;
    ssline.clear();
    ssline.str(trim(line));
    first = ssline.peek();
    if( first == "#"){//skip headers and comments
      continue;
    }    
    ssline >> stype;
    if (ssline.fail()) {
      err << "Invalid Line!:" 
	  << line << std::endl;
      continue;
    }
    tstring="";
    
    //#############################//
    //  detector parameterization
    if(stype=="D"){
      ssline >> px >> py >> pz >> C;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Detector Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      if(C < 0.){
        err<<"Invalid detector size request."<<std::endl;
        err<<"Detector model will default to detR = 0.1 kpc."<<std::endl;
        C = 0.1;
      }
      det->SetP0(px,py,pz);
      det->SetDetR(C);
      tstring=det->GetDescription();
    }
    //#############################//
    //  source parameterizations
    else if(stype=="I"){
      ssline >> Np >> C >> index >> emin >> emax >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      if(index == -1){
        err<<"Unable to calculate spectrum; "
	   << "spectral index = (-1) which introduces divergences. "
	   << "This source will be skipped."<<std::endl;
        continue;
      }  
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new IsotropicSrc(Np, C, r, 
                                 new PLawSpectrum(emin*EeV,emax*EeV,index),
                                 backtrk, 0., m, q, 0);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
    }//end if 'i'type
    else if(stype=="X"){
      ssline >> Np >> C >> index >> emin >> emax >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      if(index == -1){
        err<<"Unable to calculate spectrum; "
	   << "spectral index = (-1) which introduces divergences. "
	   << "This source will be skipped."<<std::endl;
        continue;
      }  
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new IsotropicSrc(Np, C, r, 
                                 new PLawSpectrum(emin*EeV,emax*EeV,index),
                                 backtrk, 0., m, q, 1);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
    }//end if 'x'type
    else if(stype == "P"){
      ssline >> Np >> C >> index >> emin >> emax >> d >> tL >> tB >> beamangle 
	     >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      if(index == -1){
        err<<"Unable to calculate spectrum; "
	   << "spectral index = (-1) which introduces divergences. "
	   << "This source will be skipped."<<std::endl;
        continue;
      }
      beamangle *= DEG2RAD;
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      if(d<=0)
        d=SOURCEDISTANCE;
      tL *= DEG2RAD;
      tB *= DEG2RAD;
      src_ptr = new  PointSrc(Np, C, r,
                              new PLawSpectrum(emin*EeV,emax*EeV,index)
                              , backtrk,d, tL, tB, beamangle, m, q);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
      
    }//end if 'p'type
    else if(stype=="L"){
      ssline >> Np >> C >> tE >> tL >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      //C *= KPC;
      tL *= DEG2RAD;
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new LongSrc(Np, C, r,
                            new MonoSpectrum(tE*EeV),
                            tL,backtrk,0., m, q, 0);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
      
    }//end if 'l'type
    else if(stype=="Y"){
      ssline >> Np >> C >> tE >> tL >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      tL *= DEG2RAD;
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new LongSrc(Np, C, r,
                            new MonoSpectrum(tE*EeV),
                            tL,backtrk,0., m, q, 1);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
      
    }//end if 'y'type
    else if(stype=="B"){
      ssline >> Np >> C >> tE >> tB >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
      }
      tB *= DEG2RAD;
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new LatSrc(Np, C, r,
                           new MonoSpectrum(tE*EeV),
                           tB,backtrk,0., m, q, 0);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
      
    }//end if 'b'type
    else if(stype=="Z"){
      ssline >> Np >> C >> tE >> tB >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
      }
      tB *= DEG2RAD;
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new LatSrc(Np, C, r,
                           new MonoSpectrum(tE*EeV),
                           tB,backtrk,0., m, q, 1);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
      
    }//end if 'z'type
    else if(stype=="C"){
      if (backtrk!=true) {
        err << "CRSrcs are only usable with -backtrack. Ignoring source." 
	    << line << std::endl;
        continue;
      }
      ssline >> tE >> tL >> tB >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      tL*=DEG2RAD;
      tB*=DEG2RAD;
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new CRSrc(0.,r,new MonoSpectrum(tE*EeV),true,tL,tB,m,q);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
    }//end of type 'c' 
    else if(stype=="CS"){
      if (backtrk!=true) {
        err << "CRSrcSmears are only usable with -backtrack. Ignoring source." 
	    << line << std::endl;
        continue;
      }
      ssline >> tE >> tL >> tB >> m >> q >> Np >> Eres >> Ares;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      tL*=DEG2RAD;
      tB*=DEG2RAD;
      Ares*=DEG2RAD;
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new CRSrcSmear(Np,0.,r,new MonoSpectrum(tE*EeV),true,tL,tB,m,q, Eres,Ares);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
    }//end of type 'cs' 
    else if(stype == "A"){
      ssline >> tE >> px >> py >> pz >> vx >> vy >> vz
	     >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:"
	    << line << std::endl;
        continue;
        
      }
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new  ArbTrajSrc(r, new MonoSpectrum(tE*EeV),
                                backtrk,
                                px,  py,  pz,
                                vx,  vy,  vz,
                                m,  q);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
      
    }//end if 'A' type
    else if(stype == "E"){
      ssline >> Np >> C >> index >> emin >> emax >> d >> tL >> tB >>  Brand 
	     >> Lcorr >> Leg >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      if(index == -1){
        err<<"Unable to calculate spectrum; "
	   << "spectral index = (-1) which introduces divergences. "
	   << "This source will be skipped."<<std::endl;
        continue;
      }
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      if(d <0)
        d=SOURCEDISTANCE;
      tL *= DEG2RAD;
      tB *= DEG2RAD;
      src_ptr = new  EGSrc(Np, C, r,
                           new PLawSpectrum(emin*EeV,emax*EeV,index)
                           , backtrk,d, tL, tB, Brand,Lcorr,Leg, m, q);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
    }//end if 'e'type
    else     if(stype=="BPI"){
      ssline >> Np >> C >> index >> index1 >> ebreak >> emin >> emax >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      if(index == -1 || index1==-1){
        err<<"Unable to calculate spectrum; "
	   << "spectral index = (-1) which introduces divergences. "
	   << "This source will be skipped."<<std::endl;
        continue;
      }  
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new IsotropicSrc(Np, C, r, 
                                 new BPLawSpectrum(emin*EeV,emax*EeV,
                                                   index,index1,ebreak),
                                 backtrk, 0., m, q, 0);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
    }//end if 'BPI'type
    else     if(stype=="BPX"){
      ssline >> Np >> C >> index >> index1 >> ebreak >> emin >> emax >> m >> q;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Source Definition!:" 
	    << line << std::endl;
        continue;
        
      }
      if(index == -1 || index1==-1){
        err<<"Unable to calculate spectrum; "
	   << "spectral index = (-1) which introduces divergences. "
	   << "This source will be skipped."<<std::endl;
        continue;
      }  
      m *= PROTONMASS;
      q *= PROTONCHARGE;
      src_ptr = new IsotropicSrc(Np, C, r, 
                                 new BPLawSpectrum(emin*EeV,emax*EeV,
                                                   index,index1,ebreak),
                                 backtrk, 0., m, q, 1);
      tstring=src_ptr->GetDescription();
      srcs.push_back(src_ptr);
    }//end if 'BPX'type
    
    //#############################//
    //  magnetic field parameterizations
    else if(stype=="F"){
      ssline >> fieldtype;
      if (ssline.fail()) {
        ssline.clear();
        err << "Invalid Field Definition!:" 
	    << line << std::endl;
        continue;
      }
      if(fieldtype=="ass_a"){
        ssline >> model >> norm >> pitch >> s1 >> s2 >> s3 >> s4 >> zc;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
          
        }
        ASS_A *field = new ASS_A();
        if(model=="stanev")
          field->SetModel(0);
        else if(model=="hmr")
          field->SetModel(1);
        else{
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        
        if(CheckSpiralParameters(pitch, s1, s2, s3, s4, zc)){
          ssline.clear();
          err << "Invalid parameter value(s)!:" 
	      << line << std::endl;
          continue;
        }
        SetSpiralParameters(field, pitch, s1, s2, s3, s4, zc, norm);
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="ass_s"){
        ssline >> model >> norm >> pitch >> s1 >> s2 >> s3 >> s4 >> zc;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
          
        }
        ASS_S *field = new ASS_S();
        if(model=="stanev")
          field->SetModel(0);
        else if(model=="hmr")
          field->SetModel(1);
        else{
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        
        if(CheckSpiralParameters(pitch, s1, s2, s3, s4, zc)){
          ssline.clear();
          err << "Invalid parameter value(s)!:" 
	      << line << std::endl;
          continue;
        }
        SetSpiralParameters(field, pitch, s1, s2, s3, s4, zc, norm);
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="bss_a"){
        ssline >> model >> norm >> pitch >> s1 >> s2 >> s3 >> s4 >> zc;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
          
        }
        BSS_A *field = new BSS_A();
        if(model=="stanev")
          field->SetModel(0);
        else if(model=="hmr")
          field->SetModel(1);
        else{
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        
        if(CheckSpiralParameters(pitch, s1, s2, s3, s4, zc)){
          ssline.clear();
          err << "Invalid parameter value(s)!:" 
	      << line << std::endl;
          continue;
        }
        SetSpiralParameters(field, pitch, s1, s2, s3, s4, zc, norm);
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="bss_s"){
        ssline >> model >> norm >> pitch >> s1 >> s2 >> s3 >> s4 >> zc;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
          
        }
        BSS_S *field = new BSS_S();
        if(model=="stanev")
          field->SetModel(0);
        else if(model=="hmr")
          field->SetModel(1);
        else{
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        
        if(CheckSpiralParameters(pitch, s1, s2, s3, s4, zc)){
          ssline.clear();
          err << "Invalid parameter value(s)!:" 
	      << line << std::endl;
          continue;
        }
        SetSpiralParameters(field, pitch, s1, s2, s3, s4, zc, norm);
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="uniform"){
        ssline >> bx >> by >> bz;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
          
        }
        UniformBField* field =  new UniformBField();
        field->SetMagBx(bx);
        field->SetMagBy(by);
        field->SetMagBz(bz);
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if (fieldtype=="dipole"){
        ssline >> norm;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
          
        }
        DipoleBField* field = new DipoleBField();
        field->SetNorm(norm*baseDipoleNorm);
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="simprand"){
        ssline >> norm >> zc >> s1 >> s2;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
          
        }
        SimpleRndmField* field =  new SimpleRndmField(r);
        field->SetNorm(norm);
        field->SetCorrLength(zc);
        field->SetSigmaLen(s1);
        field->SetSigmaNorm(s2);
        field->Init();
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="ring"){
        ssline >> norm >> s1 >> s2 >> zc;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        RingBField* field = new RingBField();
        field->SetNorm(norm);
        field->SetInnerEdge(s1);
        field->SetOuterEdge(s2);
        field->SetScaleHeight(zc);
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="toroidal"){
        ssline >> s1 >> zc >> s2 >> norm;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        ToroidalBField* field = new ToroidalBField();
        field->SetToroidalRadius(s1);
        field->SetHeight(zc);
        field->SetHalfWidth(s2);
        field->SetBmax(norm);
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="r^N"){
        ssline >> norm >> index >> s1 >> pitch >> model >> zc;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        R2BField* field = new R2BField();
        field->SetNorm(norm);
	field->SetIndex(index);
        field->SetInnerEdge(s1);
	field->SetPitchAngle(pitch*DEG2RAD);
	if(model=="none"){
	  field->SetModel(0);
	  field->SetScaleHeight(1.);
	}
	else if(model=="expo"){
	  field->SetModel(1);
	  field->SetScaleHeight(zc);
	}
	else if(model=="smooth"){
	  field->SetModel(2);
	  field->SetScaleHeight(zc);
	}
	else{
	  ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
	}
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="sun2008"){
        ssline >> model;
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        SUN2008 *field = new SUN2008();
        if(model=="ass+ring"){
          ssline >> p >> B0 >> R0 >> Bc >> Rc >> z0 >> R1 >> R2 >> R3 >> BH0 >> zH0 >> zH1 >> zH1a >> RH0;
          if (ssline.fail()) {
            ssline.clear();
            err << "Invalid Field Definition!:" 
		<< line << std::endl;
            continue;
          }
          field->SetModel(1);
          field->SetPitchAngle(p*DEG2RAD);
          field->SetR0(R0);
          field->SetBc(Bc);
          field->SetRc(Rc);
          field->Setz0(z0);
          field->SetR1(R1);
          field->SetR2(R2);
          field->SetR3(R3);
          field->SetzH0(zH0);
          field->SetzH1(zH1);
          field->SetzH1a(zH1a);
          field->SetRH0(RH0);
          field->SetNorm(B0);
          field->SetHaloNorm(BH0);
        } else if(model=="bss") {
          ssline >> B0 >> R0 >> Bc >> Rc >> z0 >> Rb >> p1 >> Rb1 >> p2 >> Rb2 >> BH0 >> zH0 >> zH1 >> zH1a >> RH0;
          if (ssline.fail()) {
            ssline.clear();
            err << "Invalid Field Definition!:" 
		<< line << std::endl;
            continue;
          }
          field->SetModel(0);
          field->SetR0(R0);
          field->SetBc(Bc);
          field->SetRc(Rc);
          field->Setz0(z0);
          field->SetRb(Rb);
          field->SetRb1(Rb1);
          field->SetRb2(Rb2);
          field->SetPitchAngle1(p1*DEG2RAD);
          field->SetPitchAngle2(p2*DEG2RAD);
          field->SetzH0(zH0);
          field->SetzH1(zH1);
          field->SetzH1a(zH1a);
          field->SetRH0(RH0);
          field->SetNorm(B0);
          field->SetHaloNorm(BH0);
        } else {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if(fieldtype=="pshirkov-disk"){
	ssline >> model;
	if (ssline.fail()) {
	  ssline.clear();
	  err << "Invalid Field Definition!:" 
	      << line << std::endl;
	  continue;
	}
	PSHIRKOVDISK *field = new PSHIRKOVDISK();
		  
	if(model=="bss_s"){
	  ssline >> B0 >> p >> Rc >> d >> z0;
	  if (ssline.fail()) {
	    ssline.clear();
	    err << "Invalid Field Definition!:" 
		<< line << std::endl;
	    continue;
	  }
	  field->SetModel(1);
	  field->SetPitchAngle(p*DEG2RAD);
	  field->SetRc(Rc);
	  field->Setz0(z0);
	  field->Setd(d);
	  field->SetNorm(B0);
	  tstring=field->GetDescription();
	  totalField.addfield(field);
	}
	else if(model=="ass_s") {
	  ssline >> B0 >> p >> Rc >> d >> z0;
	  if (ssline.fail()) {
	    ssline.clear();
	    err << "Invalid Field Definition!:" 
		<< line << std::endl;
	    continue;
	  }
	  field->SetModel(0);
	  field->SetPitchAngle(p*DEG2RAD);
	  field->SetRc(Rc);
	  field->Setz0(z0);
	  field->Setd(d);
	  field->SetNorm(B0);
	  tstring=field->GetDescription();
	  totalField.addfield(field);
	}
	  
      }
      else if(fieldtype=="pshirkov-halo"){ssline >> model;
	if (ssline.fail()) {
	  ssline.clear();
	  err << "Invalid Field Definition!:" 
	      << line << std::endl;
	  continue;
	}
	PSHIRKOVHALO *field = new PSHIRKOVHALO();
		  
	if(model=="sq"){
	  ssline >> BH0 >> zH0 >> zH1 >> zH1a >> RH0;
	  if (ssline.fail()) {
	    ssline.clear();
	    err << "Invalid Field Definition!:" 
		<< line << std::endl;
	    continue;
	  }
	  field->SetModel(1);
	  field->SetHaloNorm(BH0);
	  field->SetzH0(zH0);
	  field->SetzH1(zH1);
	  field->SetzH1a(zH1a);
	  field->SetRH0(RH0);
	  tstring=field->GetDescription();
	  totalField.addfield(field);
	}
	else if(model=="abs") {
	  ssline >> BH0 >> zH0 >> zH1 >> zH1a >> RH0;
	  if (ssline.fail()) {
	    ssline.clear();
	    err << "Invalid Field Definition!:" 
		<< line << std::endl;
	    continue;
	  }
	  field->SetModel(0);
	  field->SetHaloNorm(BH0);
	  field->SetzH0(zH0);
	  field->SetzH1(zH1);
	  field->SetzH1a(zH1a);
	  field->SetRH0(RH0);
	  tstring=field->GetDescription();
	  totalField.addfield(field);
	}
		  
      }
      else if (fieldtype=="jf2012"){
        ssline >> index1 >> stri;
        JF2012* field;
        if(stri==0){
          ssline >> JF12b1 >> JF12b2 >> JF12b3 >> JF12b4 >> JF12b5 >>
          JF12b6 >> JF12b7 >> JF12bring >> JF12hdisk >> JF12wdisk >>
          JF12Bn >> JF12Bs >> JF12rn >> JF12rs >> JF12wh >> JF12z0 >>
          JF12Bx >> JF12ThX >> JF12rXc >> JF12rX;
          field = new JF2012(r,
                      index1,
                      stri,
                      JF12b1,
                      JF12b2,
                      JF12b3,
                      JF12b4,
                      JF12b5,
                      JF12b6,
                      JF12b7,
                      JF12bring,
                      JF12hdisk,
                      JF12wdisk,
                      JF12Bn,
                      JF12Bs,
                      JF12rn,
                      JF12rs,
                      JF12wh,
                      JF12z0,
                      JF12Bx,
                      JF12ThX,
                      JF12rXc,
                      JF12rX);
        }
        else{
          ssline >> flag >> tstring;
          field = new JF2012(r, index1, stri, flag, tstring.c_str());
        }
        if (ssline.fail()) {
          ssline.clear();
          err << "Invalid Field Definition!:" 
	      << line << std::endl;
          continue;
        }
        tstring=field->GetDescription();
        totalField.addfield(field);
      }
      else if (fieldtype=="jf2012_random"){
	ssline >> flag >> Np >> tstring;
	ssline.clear();
        JF2012Rndm* field = new JF2012Rndm(r, flag, Np, tstring.c_str());
        tstring=field->GetDescription();
        totalField.addfield(field);
      }

    }//end if stype==F

    // deal with incorrect input lines
    else{
      err << "Incorrect source or field read from file! Ignoring line!"
	  << std::endl;
      continue;
    }
    std::string rmline = clean(line);
    if (tstring!="" && tstring!=rmline) {
      err << "Description mismatch!" << std::endl
	  << "---" << tstring << "---" << std::endl
	  << "---" << rmline << "---" << std::endl
	  << "Please check that this input entry matches the entry used by CRT." 
	  << std::endl;
    }
    first="";
    stype="";
    fieldtype="";
  }//end while EOF  
  
}

std::string trim(std::string str)
{
  char const* delims = " \t\r\n";
  // trim leading whitespace
  std::string::size_type  notwhite = str.find_first_not_of(delims);
  str.erase(0,notwhite);
  
  // trim trailing whitespace
  notwhite = str.find_last_not_of(delims); 
  str.erase(notwhite+1);
  return str;
}

std::string rmwhite(std::string str)
{
  char const* delims = " \t";
  // trim leading whitespace
  std::string::size_type  white = str.find_first_of(delims);
  while(white!=std::string::npos && str.length()!=1) {
    str.erase(white,1);
    white = str.find_first_of(delims);
  }
  
  return str;
}

std::string clean(std::string str)
{
  char const* t = "\t";
  // trim leading whitespace
  std::string::size_type  tab = str.find_first_of(t);
  while(tab!=std::string::npos && str.length()!=1) {
    str.replace(tab,1," ");
    tab = str.find_first_of(t);
  }
  char const* ds = "  ";
  // trim leading whitespace
  std::string::size_type  dspace = str.find(ds);
  while(dspace!=std::string::npos && str.length()!=1) {
    str.replace(dspace,2," ");
    dspace = str.find(ds);
  }
  return str;
}

bool CheckSpiralParameters(double p, double l1, double l2,
                           double l3, double l4, double zcut)
{
  if( (p==0)
      || (l1<=0)
      || (l2<=0)
      || (l3<=0)
      || (l4<=0)
      || (zcut<=0))
    return true; // a bad parameter
  else
    return false;
}

void SetSpiralParameters(SSFIELD *bf, double p, double l1, double l2,
                         double l3, double l4, double zcut, double N){
  bf->SetPitchAngle(p*DEG2RAD);
  bf->SetScale1(l1);
  bf->SetScale2(l2);
  bf->SetScale3(l3);
  bf->SetScale4(l4);
  bf->SetZcut(zcut);
  bf->SetNorm(N);
}
