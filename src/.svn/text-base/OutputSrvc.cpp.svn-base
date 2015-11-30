//////////////////////////////////////////////////////////
// Filename: OutputSrvc.cpp
// Authors:  Brian Baughman 
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the output service class
//
//////////////////////////////////////////////////////////

#include "OutputSrvc.h"

OutputSrvc* OutputSrvc::pinstance = 0;// initialize pointer
OutputSrvc* OutputSrvc::Instance () 
{
  if (pinstance == 0)  // is it the first call?
  {  
    pinstance = new OutputSrvc; // create sole instance
  }
  return pinstance; // address of sole instance
}


OutputSrvc::OutputSrvc() 
{
  // Set all Outputs to standard Outputs
  ofile_Log="";
  setLog(std::cout.rdbuf());
  ofile_Err="";
  setErr(std::cerr.rdbuf());
  ofile_Out="";
  setOut(std::cout.rdbuf());
}

void OutputSrvc::Close(void) 
{
  std::vector<std::ofstream*>::iterator ofiter=openfiles.begin();
  for(;ofiter!=openfiles.end();++ofiter) {
    (*ofiter)->flush();
    (*ofiter)->close();
  }
}

bool OutputSrvc::setLog(std::streambuf* strbuf) {
  Log = strbuf;
  return true;
}
bool OutputSrvc::setLog(std::string ofile) {
  Log=openfile(ofile);
  if (Log!=0)
    return true;
  else {
    setLog(std::cout.rdbuf());
    return false;
  }
    
}
bool OutputSrvc::setErr(std::streambuf* strbuf) {
  Err=strbuf;
  return true;
}
bool OutputSrvc::setErr(std::string ofile) {
  Err=openfile(ofile);
  if (Err!=0)
    return true;
  else {
    setErr(std::cerr.rdbuf());
    return false;
  }
}
bool OutputSrvc::setOut(std::streambuf* strbuf) {
  Out=strbuf;
  return true;
}
bool OutputSrvc::setOut(std::string ofile) {
  Out=openfile(ofile);
  if (Out!=0)
    return true;
  else {
    setOut(std::cout.rdbuf());
    return false;
  }
}
std::streambuf* OutputSrvc::openfile(std::string ofile) {
  /**************************/
  // Setup Output stream
  /**************************/
  openfiles.push_back(new std::ofstream());
  std::ofstream* ofilestream = openfiles.at(openfiles.size()-1);
  // use file as Output
  ofilestream->open(ofile.c_str(),std::ios::out);
  if( !ofilestream ) {
    std::cerr << "Error opening Output file stream: " 
    << ofile << std::endl;
    return 0;
  }
  
  return ofilestream->rdbuf();
}






