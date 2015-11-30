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

#ifndef _OUTPUTSRVC_H
#define _OUTPUTSRVC_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class OutputSrvc 
{
public:
  static OutputSrvc* Instance();
  bool setLog(std::streambuf* strbuf);
  bool setLog(std::string ofile);
  bool setErr(std::streambuf* strbuf);
  bool setErr(std::string ofile);
  bool setOut(std::streambuf* strbuf);
  bool setOut(std::string ofile);
  std::streambuf* getLog(void) {
    return Log;
  }
  void getLog(std::ostream& fLog) {
    fLog.rdbuf(Log);
  }
  std::streambuf* getErr(void) {
    return Err;
  }
  void getErr(std::ostream& fErr) {
    fErr.rdbuf(Err);
  }
  std::streambuf* getOut(void) {
    return Out;
  }
  void getOut(std::ostream& fOut) {
    fOut.rdbuf(Out);
  }
  void Close(void);
  
protected:
  OutputSrvc();
  OutputSrvc(const OutputSrvc&);
  OutputSrvc& operator= (const OutputSrvc&);
  // Log streambuf
  std::streambuf* Log;
  // Error streambuf
  std::streambuf* Err;
  // Standard streambuf
  std::streambuf* Out;
private:
  static OutputSrvc* pinstance;
  std::string ofile_Log;
  std::string ofile_Err;
  std::string ofile_Out;
  std::streambuf* openfile(std::string ofile);
  std::vector<std::ofstream*> openfiles;
  
};
#endif //_OUTPUTSRVC_H

