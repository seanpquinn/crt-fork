//////////////////////////////////////////////////////////
// Filename: bfield_ass_s.h
// Authors:  Michael Sutherland
//           Brian Baughman 
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the ASS_S spiral magnetic field
//
//////////////////////////////////////////////////////////

#ifndef _ASS_S_H
#define _ASS_S_H

#include "bfield_ssfield.h"
#include <string>

class ASS_S: public SSFIELD
{

 private:


 public:
  ASS_S():SSFIELD(1,true) {

  };
  virtual ~ASS_S(){
    
  };

  
  std::string GetType(void) const
  {
    return "ass_s";
  };

};


#endif
