//////////////////////////////////////////////////////////
// Filename: bfield_bss_s.h
// Authors:  Michael Sutherland
//           Brian Baughman
//
// Copyright 2007 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                the BSS_S spiral magnetic field
//
//////////////////////////////////////////////////////////

#ifndef _BSS_S_H
#define _BSS_S_H

#include "bfield_ssfield.h"
#include <string>

class BSS_S: public SSFIELD
{

 private:


 public:
  BSS_S():SSFIELD(1,false){
    
  };
  virtual ~BSS_S(){
    
  };

 
  std::string GetType(void) const
  {
    return "bss_s";
  };

};


#endif
