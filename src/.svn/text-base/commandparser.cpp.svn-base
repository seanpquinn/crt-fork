//////////////////////////////////////////////////////////
// Filename: commandparser.cpp
// Authors:  Brian Baughman
//
// Copyright 2005 __MyCompanyName__. All rights reserved.
//
// Description: This file contains the implementation of
//                of the command parser class
//
//////////////////////////////////////////////////////////

#include "commandparser.h"
#include <string.h>
#include <iostream>
void CommandParser::parseArgs(void)
{
  if(numArgs_>1) {
    
    std::string::size_type equal;
    std::string::size_type neg;
    std::string cur_arg;
    
    std::string ch_cur_var;
    std::string ch_cur_value;
    
    for(int i=1;i<numArgs_;++i) {
      cur_arg = chToString(parseArgs_[i]);
      equal   = cur_arg.find_first_of('=');
      neg     = cur_arg.find_first_of('-');
      if(equal != std::string::npos && neg!=0) {
        //split the arguments
        ch_cur_value.assign(cur_arg,equal+1,cur_arg.length()-equal-1);
        ch_cur_var.assign(cur_arg,0,equal);
      } else if (neg==0) {
        // -B[it] set some bit option ALWAYS alone!
				if (cur_arg[neg+1]=='-') {
					neg+=1;
				}
        ch_cur_var.assign(cur_arg,neg+1,cur_arg.length()-neg);
        ch_cur_value = "t";
      } else 
        continue;
			
			argMap_.insert(
										 std::pair<std::string,ArgInfo*> 
										 (ch_cur_var,new ArgInfo(ch_cur_var,ch_cur_value,"Unused.")));
			submittedArgs_.push_back(ch_cur_var);
			
    }
  }
}

std::string CommandParser::getHelp(bool showval)
{
  argConstIter_ = argMap_.begin();
  std::ostringstream  retval;
	std::string descrip;
	size_t position = 0;
	size_t replsize = 0;
	std::string replstr = "\n\t";
	std::vector<std::string>::const_iterator nIter = requestedArgs_.begin();
	std::string label;
	for (;nIter!=requestedArgs_.end();++nIter) {
		argConstIter_ = argMap_.find((*nIter));
		descrip = argConstIter_->second->getDescription(showval);
		position = descrip.find( "\n" );
		label = argConstIter_->second->getLabel();
		replsize = label.size()+3;
		replstr = "\n\t";
		for (size_t j=0;j<replsize;++j)
			replstr+=" ";
		while ( position != std::string::npos ) 
		{
			descrip.replace( position, 1, replstr);
			position = descrip.find( "\n", position+replstr.size()+1 );
			if (position==descrip.size()-1)
				position = std::string::npos;
		} 
		retval << '\t' << label << " - " 
		<<   descrip << std::endl << std::endl; 
  }
  return retval.str();
}

std::string CommandParser::getSet()
{
  argConstIter_ = argMap_.begin();
  std::ostringstream  retval;
	std::string descrip;
	std::string replstr = "\n\t";
  for(;argConstIter_!=argMap_.end();++argConstIter_) {
    if (argConstIter_->second->description!="Unused.") {
      retval << '\t' << argConstIter_->second->name << " = " 
      <<   argConstIter_->second->value << std::endl; 
    }
  }
  return retval.str();
}


void CommandParser::checkArg(const std::string& f_arg,
														 const std::string& defaultVal,
                             const std::string& description)
{
	argIter_ = argMap_.find(f_arg);
	if(argIter_ != argMap_.end() ) {
		argIter_->second->setDescription(description);
		argIter_->second->setDefault(defaultVal);
	}
	else {
		argMap_.insert(
									 std::pair<std::string,ArgInfo*> 
									 (f_arg,new ArgInfo(f_arg,defaultVal,description,defaultVal)));
  }
}

std::string CommandParser::getArg(const std::string& f_arg,
																	const std::string& defaultVal,
                                  const std::string& description)
{
	std::vector<std::string> names;
	tokenize(f_arg,names,",;:");
	std::vector<std::string>::const_iterator namesIter = names.begin();
	std::string found = (*(names.begin()));
	for (;namesIter!=names.end();++namesIter) {
		argIter_ = argMap_.find((*namesIter));
		if(argIter_ != argMap_.end() ) {
			found = (*namesIter);
		}
	}
	checkArg(found,defaultVal,description);
	requestedArgs_.push_back((*(names.begin())));
	argIter_ = argMap_.find(found);
	for (namesIter = names.begin();namesIter!=names.end();++namesIter) {
		argIter_->second->addAlt(*namesIter);
		if (found.compare((*namesIter)) != 0) {
			argMap_.insert(std::pair<std::string,ArgInfo*> 
										 ((*namesIter),argIter_->second));
		}
	}
  return argIter_->second->value;
}
std::string CommandParser::chToString(const char* f_ch)
{
  std::string retval;
  if(f_ch==NULL) retval = "";
  else {
    size_t chlen = strlen(f_ch);
    retval.assign(f_ch,chlen);
  }
  return retval;
}

void CommandParser::tokenize(const std::string& str, 
                             std::vector<std::string> &tokens,
                             const std::string& delimiters)
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

bool CommandParser::checkArgs(void)
{
  argConstIter_ = argMap_.begin();
  std::vector<std::string>::iterator reqiter;
  for (;argConstIter_!=argMap_.end();++argConstIter_) {
		if (argConstIter_->second->description.compare("Unused.")==0) {
			unfoundArgs_.push_back(argConstIter_->first);
		}
	}
	
	if (unfoundArgs_.size()>0)
		return false;
	else
		return true;
}
