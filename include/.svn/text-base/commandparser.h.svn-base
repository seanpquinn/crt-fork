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

#ifndef _H_COMMANDPARSER
#define _H_COMMANDPARSER 1
#include <map>
#include <string>
#include <vector>
#include <sstream>

class ArgInfo 
{
public:
  ArgInfo(void){
  };
  ArgInfo(std::string fname,std::string fvalue="", std::string fdescription="",
					std::string fdefaultVal="") :
	name(fname), 
	value(fvalue), 
	description(fdescription), 
	defaultVal(fdefaultVal) {
    
  };
  std::string name;
  std::string value;
  std::string description;
  std::string defaultVal;
	std::vector<std::string> altNames;
  void setDefault(std::string fdefaultVal) {
    defaultVal=fdefaultVal;
  };
	void setValue(std::string fVal) {
    value=fVal;
  };
  void setDescription(std::string fdescription) {
    description=fdescription;
  };
	std::string getDescription(bool showval=false) const {
		if (defaultVal.compare("")==0)
			return description;
		else 
			if (showval)
				return description+"\nDefault: "+defaultVal+"\nValue: "+value;
			else
				return description+"\nDefault: "+defaultVal;
		
  };
	void addAlt(std::string fAlt) {
    altNames.push_back(fAlt);
  };
	std::string getLabel(void) {
		std::string rv;
		std::vector<std::string>::const_iterator iter = altNames.begin();
		for(;iter!=altNames.end();++iter) {
			if (iter!=altNames.begin())
				rv+=";";
			rv+=(*iter);
		}
		return rv;
	};
};

class CommandParser
{
public:
  CommandParser(int nargs,  char * const args[]): 
	numArgs_(nargs), parseArgs_(args)
  {
    parseArgs();
  };
  virtual ~CommandParser(void){};
  std::string getArg(const std::string& f_arg,
										 const std::string& defaultVal="",
                     const std::string& description="");
  std::vector< std::string > getRequeted(void) {
    return requestedArgs_;
  };
  std::vector< std::string > getUnfound(void) {
    return unfoundArgs_;
  };
  std::string getHelp(bool showval=false);
  std::string getSet(void);
  std::string chToString(const char* f_ch);
  bool checkArgs(void);
private:
  int numArgs_;
  char* const *parseArgs_;
  std::map< std::string , ArgInfo* > argMap_;
	std::map<std::string,ArgInfo*>::const_iterator argConstIter_;
	std::map<std::string,ArgInfo*>::iterator argIter_;
  std::vector< std::string > requestedArgs_;
  std::vector< std::string > submittedArgs_;
  std::vector< std::string > unfoundArgs_;
	void parseArgs(void);
  void checkArg(const std::string& f_arg,
		const std::string& defaultVal="",
                const std::string& description="");
  void tokenize(const std::string& str, std::vector<std::string> &tokens,
                const std::string& delimiters = " ");
	
};



#endif // _H_COMMANDPARSER

