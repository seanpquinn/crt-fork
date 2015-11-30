//////////////////////////////////////////////////////////
// Filename: converter.h
// Authors:  Brian Baughman 
//
// Copyright 2008 __MyCompanyName__. All rights reserved.
//
// Description: This file contains basic codes found on the web
//                and in the commons
//
//////////////////////////////////////////////////////////

//#ifndef _H_CONVERTER
//#define _H_CONVERTER 1
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
//#include <stdexcept>

// Class thrown when conversion error occurs
class BadConversion {
  //: public std::runtime_error {
public:
  std::string message;
  BadConversion(const std::string& s) : message(s)
  // : std::runtime_error(s)
  { }
  
};
// Function to convert any type understandable by stringstream to a string
// Usage: 
// commontype mytype=....;
// std::string mystring = stringify(mytype);
template<typename T>
inline std::string stringify(const T& x)
{
  std::ostringstream o;
  if (!(o << x))
    throw BadConversion(std::string("stringify(")
                        + typeid(x).name() + ")");
  return o.str();
}
// Function to convert a std::string to any type understandable by stringstream 
// Usage: 
// std::string mystring = "8675309";
// commontype mytype;
// convert(mystring,mytype,true);
// The last arguement is a boolean which tells convert to fail if the conversion
// worked but had leftover characters
template<typename T>
inline void convert(const std::string& s, T& x,
                    bool failIfLeftoverChars = true)
{
  std::istringstream i(s);
  char c;
  if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
    throw BadConversion(s);
}
// Return value function to convert a std::string to any type understandable by 
// stringstream. This is just an interface to convert
// Usage: 
// std::string mystring = "8675309";
// commontype mytype = convert<commontype>(mystring,true);
// The last arguement is a boolean which tells convert to fail if the conversion
// worked but had leftover characters
template<typename T>
inline T convertTo(const std::string& s,
                   bool failIfLeftoverChars = true)
{
  T x;
  convert(s, x, failIfLeftoverChars);
  return x;
}
// Function to converts a std::string to any type understandable by 
// stringstream but does NOT throw an exception if the conversion fails.
// Usage: 
// std::string mystring = "8675309";
// commontype mytype = convert<commontype>(mystring,true);
// The third arguement allows the user to define the ios_base so convertion to
// hex or binary can be accomplished
// The last arguement is a boolean which tells convert to fail if the conversion
// worked but had leftover characters
template <class T>
bool from_string(T& x, const std::string& s
                 ,std::ios_base& (*f)(std::ios_base&)=std::dec,
                 bool failIfLeftoverChars = true)
{
  std::istringstream i(s);
  char c;
  if (!(i >> f >> x) || (failIfLeftoverChars && i.get(c)))
    return false;
  else
    return true;
};

//#endif //_H_CONVERTER

