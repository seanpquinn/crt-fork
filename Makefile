#//////////////////////////////////////////////////////////
#// Filename: Makefile
#// Authors:  Brian Baughman        bbaugh@mps.ohio-state.edu
#//
#// Copyright 2007 __MyCompanyName__. All rights reserved.
#//
#// Description: This file compiles the executables 'CRT',
#//                'getfield', and 'test'
#//
#//////////////////////////////////////////////////////////

# Name of the Executable
EXEC     = CRT
# List of directories to include in the search path for header files.
INCLUDES = -I./include

# Set path to look for sources in
VPATH = src ./

# Directory to place object files
BUILDDIR = build


# Source file suffix
SRCSUF   = cpp
SRCS       := $(wildcard src/*.$(SRCSUF))


# Object file suffix
OBJSUF   = o
OBJS       := $(patsubst src/%.$(SRCSUF),$(BUILDDIR)/%.$(OBJSUF),$(SRCS))

# Buildable file suffix
SRCEXSUF = cc
EXSRCS   = $(wildcard *.$(SRCEXSUF))
EXOBJS  := $(patsubst %.$(SRCEXSUF),$(BUILDDIR)/%.$(OBJSUF),$(EXSRCS))


# REMOVE --debug for faster running code.
OPT      = -O2 -Wall
# Command line options for the C++ compiler  
CXXFLAGS = $(OPT) -fPIC $(INCLUDES) \
$(shell gsl-config --cflags)

# Specify the c++ copiler to use
CXX      = g++
# Specify the linker
LD       = g++
# list any libraries needed (other than the standard set in /usr/lib)
LIBS     = $(shell gsl-config --libs)


# Compile all objects
$(BUILDDIR)/%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	/bin/mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c  $< -o $@

$(BUILDDIR)/%.$(OBJSUF) : %.$(SRCEXSUF)
	@echo "<**Compiling**> "$<
	/bin/mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c  $< -o $@


# Compile main
$(EXEC) : $(OBJS) $(BUILDDIR)/$(EXEC).$(OBJSUF)
	/bin/mkdir -p bin
	$(LD) $(LDFLAGS) $(OBJS) $(BUILDDIR)/$(EXEC).$(OBJSUF) $(LIBS) -o ./bin/$@

test : $(OBJS) $(BUILDDIR)/test.$(OBJSUF)
	/bin/mkdir -p bin
	$(LD) $(LDFLAGS) $(OBJS) $(BUILDDIR)/test.$(OBJSUF) $(LIBS) -o ./bin/$@

getfield : $(OBJS) $(BUILDDIR)/getfield.$(OBJSUF)
	/bin/mkdir -p bin
	$(LD) $(LDFLAGS) $(OBJS) $(BUILDDIR)/getfield.$(OBJSUF) $(LIBS)-o ./bin/$@


# build everything
all: $(OBJS) $(EXEC) getfield test

# Remove all objects and built executable
clean:
	/bin/rm -f ./bin/$(EXEC) ./bin/getfield ./bin/test build/*.$(OBJSUF)  core.*


# DO NOT DELETE
