# Introduction

Welcome to the repository for the modified version of CRT. This page includes instructions for installation and running the code. A full description of the analysis that uses this code is available in this [GAP note] (GAP2016_012). Please note development is in early stages, and things are in a state of flux.

## Installation

The code has been developed and tested on systems running Ubuntu 15 and Debian 8. Assuming one has the necessary dependencies on a Windows or Mac machine, it should run. However, we won't be able to offer any support.

### Dependencies

Before we download and compile, please ensure you have the following libraries/packages:

  - GNU scientific library
```sh
sudo apt-get install libgsl0-dev
```
 - Python 3
  * SciPy

We've attempted to make use of Python's built in packages as much as possible.

### Compiling

The CRT backend is written in C++ and needs to be compiled. We've been able to do this using ```g++``` without any issues.

```sh
make
```

in the root directory.

## Usage

Some minimal setup is required by the user. In order to generate Monte Carlo samples of potential galactic entry positions, the program must be given observed:

 - Auger event ID
 - Arrival directions
 - Energies

in a columnated text file. The ```MC.py``` script is expecting 

  [GAP note]: <https://www.auger.org/index.php/document-centre/finish/139-gap-notes-2016/3994-gap2016-012>