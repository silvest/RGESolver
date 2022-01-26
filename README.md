# `RGESolver`
A `C++` library to perform renormalization group evolution of SMEFT coefficients, both numerically and with the leading-log approximation. 
  The general flavour case at dimension-six level is considered. Operators that violate lepton and/or baryon number conservation are not considered.
## Dependencies
* `BOOST`  : BOOST is a C++ library which can be obtained from the BOOST website or from Linux package managers or Mac ports. RGESolver only requires the BOOST headers, not the full libraries, so a header-only installation is sufficient.
* `GSL` : The GNU Scientific Library (GSL) is a C library for numerical computations. It can be found on the GSL website. Most Linux package managers will have a stable version as will any ports for Mac. 


## Installation

The installation can be performed with:
```
mkdir build && cd $_
cmake ..
make
make install
```
### Command line options
`-DLOCAL_INSTALL:BOOL=<ON or OFF>`
`-DCMAKE_INSTALL_PREFIX:PATH=<RGESolver installation directory>`
`-DDEBUG_MODE:BOOL=<ON or OFF>`
`-DBOOST_INCLUDE_DIR:PATH=<include path>/boost/`
`-DGSL_CONFIG_DIR:PATH=<gsl-config directory>`

 Note that depending on the setting of installation prefix you might need root privileges to be able to install `RGESolver` with `sudo make install` instead of `just make install`.
