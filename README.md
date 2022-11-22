
# `RGESolver`

A `C++` library to perform renormalization group evolution of SMEFT coefficients numerically. A faster, approximate solution that neglects the scale dependence of the anomalous dimension matrix is also available.
  The general flavour case at dimension-six level is considered. Operators that violate lepton and/or baryon number conservation are not considered. The documentation for this library can be found [here](https://silvest.github.io/RGESolver/annotated.html)
  
`RGESolver` is a free software under the copyright of the GNU General Public License.

## Dependencies

* `BOOST`  : BOOST is a C++ library which can be obtained from the [`BOOST` website](https://www.boost.org/) or from Linux package managers or Mac ports. RGESolver only requires the BOOST headers, not the full libraries, so a header-only installation is sufficient.
* `GSL` : The GNU Scientific Library (GSL) is a C library for numerical computations. It can be found on the [`GSL` website](https://www.gnu.org/software/gsl/). Most Linux package managers will have a stable version as will any ports for Mac. 
* `C++11` : A compiler that supports at least `C++11` standard is required.
## Installation
The installation of `RGESolver` requires the availability of `CMake` in the system (version `3.1` or greater). A description of `CMake` and the instructions for its installation can be found in the [`CMake`website](https://cmake.org/).
Clone the repository with 
```
git clone https://github.com/silvest/RGESolver --recursive
```
The installation can be performed writhing the following lines in a terminal session (in the `RGESolver` directory):
```
mkdir build && cd $_
cmake .. <options>
cmake --build .
cmake --install .
```
Note that depending on the setting of installation prefix (see below) the user might need root privileges to be able to install `RGESolver`.

### Command line options for the installation

* `-DLOCAL_INSTALL:BOOL=<ON or OFF>` : to install `RGESolver` in the directory `build/install` (default: `OFF`).
* `-DCMAKE_INSTALL_PREFIX:PATH=<RGESolver installation directory>` : the directory in which `RGESolver`	will be installed (default: `/usr/local`). This variable cannot be modified when `-DLOCAL INSTALL ALL=ON` is set.
* `-DDEBUG_MODE:BOOL=<ON or OFF>` : to enable the debug mode (default: `OFF`).
* `-DBOOST_INCLUDE_DIR:PATH=<include path>/boost/` : `CMake`checks for `BOOST` headers availability in the system and fails if they are not installed. Thus, if  `BOOST` is not installed 	in the predefined search path, the user can specify where it is with this option. The path must end with the `boost/`directory which contains the headers.
* `-DGSL_CONFIG_DIR:PATH=<gsl-config directory>` :  `RGESolver` uses `gsl-config` to get the `GSL` parameters. If this is not in the predefined search path, the user can specify it with this option.

## Usage

The `rgesolver-config` script is available in the `<CMAKE_INSTALL_PREFIX>/bin` directory (default: `/usr/local`), which can be invoked with the following options:
* `--cflags`: to obtain the include path needed for compilation against the `RGESolver`.
* `--libs`:  to obtain the flags needed for linking against the `RGESolver`.

If the path `<CMAKE_INSTALL_PREFIX>/bin` is not in the predefined search path, the compilation will (most likely) fail. if the user wants to use the compilation command above, it is suggested to add `<CMAKE_INSTALL_PREFIX>/bin` to the `$PATH` variable. 
Alternatively, the script can be invoked from a terminal session in `<CMAKE_INSTALL_PREFIX>/bin` to visualize the paths to the library and to the headers.

After the installation, the example program `Example1.cpp` (available in the `Examples` directory) can be compiled with the command 
 ```
 g++ -o app Example1.cpp `rgesolver-config --cflags` `rgesolver-config --libs`
 ```

## Uninstall

The user can uninstall the library typing in a terminal session in the `build` directory:
 ```
cmake --build . --target uninstall
 ```
Also in this case, depending on the setting of installation prefix, the user might need root privileges to be able to uninstall `RGESolver`.


