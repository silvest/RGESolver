# Install script for directory: /home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xheaderx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/RGESolver" TYPE FILE FILES
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/Expanded.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_complex.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_function_adapter.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_matrix_base.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_matrix_complex.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_matrix_double.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_rgerunner.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_special_functions.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_vector_base.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_vector_complex.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/gslpp_vector_double.h"
    "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/src/std_make_vector.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/oem/EvolutoreStandalone/Programmi/V10_3/RGESolver/gslpp/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
