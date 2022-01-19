cmake_minimum_required(VERSION 3.10)
project(SMEFTRunner VERSION 1.0 DESCRIPTION "Some library")

add_library(SMEFTRunner RGESolver.cc)

  
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/SMEFTRunner.pc.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/SMEFTRunner.pc" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/smeftrunner-config.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/smeftrunner-config" @ONLY)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/smeftrunner-config ${CMAKE_CURRENT_BINARY_DIR}/SMEFTRunner.pc 
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR} COMPONENT smeftrunner-config)

#add_executable(RGESolver )
#target_link_libraries()

#include(GNUInstallDirs)

# specify the C++ standard
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED True)

##########  BOOST  ##########

# find_package(Boost) is not used.

set(BOOST_ERROR_MESSAGE "Boost headers could not be found. Please specify the Boost include path with ${RedBold}-DBOOST_INCLUDE_DIR=<include path>/boost/${ColReset}.")

if(NOT BOOST_INCLUDE_DIR)
  find_path(BOOST_INCLUDE_DIR "boost/tokenizer.hpp" HINTS ENV C_INCLUDE_PATH ENV CPLUS_INCLUDE_PATH DOC "Boost include path")
  if(BOOST_INCLUDE_DIR)
    set(BOOST_INCLUDE_DIR "${BOOST_INCLUDE_DIR}/boost/" CACHE PATH "Boost include path" FORCE)
  else()
    message(FATAL_ERROR "${BOOST_ERROR_MESSAGE}")
  endif()
else()
  set(BOOST_INCLUDE_DIR ${BOOST_INCLUDE_DIR} CACHE PATH "Boost include path" FORCE)
  find_path(TMP_BOOST_INC "tokenizer.hpp" PATHS ${BOOST_INCLUDE_DIR} NO_DEFAULT_PATH)
  if(NOT TMP_BOOST_INC)
    unset(TMP_BOOST_INC CACHE)
    message(FATAL_ERROR "${BOOST_ERROR_MESSAGE}")
  else()
    unset(TMP_BOOST_INC CACHE)
  endif()
endif()

message(STATUS "Found Boost headers in ${BOOST_INCLUDE_DIR}")
string(REGEX REPLACE "/boost$" "" BOOST_INC ${BOOST_INCLUDE_DIR})




##########  GSL  ##########

if(NOT GSL_CONFIG_DIR)
  find_program(GSL_CONFIG gsl-config)
  if(GSL_CONFIG)
    string(REGEX REPLACE "/gsl-config$" "" TMP_PATH ${GSL_CONFIG})
    set(GSL_CONFIG_DIR ${TMP_PATH} CACHE PATH "Path to gsl-config")
    unset(TMP_PATH)
  endif()
else()
  set(GSL_CONFIG_DIR ${GSL_CONFIG_DIR} CACHE PATH "Path to gsl-config" FORCE)
  find_program(GSL_CONFIG gsl-config PATHS ${GSL_CONFIG_DIR} NO_DEFAULT_PATH)
endif()
if(NOT GSL_CONFIG)
  unset(GSL_CONFIG CACHE)
  message(FATAL_ERROR "gsl-confg could not be found. Please specify the directory of gsl-config with ${RedBold}-DGSL_CONFIG_DIR=<directory>${ColReset}.")
else()  
  message(STATUS "Found gsl-config in ${GSL_CONFIG_DIR}")
  execute_process(COMMAND ${GSL_CONFIG} --version
                  OUTPUT_VARIABLE GSL_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "GSL version: ${GSL_VERSION}")
  execute_process(COMMAND ${GSL_CONFIG} --cflags
                  OUTPUT_VARIABLE GSL_CFLAGS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "GSL cflags: ${GSL_CFLAGS}")
  execute_process(COMMAND ${GSL_CONFIG} --libs
                  OUTPUT_VARIABLE GSL_LIBS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "GSL libs: ${GSL_LIBS}")
  string(REGEX MATCH "^-L[^ ]+" GSL_LIB_PATH ${GSL_LIBS})
  if(NOT GSL_LIB_PATH)
    string(REGEX REPLACE "bin$" "lib" GSL_LIB_PATH ${GSL_CONFIG_DIR})
  else()
    string(REGEX REPLACE "^-L" "" GSL_LIB_PATH ${GSL_LIB_PATH})
  endif()
  unset(GSL_CONFIG CACHE)
endif()


#Optimization and warnings 
set ( CMAKE_CXX_FLAGS "-O3 -Wall ${GSL_CFLAGS}")
include_directories(${BOOST_INC})
 
 
 
 
###### make uninstall ###############
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
 #              "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" @ONLY)
#add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} 
#                 -P ${CMAKE_CURRENT_BINARY_DIR}/#cmake_uninstall.cmake)

 
 
 
 


