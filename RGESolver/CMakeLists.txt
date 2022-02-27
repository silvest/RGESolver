#--------------------------------------------------------------------
# List of available command-line options:
#
#   -DLOCAL_INSTALL:BOOL=<ON or OFF>
#    If set to on, install the library in the folder 
#    BUILD_DIR/install
#    
#   -DCMAKE_INSTALL_PREFIX:PATH=<RGESolver installation directory>
#    This variable cannot be modified 
#    when -DLOCAL INSTALL ALL=ON is set.
#   
#   -DDEBUG_MODE:BOOL=<ON or OFF>
#    Prints on-screen cache and normal variables
#
#   -DBOOST_INCLUDE_DIR:PATH=<include path>/boost/
#    Path to boost headers
#
#   -DGSL_CONFIG_DIR:PATH=<gsl-config directory>
#    Path to gsl-config
#
#--------------------------------------------------------------------

cmake_minimum_required(VERSION 3.1)
set(LIBTYPE STATIC)


cmake_policy(SET CMP0048 NEW)

project(RGESolver VERSION 1.0 DESCRIPTION "C++ library to solve RGEs in the context of the SMEFT")
#add_library(RGESolver RGESolver.cc)

#set(LIBTYPE OBJECT)



# default boolean flags
option(LOCAL_INSTALL "Local installation" OFF)
option(DEBUG_MODE "Debug mode" OFF)



#Local installation option
if(LOCAL_INSTALL)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/install
      CACHE PATH "RGESolver installation directory" FORCE)
endif()


#User provided RGESolver installation path
if(NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX})
   set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}
       CACHE PATH "RGESolver installation directory" FORCE)
endif()
string(REGEX REPLACE "^~" "$ENV{HOME}" TMP_DIR ${CMAKE_INSTALL_PREFIX})
set(CMAKE_INSTALL_PREFIX ${TMP_DIR} CACHE PATH "RGESolver installation directory" FORCE)
unset(TMP_DIR)

#Printing the installation directory
message(STATUS "Installation directory : ${CMAKE_INSTALL_PREFIX}$")


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
  endif()#CXX_FLAGS
  unset(GSL_CONFIG CACHE)
endif()


#Compiler flags: optimization and warnings 
set ( CMAKE_CXX_FLAGS "-O3 -Wall ${GSL_CFLAGS}")
include_directories(${BOOST_INC})

#Searches if "-std=c++11" flag is supported. If yes, 
#is added to the compilation 
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_FLAG)
if (HAS_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    #set (CMAKE_CXX_STANDARD 11)
    #set (CXX_STANDARD_REQUIRED ON)
else()
    message(FATAL_ERROR "Required at least c++11") 
endif()

# Fix behavior of CMAKE_CXX_STANDARD when targeting macOS.
if (POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)
endif ()	
	

#Printing the compiler flags
message(STATUS "Compiler flags : ${CMAKE_CXX_FLAGS}")
 
 
 
######### make install ############### 


INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/RGESolver.h 
        DESTINATION include COMPONENT header)

 configure_file("${CMAKE_CURRENT_SOURCE_DIR}/RGESolver.pc.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/RGESolver.pc" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/rgesolver-config.in" 
                 "${CMAKE_CURRENT_BINARY_DIR}/rgesolver-config" @ONLY)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/rgesolver-config ${CMAKE_CURRENT_BINARY_DIR}/RGESolver.pc 
          DESTINATION bin COMPONENT rgesolver-config)

 
##########  DEBUG  ##########
if(DEBUG_MODE)
  message("\nCache variables:")
  message("  CMAKE_CXX_COMPILER:FILEPATH = ${CMAKE_CXX_COMPILER}")
  message("  CMAKE_INSTALL_PREFIX:PATH = ${CMAKE_INSTALL_PREFIX}")
  message("  DEBUG_MODE:BOOL = ${DEBUG_MODE}")
  message("  LOCAL_INSTALL:BOOL = ${LOCAL_INSTALL}")
  message("  BOOST_INCLUDE_DIR:PATH = ${BOOST_INCLUDE_DIR}")
  message("  GSL_CONFIG_DIR:PATH = ${GSL_CONFIG_DIR}")
  message("\nNormal variables:")
  #message("  LIBTYPE = ${LIBTYPE}")
  message("  BOOST_INC = ${BOOST_INC}")
  message("  GSL_VERSION = ${GSL_VERSION}")
  message("  GSL_CFLAGS = ${GSL_CFLAGS}")
  message("  GSL_LIBS = ${GSL_LIBS}")
  message("  GSL_LIB_PATH = ${GSL_LIB_PATH}")
  message("")

endif(DEBUG_MODE)
 
#############################

set(PROJECTLIST gslpp) 
 
# sub-directories

foreach(_project ${PROJECTLIST})
   add_subdirectory(${_project})
endforeach(_project)

#RGESolver library libRGESolver.a
if(${LIBTYPE} STREQUAL "OBJECT") 
  # for CMake v2.8.8 or later

  foreach(_project ${PROJECTLIST})
    set(OBJECTLIST ${OBJECTLIST} $<TARGET_OBJECTS:${_project}>)
  endforeach(_project)
  add_library(RGESolver STATIC ${OBJECTLIST})
else()
  # for older CMake
  set(CMAKE_CXX_FLAGS "-Wall ${GSL_CFLAGS}")
  include_directories(${BOOST_INC})
  #set(SOURCES RGESolver.cpp)

  foreach(_project ${PROJECTLIST})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${_project}/src)
    file(GLOB SOURCES ${_project}/src/*.cpp)

    set(SRCLIST ${SRCLIST} ${SOURCES})
  endforeach(_project)
  add_library(RGESolver STATIC RGESolver.cc ${SRCLIST})
endif()


INSTALL(TARGETS RGESolver DESTINATION lib COMPONENT library)

# RGESolver header file RGESolver.h 
list(REVERSE PROJECTLIST)
foreach(_project ${PROJECTLIST})
  file(GLOB HEADERS ${_project}/src/*.h)
  foreach(_header ${HEADERS})
    get_filename_component(_tmp ${_header} NAME)
  endforeach(_header)
endforeach(_project)

#INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/RGESolver.h 
 #       DESTINATION include#/RGESolver
  #       COMPONENT header)
         
#message(STATUS "{CMAKE_CURRENT_BINARY_DIR}: ${CMAKE_CURRENT_BINARY_DIR}}$") 
# message(STATUS "{CMAKE_CURRENT_SOURCE_DIR}: ${CMAKE_CURRENT_SOURCE_DIR}$") 
 
###### make uninstall ###############
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
             "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} 
                -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

 
 
 
