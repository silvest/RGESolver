cmake_minimum_required(VERSION 3.1)

project(Solver VERSION 1.0 DESCRIPTION "C++ library to solve RGEs in the context of the SMEFT")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS "-Wall ${GSL_CFLAGS}  -DBOOST_NO_CXX11_RVALUE_REFERENCES")

include_directories(${BOOST_INC})

file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/RGESolver.cpp")
add_library(Solver ${LIBTYPE} ${srcs})

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
INSTALL(FILES ${headers} DESTINATION include/Solver COMPONENT header)

