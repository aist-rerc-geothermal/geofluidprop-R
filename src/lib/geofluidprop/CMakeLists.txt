#-----------------------------------------------------------
# CMake project setting
#-----------------------------------------------------------
cmake_minimum_required (VERSION 2.8.3)
enable_testing()

project(GEOFLUIDPROP)
set(EOSLIBNAME "geofluidprop")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/scripts/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/third-party/VTK_v8.2.0)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Set CMAKE_BUILD_TYPE = '${default_build_type}' as default")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_CONFIGURATION_TYPES Debug Release)
if(NOT MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
endif()
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#-----------------------------------------------------------
# User options
#-----------------------------------------------------------
option(BUILD_EXAMPLES "Build examples" ON)
option(USE_PROST "Use PROST library" OFF)
option(USE_FREESTEAM "Use freesteam library" OFF)

#-----------------------------------------------------------
# User options (experimental)
#-----------------------------------------------------------
option(USE_LONGDOUBLE "Use long double in some calculations" OFF)
if(CMAKE_COMPILER_IS_GNUCC)
  option(USE_QUAD "Use quad precision in some calculations" OFF)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  option(USE_QUAD "Use quad precision in some calculations" OFF)
else()
  option(USE_QUAD "Use quad precision in some calculations" OFF)
endif()

option(BUILD_SHARED_LIB "Build shared lib" OFF) # for Windows DLL
option(BUILD_R_LIB "Build R lib" OFF) # requires SWIG, R
option(BUILD_PYTHON_LIB "Build python3 lib" OFF) # requires SWIG


#-----------------------------------------------------------
# Compilier setting
#-----------------------------------------------------------
if (BUILD_SHARED_LIB OR BUILD_R_LIB OR BUILD_PYTHON_LIB)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

if(USE_PROST)
  find_package(PROST REQUIRED)
  add_definitions(-DUSE_PROST)
endif()
if(USE_FREESTEAM)
  find_package(FreeSteam REQUIRED)
  add_definitions(-DUSE_FREESTEAM)
endif()
if(USE_QUAD)
  add_definitions(-DUSE_QUAD)
endif()
if(USE_LONGDOUBLE)
  add_definitions(-DUSE_LONGDOUBLE)
endif()

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --std=gnu11")
  # use gnu11 to compile utest
  #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --std=c11")

  # if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  #   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
  # else()
  #   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
  # endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -diag-disable 161")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS_DEBUG "-fp-model precise -prec-div -g")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    #set(CMAKE_C_FLAGS_RELEASE "-O3 -fp-model precise -prec-div -xHost -DNDEBUG")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -DNDEBUG")
    #set(CMAKE_C_FLAGS_RELEASE "-ipo -O3 -fp-model precise -prec-div -xHost -DNDEBUG")
  endif()
elseif(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_compile_options(/wd4068)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --std=gnu11")
endif()


#-----------------------------------------------------------
# Set source code files
#-----------------------------------------------------------
include_directories(third-party/utest)
add_subdirectory( src )
add_subdirectory( test )
add_subdirectory( wrapper )
if (BUILD_EXAMPLES)
  add_subdirectory( examples )
endif()
