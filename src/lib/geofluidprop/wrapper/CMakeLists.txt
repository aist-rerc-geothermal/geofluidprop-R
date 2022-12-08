
IF(POLICY CMP0078)
    cmake_policy(SET CMP0078 NEW)
ENDIF()
IF(POLICY CMP0086)
    cmake_policy(SET CMP0086 NEW)
ENDIF()

if (BUILD_R_LIB)
  find_package(SWIG REQUIRED)
  find_package(R REQUIRED)
  include(${SWIG_USE_FILE})
  include_directories(${R_INCLUDE_DIR})
  set(INTERFACE_FILES swig/geofluidprop.i)
  include_directories("${PROJECT_SOURCE_DIR}/src")
  set(CMAKE_SWIG_OUTDIR "${PROJECT_BINARY_DIR}/lib/r")
  if (MSVC)
    set(CMAKE_SWIG_OUTDIR "${PROJECT_BINARY_DIR}/lib/r/${CMAKE_BUILD_TYPE}")
  endif()
  swig_add_library(geofluidprop_r  TYPE SHARED LANGUAGE R
  #    OUTPUT_DIR "${PROJECT_BINARY_DIR}/lib"
  #    OUTFILE_DIR "${PROJECT_BINARY_DIR}/lib"
      SOURCES ${INTERFACE_FILES} ${SOURCES}
  )
  set_target_properties(geofluidprop_r PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib/r")
  set_target_properties(geofluidprop_r PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib/r")
  set_target_properties(geofluidprop_r PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib/r")
  set_target_properties(geofluidprop_r PROPERTIES OUTPUT_NAME "geofluidprop")
  #set_property(TARGET geofluidprop_r PROPERTY SWIG_USE_LIBRARY_INCLUDE_DIRECTORIES TRUE)
  swig_link_libraries(geofluidprop_r ${EOSLIBNAME})
  if(WIN32)
    # On windows, it seems we need to link R.lib
    if (DEFINED ENV{R_LIB})
      message(STATUS "found R_LIB path in env: $ENV{R_LIB}")
      set(tmp_R_LIB $ENV{R_LIB})
    endif()
    if (DEFINED R_LIB)
      message(STATUS "R_LIB path is provided via cmake arguments. use this path: ${R_LIB}")
      set(tmp_R_LIB ${R_LIB})
    elseif(NOT DEFINED tmp_R_LIB)
      message(ERROR "cannot find R_LIB variable. The variable should be passed to cmake")
    endif()
    swig_link_libraries(geofluidprop_r ${tmp_R_LIB})
  endif()
endif()

if (BUILD_PYTHON_LIB)
  find_package(SWIG REQUIRED)
  find_package(PythonLibs REQUIRED)
  include(${SWIG_USE_FILE})
  include_directories(${PYTHON_INCLUDE_PATH})
  set(INTERFACE_FILES swig/geofluidprop.i)
  set(CMAKE_SWIG_OUTDIR "${PROJECT_BINARY_DIR}/lib/python")
  swig_add_library(geofluidprop_py TYPE SHARED LANGUAGE python
      OUTPUT_DIR "${PROJECT_BINARY_DIR}/lib/python"
      OUTFILE_DIR "${PROJECT_BINARY_DIR}/lib/python"
      SOURCES ${INTERFACE_FILES} ${SOURCES}
  )
  set_target_properties(geofluidprop_py PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib/python")
  swig_link_libraries(geofluidprop_py ${EOSLIBNAME} ${PYTHON_LIBRARIES})
endif()

