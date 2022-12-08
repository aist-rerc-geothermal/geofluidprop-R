#########################################################################
# - Find PROST
# Find the PROST include and library
#
# This module defines
#  PROST_INCLUDE_DIR
#  PROST_LIBRARY
#  PROST_FOUND
#########################################################################

find_path(PROST_INCLUDE_DIR steam4.h
  PATHS
    ENV PROST_DIR
    ${PROST_DIR}
    /usr/local
  PATH_SUFFIXES
    include
)

find_library(PROST_LIBRARY
  NAMES
    steam4
  PATHS
    ENV PROST_DIR
    ${PROST_DIR}
    /usr/local
  PATH_SUFFIXES
    lib
)

mark_as_advanced(
  PROST_INCLUDE_DIR
  PROST_LIBRARY
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PROST
  REQUIRED_VARS
    PROST_INCLUDE_DIR
    PROST_LIBRARY
  )

if(PROST_FOUND AND NOT TARGET PROST::PROST)
  add_library(PROST::PROST UNKNOWN IMPORTED)
  set_target_properties(PROST::PROST PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "C"
    IMPORTED_LOCATION "${PROST_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PROST_INCLUDE_DIR}"
    )
endif()
