#########################################################################
# - Find freesteam
# Find the freesteam include and library
#
# This module defines
#  FREESTEAM_INCLUDE_DIR
#  FREESTEAM_LIBRARY
#  FREESTEAM_FOUND
#########################################################################

find_path(FREESTEAM_INCLUDE_DIR steam.h
  PATHS
    ENV FREESTEAM_DIR
    ${FREESTEAM_DIR}
    /usr/local
  PATH_SUFFIXES
    include
    include/freesteam
)

find_library(FREESTEAM_LIBRARY
  NAMES
    freesteam
  PATHS
    ENV FREESTEAM_DIR
    ${FREESTEAM_DIR}
    /usr/local
  PATH_SUFFIXES
    lib
)

mark_as_advanced(
  FREESTEAM_INCLUDE_DIR
  FREESTEAM_LIBRARY
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FREESTEAM
  REQUIRED_VARS
    FREESTEAM_INCLUDE_DIR
    FREESTEAM_LIBRARY
  )

if(FREESTEAM_FOUND AND NOT TARGET FREESTEAM::FREESTEAM)
  add_library(FREESTEAM::FREESTEAM UNKNOWN IMPORTED)
  set_target_properties(FREESTEAM::FREESTEAM PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "C"
    IMPORTED_LOCATION "${FREESTEAM_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FREESTEAM_INCLUDE_DIR}"
    )
endif()
