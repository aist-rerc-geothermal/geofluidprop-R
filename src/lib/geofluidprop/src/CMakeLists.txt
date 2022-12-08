
include_directories( . util )

add_subdirectory( util )
add_subdirectory( spline )
add_subdirectory( sbtl )
add_subdirectory( model )
add_subdirectory( interface )

add_library(${EOSLIBNAME} ${SOURCES} )
set_target_properties(${EOSLIBNAME} PROPERTIES OUTPUT_NAME ${EOSLIBNAME})

if(CMAKE_COMPILER_IS_GNUCC)
    target_link_libraries(${EOSLIBNAME} m)
endif()
if (USE_QUAD)
  target_link_libraries(${EOSLIBNAME} quadmath)
endif()

if (USE_PROST)
  include_directories("${PROST_INCLUDE_DIR}")
  target_link_libraries(${EOSLIBNAME} "${PROST_LIBRARY}")
endif()

if(USE_FREESTEAM)
  include_directories("${FREESTEAM_INCLUDE_DIR}")
  target_link_libraries(${EOSLIBNAME} "${FREESTEAM_LIBRARY}")
endif()

set_target_properties(${EOSLIBNAME} PROPERTIES LINKER_LANGUAGE C)

set_target_properties(${EOSLIBNAME} PROPERTIES PUBLIC_HEADER
        "${EOSLIBNAME}/eos.h;${EOSLIBNAME}/eos_args.h;${EOSLIBNAME}/eos_type.h")
install(TARGETS ${EOSLIBNAME}
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
)

