
include_directories( . )

add_executable( time-rho time_rho_ph.c )
target_link_libraries(time-rho ${EOSLIBNAME})
set_target_properties(time-rho PROPERTIES FOLDER "test/interface")

set(TMP_SOURCES
    test-eos.c)

list(TRANSFORM TMP_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(SOURCES
    ${SOURCES}
    ${TMP_SOURCES}
    PARENT_SCOPE)
