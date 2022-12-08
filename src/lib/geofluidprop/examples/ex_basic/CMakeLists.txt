
include_directories(${PROJECT_SOURCE_DIR}/src/interface)
add_executable( ex_basic ex_basic.c)
target_link_libraries(ex_basic ${EOSLIBNAME})
