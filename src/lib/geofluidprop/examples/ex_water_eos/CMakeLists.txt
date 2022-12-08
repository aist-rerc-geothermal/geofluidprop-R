
include_directories(${PROJECT_SOURCE_DIR}/src/interface)
add_executable( ex_water_eos ex_water_eos.c)
target_link_libraries(ex_water_eos ${EOSLIBNAME})
