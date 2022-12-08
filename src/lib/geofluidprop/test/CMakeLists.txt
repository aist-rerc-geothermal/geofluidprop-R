

include_directories( ../src )
include_directories( ../src/interface )
if (USE_PROST)
  include_directories("${PROST_INCLUDE_DIR}")
endif()

add_subdirectory( interface )
add_subdirectory( model )
add_subdirectory( util )

add_executable( test-all test_main.c ${SOURCES})
target_link_libraries(test-all ${EOSLIBNAME})
#set_target_properties(test-all PROPERTIES FOLDER "test/util")

add_test(
    NAME run_test_all
    COMMAND COMMAND $<TARGET_FILE:test-all>
)
