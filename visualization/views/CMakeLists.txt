project(visualization-test-views)
add_executable(${PROJECT_NAME} definitions.cpp main.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME} ${BIN})

project(memory-leaks-visualization)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

IF(NOT MSVC AND NOT(NO_VALGRIND_TESTS))
  add_test(memory-leaks-visualization "sudo" "valgrind" "--leak-check=full" "--track-origins=yes" "--log-file=/var/www/valgrind-logfiles/${PROJECT_NAME}-logfile" ${BIN})
ENDIF()

project(threads-visualization)
set(CMAKE_CXX_FLAGS "-DTHREAD_TESTING ${CMAKE_CXX_FLAGS}")
add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")
set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

IF(NOT MSVC AND NOT(NO_VALGRIND_TESTS))
  add_test(threads-assembling "sudo" "valgrind" "--tool=drd" "--log-file=/var/www/valgrind-logfiles/${PROJECT_NAME}-logfile" ${BIN})
ENDIF()
