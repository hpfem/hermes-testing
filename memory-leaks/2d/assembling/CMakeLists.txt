project(memory-leaks-assembling)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

IF(NOT MSVC AND NOT(NO_VALGRIND_TESTS))
  add_test(memory-leaks-assembling "sudo" "valgrind" "--leak-check=full" "--track-origins=yes" "--log-file=/var/www/valgrind-logfiles/${PROJECT_NAME}-logfile" ${BIN})
ENDIF()

project(threads-assembling)
set(CMAKE_CXX_FLAGS "-DTHREAD_TESTING ${CMAKE_CXX_FLAGS}")
add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")
set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

IF(NOT MSVC AND NOT(NO_VALGRIND_TESTS))
  add_test(threads-assembling "sudo" "valgrind" "--tool=drd" "--log-file=/var/www/valgrind-logfiles/${PROJECT_NAME}-logfile" ${BIN})
ENDIF()