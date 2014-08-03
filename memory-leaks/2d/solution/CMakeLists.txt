project(memory-leaks-solution)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

IF(NOT MSVC AND WITH_VALGRIND_TESTS)
  add_test(memory-leaks-solution "sudo" "valgrind" "--leak-check=full" "--track-origins=yes" "--log-file=/var/www/valgrind-logfiles/${PROJECT_NAME}-logfile" ${BIN})
ENDIF()