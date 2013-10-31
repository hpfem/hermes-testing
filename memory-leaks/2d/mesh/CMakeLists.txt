project(memory-leaks-mesh)
add_executable(${PROJECT_NAME} main.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

IF(NOT MSVC)
  add_test(memory-leaks-mesh "valgrind" "--leak-check=full" "--log-file=valgrind-logfile" ${BIN})
ENDIF(NOT MSVC)