project(test-matrix-export-real)

add_executable(${PROJECT_NAME} main.cpp)

set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

add_test(test-umfpack-solver-1 ${BIN} umfpack 1)
add_test(test-umfpack-solver-2 ${BIN} umfpack 2)
add_test(test-umfpack-solver-3 ${BIN} umfpack 3)

add_test(test-umfpack-solver-b-1 ${BIN} umfpack-block 1)
add_test(test-umfpack-solver-b-2 ${BIN} umfpack-block 2)
add_test(test-umfpack-solver-b-3 ${BIN} umfpack-block 3)

if(HAVE_AZTECOO)
  add_test(test-aztecoo-solver-1 ${BIN} aztecoo 1)
  add_test(test-aztecoo-solver-2 ${BIN} aztecoo 2)
  add_test(test-aztecoo-solver-3 ${BIN} aztecoo 3)

  add_test(test-aztecoo-solver-b-1 ${BIN} aztecoo-block 1)
  add_test(test-aztecoo-solver-b-2 ${BIN} aztecoo-block 2)
  add_test(test-aztecoo-solver-b-3 ${BIN} aztecoo-block 3)
endif(HAVE_AZTECOO)

if(HAVE_AMESOS)
  add_test(test-amesos-solver-1 ${BIN} amesos 1)
  add_test(test-amesos-solver-2 ${BIN} amesos 2)
  add_test(test-amesos-solver-3 ${BIN} amesos 3)

  add_test(test-amesos-solver-b-1 ${BIN} amesos-block 1)
  add_test(test-amesos-solver-b-2 ${BIN} amesos-block 2)
  add_test(test-amesos-solver-b-3 ${BIN} amesos-block 3)
endif(HAVE_AMESOS)

if(WITH_MUMPS)
add_test(test-mumps-solver-1 ${BIN} mumps 1)
add_test(test-mumps-solver-2 ${BIN} mumps 2)
add_test(test-mumps-solver-3 ${BIN} mumps 3)

add_test(test-mumps-solver-b-1 ${BIN} mumps-block 1)
add_test(test-mumps-solver-b-2 ${BIN} mumps-block 2)
add_test(test-mumps-solver-b-3 ${BIN} mumps-block 3)
endif(WITH_MUMPS)
