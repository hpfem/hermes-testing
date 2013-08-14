project(test-solvers-real)

add_executable(${PROJECT_NAME} main.cpp)

set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

add_test(test-superlu-solver-1 ${BIN} superlu 1)
add_test(test-superlu-solver-2 ${BIN} superlu 2)
add_test(test-superlu-solver-3 ${BIN} superlu 3)

add_test(test-superlu-solver-block-1 ${BIN} superlu-block 1)
add_test(test-superlu-solver-block-2 ${BIN} superlu-block 2)
add_test(test-superlu-solver-block-3 ${BIN} superlu-block 3)

add_test(test-paralution-solver-1 ${BIN} paralution 1)
add_test(test-paralution-solver-2 ${BIN} paralution 2)
add_test(test-paralution-solver-3 ${BIN} paralution 3)

add_test(test-paralution-solver-block-1 ${BIN} paralution-block 1)
add_test(test-paralution-solver-block-2 ${BIN} paralution-block 2)
add_test(test-paralution-solver-block-3 ${BIN} paralution-block 3)

add_test(test-umfpack-solver-1 ${BIN} umfpack 1)
add_test(test-umfpack-solver-2 ${BIN} umfpack 2)
add_test(test-umfpack-solver-3 ${BIN} umfpack 3)

add_test(test-umfpack-solver-block-1 ${BIN} umfpack-block 1)
add_test(test-umfpack-solver-block-2 ${BIN} umfpack-block 2)
add_test(test-umfpack-solver-block-3 ${BIN} umfpack-block 3)

add_test(test-petsc-solver-1 ${BIN} petsc 1)
add_test(test-petsc-solver-2 ${BIN} petsc 2)
add_test(test-petsc-solver-3 ${BIN} petsc 3)

add_test(test-petsc-solver-block-1 ${BIN} petsc-block 1)
add_test(test-petsc-solver-block-2 ${BIN} petsc-block 2)
add_test(test-petsc-solver-block-3 ${BIN} petsc-block 3)

add_test(test-aztecoo-solver-1 ${BIN} aztecoo 1)
add_test(test-aztecoo-solver-2 ${BIN} aztecoo 2)
add_test(test-aztecoo-solver-3 ${BIN} aztecoo 3)

add_test(test-aztecoo-solver-block-1 ${BIN} aztecoo-block 1)
add_test(test-aztecoo-solver-block-2 ${BIN} aztecoo-block 2)
add_test(test-aztecoo-solver-block-3 ${BIN} aztecoo-block 3)

add_test(test-amesos-solver-1 ${BIN} amesos 1)
add_test(test-amesos-solver-2 ${BIN} amesos 2)
add_test(test-amesos-solver-3 ${BIN} amesos 3)

add_test(test-amesos-solver-block-1 ${BIN} amesos-block 1)
add_test(test-amesos-solver-block-2 ${BIN} amesos-block 2)
add_test(test-amesos-solver-block-3 ${BIN} amesos-block 3)

add_test(test-mumps-solver-1 ${BIN} mumps 1)
add_test(test-mumps-solver-2 ${BIN} mumps 2)
add_test(test-mumps-solver-3 ${BIN} mumps 3)

add_test(test-mumps-solver-block-1 ${BIN} mumps-block 1)
add_test(test-mumps-solver-block-2 ${BIN} mumps-block 2)
add_test(test-mumps-solver-block-3 ${BIN} mumps-block 3)
