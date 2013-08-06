project(test-solvers-complex)

add_executable(${PROJECT_NAME} main.cpp)

set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

add_test(test-petsc-solver-cplx-1 ${BIN} petsc)
add_test(test-petsc-solver-cplx-b-1 ${BIN} petsc-block)

add_test(test-umfpack-solver-cplx-1 ${BIN} umfpack)
add_test(test-umfpack-solver-cplx-b-1 ${BIN} umfpack-block)

add_test(test-aztecoo-solver-cplx-1 ${BIN} aztecoo)
add_test(test-aztecoo-solver-cplx-b-1 ${BIN} aztecoo-block)

add_test(test-mumps-solver-cplx-1 ${BIN} mumps)
add_test(test-mumps-solver-cplx-b-1 ${BIN} mumps-block)
