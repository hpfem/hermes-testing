project(test-matrix-import-export-complex)

add_executable(${PROJECT_NAME} main.cpp)

set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

add_test(test-petsc-matrix_export-cplx-1 ${BIN} petsc)
add_test(test-petsc-matrix_export-cplx_block-1 ${BIN} petsc_block)

add_test(test-umfpack-matrix_export-cplx-1 ${BIN} umfpack)
add_test(test-umfpack-matrix_export-cplx_block-1 ${BIN} umfpack_block)

add_test(test-aztecoo-matrix_export-cplx-1 ${BIN} aztecoo)
add_test(test-aztecoo-matrix_export-cplx_block-1 ${BIN} aztecoo_block)

add_test(test-mumps-matrix_export-cplx-1 ${BIN} mumps)
add_test(test-mumps-matrix_export-cplx_block-1 ${BIN} mumps_block)
