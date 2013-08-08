project(test-matrix-export-real)

add_executable(${PROJECT_NAME} main.cpp)

set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

add_test(test-umfpack-matrix_export-1 ${BIN} umfpack 1)
add_test(test-umfpack-matrix_export-2 ${BIN} umfpack 2)
add_test(test-umfpack-matrix_export-3 ${BIN} umfpack 3)

add_test(test-umfpack-matrix_export_block-1 ${BIN} umfpack_block 1)
add_test(test-umfpack-matrix_export_block-2 ${BIN} umfpack_block 2)
add_test(test-umfpack-matrix_export_block-3 ${BIN} umfpack_block 3)

add_test(test-petsc-matrix_export-1 ${BIN} petsc 1)
add_test(test-petsc-matrix_export-2 ${BIN} petsc 2)
add_test(test-petsc-matrix_export-3 ${BIN} petsc 3)

add_test(test-petsc-matrix_export_block-1 ${BIN} petsc_block 1)
add_test(test-petsc-matrix_export_block-2 ${BIN} petsc_block 2)
add_test(test-petsc-matrix_export_block-3 ${BIN} petsc_block 3)

add_test(test-aztecoo-matrix_export-1 ${BIN} aztecoo 1)
add_test(test-aztecoo-matrix_export-2 ${BIN} aztecoo 2)
add_test(test-aztecoo-matrix_export-3 ${BIN} aztecoo 3)

add_test(test-aztecoo-matrix_export_block-1 ${BIN} aztecoo_block 1)
add_test(test-aztecoo-matrix_export_block-2 ${BIN} aztecoo_block 2)
add_test(test-aztecoo-matrix_export_block-3 ${BIN} aztecoo_block 3)

add_test(test-amesos-matrix_export-1 ${BIN} amesos 1)
add_test(test-amesos-matrix_export-2 ${BIN} amesos 2)
add_test(test-amesos-matrix_export-3 ${BIN} amesos 3)

add_test(test-amesos-matrix_export_block-1 ${BIN} amesos_block 1)
add_test(test-amesos-matrix_export_block-2 ${BIN} amesos_block 2)
add_test(test-amesos-matrix_export_block-3 ${BIN} amesos_block 3)

add_test(test-mumps-matrix_export-1 ${BIN} mumps 1)
add_test(test-mumps-matrix_export-2 ${BIN} mumps 2)
add_test(test-mumps-matrix_export-3 ${BIN} mumps 3)

add_test(test-mumps-matrix_export_block-1 ${BIN} mumps_block 1)
add_test(test-mumps-matrix_export_block-2 ${BIN} mumps_block 2)
add_test(test-mumps-matrix_export_block-3 ${BIN} mumps_block 3)
