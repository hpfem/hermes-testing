project(loop-solvers-adapt)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")
set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
add_test("test-${PROJECT_NAME}" ${BIN} a)