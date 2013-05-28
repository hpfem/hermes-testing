project(loadSave-test-01-mesh)
add_executable(${PROJECT_NAME} definitions.cpp main.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(loadSave-test-01-mesh ${BIN} a)