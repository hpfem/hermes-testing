project(test-integral-calculator)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)

set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-integral-calculator ${BIN})
