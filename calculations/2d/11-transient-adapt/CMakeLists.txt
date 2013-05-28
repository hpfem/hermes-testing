project(11-transient-adapt)
add_executable(${PROJECT_NAME} definitions.cpp main.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

add_test(${PROJECT_NAME} ${BIN})