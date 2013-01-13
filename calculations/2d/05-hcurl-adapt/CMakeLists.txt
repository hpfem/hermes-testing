project(05-hcurl-adapt)
add_executable(${PROJECT_NAME} main.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

add_test(05-hcurl-adapt ${BIN})
