project(01-poisson)
add_executable(${PROJECT_NAME} definitions.cpp main.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

add_test(01-poisson ${BIN})
