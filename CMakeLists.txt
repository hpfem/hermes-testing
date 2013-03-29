project(hermes_testing)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
      include(CommonTargetProperties)

# Where to look for the static libraries.
set(HERMES_DIRECTORY /usr/local/lib)
set(HERMES_INCLUDE_PATH /usr/local/include)
set(DEP_INCLUDE_PATHS /usr/local/include)

# Allow to override the default values in CMake.vars:
include(CMake.vars OPTIONAL)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		# Disable all warnings and turn on only important ones:
		set(CMAKE_CXX_FLAGS "-w ${CMAKE_CXX_FLAGS}")
		set(CMAKE_CXX_FLAGS "-Wuninitialized -Wvla -Wsign-compare ${CMAKE_CXX_FLAGS}")

		set(RELEASE_FLAGS "-DNDEBUG -O3")
		set(DEBUG_FLAGS   "-g")
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

# Enabling multiprocessor build on MSVC
if(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
else(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif(MSVC)

include_directories("${PTHREAD_ROOT}/include")
find_package(PTHREAD REQUIRED)
find_package(TCMALLOC REQUIRED)
include_directories(${TCMALLOC_INCLUDE_DIR})

set(HERMES2D_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes2d)
set(HERMES_COMMON_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes_common)

find_package(HERMES_COMMON REQUIRED)

include_directories(${HERMES_COMMON_INCLUDE_PATH})
include_directories(${HERMES2D_INCLUDE_PATH})
include_directories(${DEP_INCLUDE_PATHS})

enable_testing()

add_subdirectory(memory-leaks)
add_subdirectory(performance)
add_subdirectory(visualization)
add_subdirectory(load-save)
add_subdirectory(inner-functionality)
add_subdirectory(calculations)

