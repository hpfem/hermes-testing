project(hermes_testing)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake)
include(CommonTargetProperties)

# For Win64
if(${CMAKE_CL_64})
  set(WIN64 YES)
else(${CMAKE_CL_64})
  set(WIN64 NO)
endif(${CMAKE_CL_64})

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
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /openmp")
  if(64_BIT)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DWIN64")
      string (REPLACE "/DWIN32" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  else(64_BIT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DWIN32")
  endif(64_BIT)
else(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif(MSVC)

include_directories("${PTHREAD_ROOT}/include")

find_package(UMFPACK REQUIRED)
include_directories(${UMFPACK_INCLUDE_DIRS})

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

set(TESTING_CORE_LIBRARY "testing-core")
add_subdirectory(testing-core)
add_subdirectory(visualization)
add_subdirectory(load-save)
add_subdirectory(inner-functionality)
add_subdirectory(calculations)