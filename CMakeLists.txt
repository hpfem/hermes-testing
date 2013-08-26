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

# This has to be the same as in the library's CMake.
set(H2D_WITH_GLUT YES)

set(WITH_TC_MALLOC NO)
set(WITH_MATIO NO)
set(WITH_BSON NO)
set(WITH_PARALUTION YES)
set(WITH_SUPERLU NO)

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

if(${WITH_TC_MALLOC})
  find_package(TCMALLOC REQUIRED)
  include_directories(${TCMALLOC_INCLUDE_DIR})
endif()

if(${WITH_PARALUTION})
  find_package(PARALUTION REQUIRED)
  include_directories(${PARALUTION_INCLUDE_DIR})
endif()

if(${WITH_MATIO})
  find_package(MATIO REQUIRED)
  include_directories(${MATIO_INCLUDE_DIR})
endif()

if(${WITH_BSON})
  find_package(BSON REQUIRED)
  include_directories(${BSON_INCLUDE_DIR})
endif()
      
if(WITH_SUPERLU)
  find_package(SUPERLU REQUIRED)
  include_directories(${SUPERLU_INCLUDE_DIR})
  if(SUPERLU_MT)
    add_definitions(-DSLU_MT)
  endif(SUPERLU_MT)
endif(WITH_SUPERLU)
set(HERMES2D_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes2d)
set(HERMES_COMMON_INCLUDE_PATH ${HERMES_INCLUDE_PATH}/hermes_common)

find_package(HERMES_COMMON REQUIRED)

include_directories(${HERMES_COMMON_INCLUDE_PATH})
include_directories(${HERMES2D_INCLUDE_PATH})
include_directories(${DEP_INCLUDE_PATHS})

enable_testing()

set(TESTING_CORE_LIBRARY "testing-core")
add_subdirectory(testing-core)
if(H2D_WITH_GLUT)
  add_subdirectory(visualization)
endif()
add_subdirectory(load-save)
add_subdirectory(inner-functionality)
add_subdirectory(calculations)