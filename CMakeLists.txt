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
set(WITH_OPENMP YES)

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
  if(WITH_OPENMP)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
  endif(WITH_OPENMP)
endif(MSVC)

# This overrides CXX flags for MSVC
if(MSVC)
  if(WIN64)
    set(MSVC_DEFINES "/DWIN64 /D_WINDOWS /Dpopen=_popen /Dpclose=_pclose /D__value=_value /Dfinite=_finite /Dhypot=_hypot /Disatty=_isatty /Dfileno=_fileno /D_CRT_SECURE_NO_WARNINGS /DYY_NO_UNISTD_H /D_USE_MATH_DEFINES /DIMPLEMENT_C99 /wd4275 /wd4251")
  else(WIN64)
    set(MSVC_DEFINES "/DWIN32 /D_WINDOWS /Dpopen=_popen /Dpclose=_pclose /D__value=_value /Dfinite=_finite /Dhypot=_hypot /Disatty=_isatty /Dfileno=_fileno /D_CRT_SECURE_NO_WARNINGS /DYY_NO_UNISTD_H /D_USE_MATH_DEFINES /DIMPLEMENT_C99 /wd4275 /wd4251")
  endif(WIN64)
  set(CMAKE_CXX_FLAGS_DEBUG          "/D_DEBUG /Od /Ob2 /MDd /Zi ${MSVC_DEFINES}")
  set(HERMES_DEBUG_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "/O2 /Ob2 /MD ${MSVC_DEFINES}")
  set(CMAKE_CXX_FLAGS_RELEASE        "/DNDEBUG /O2 /Ob2 /MD ${MSVC_DEFINES}")
  set(HERMES_RELEASE_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "/DNDEBUG /O2 /Ob2 /MD /Zi ${MSVC_DEFINES}")
  
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SAFESEH:NO")
  set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /SAFESEH:NO")
  set (CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} /SAFESEH:NO")
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

#------------------------------------------------------------------------------
# Report.
#~~~~~~~~
message("\nHermes-testing Configuration results")
message("---------------------")
if(WIN32)
  message("Build for Win64: ${WIN64}")
endif(WIN32)

message("\n-------Features-------")
message("Build with OpenMP: ${WITH_OPENMP}")
message("Build with TCMalloc: ${WITH_TC_MALLOC}")
message("Build with BSON: ${WITH_BSON}")
message("Build with MATIO: ${WITH_MATIO}")
if(${WITH_MATIO})
  message(" MATIO with HDF5: ${MATIO_WITH_HDF5}")
endif()
if(${WITH_MPI})
  message("Build with MPI: ${WITH_MPI}")
endif()
message("Build with EXODUSII: ${WITH_EXODUSII}")

message("\n-------Solvers-------")
message("Build with UMFPACK: ${WITH_UMFPACK}")
message("Build with PARALUTION: ${WITH_PARALUTION}")
message("Build with PETSC: ${WITH_PETSC}")
message("Build with MUMPS: ${WITH_MUMPS}")
message("Build with SUPERLU${MT}: ${WITH_SUPERLU}")
message("Build with TRILINOS: ${WITH_TRILINOS}")

message("---------------------")
message("\n")