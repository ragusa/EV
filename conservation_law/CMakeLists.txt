##
#  CMake script:
##

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS
    ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  #
  # If the deal.II library cannot be found (because it is not installed at a
  # default location or your project resides at an uncommon place), you
  # can specify additional hints for search paths here, e.g.
  # "$ENV{HOME}/workspace/deal.II"
  )

IF (NOT ${deal.II_FOUND})
   MESSAGE(FATAL_ERROR
           "\n\n"
	   " *** Could not locate deal.II. *** "
	   "\n\n"
           " *** You may want to either pass the \n"
           " ***   -DDEAL_II_DIR=/path/to/deal.II\n"
           " *** flag to cmake or set an environment variable \"DEAL_II_DIR\"\n"
           " *** that contains this path.")
ENDIF ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(myproject)

SET(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
CONFIGURE_FILE(
  "${PROJECT_SOURCE_DIR}/include/other/CMakeVars.h.in"
  "${PROJECT_BINARY_DIR}/include/other/CMakeVars.h")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

# serial targets
ADD_EXECUTABLE(transport main.cc)
TARGET_COMPILE_DEFINITIONS(transport PRIVATE TRANSPORT)
DEAL_II_SETUP_TARGET(transport)

#ADD_EXECUTABLE(burgers main.cc)
#TARGET_COMPILE_DEFINITIONS(burgers PRIVATE BURGERS)
#DEAL_II_SETUP_TARGET(burgers)
#
#ADD_EXECUTABLE(euler main.cc)
#TARGET_COMPILE_DEFINITIONS(euler PRIVATE EULER)
#DEAL_II_SETUP_TARGET(euler)

ADD_EXECUTABLE(shallowwater main.cc)
TARGET_COMPILE_DEFINITIONS(shallowwater PRIVATE SHALLOWWATER)
DEAL_II_SETUP_TARGET(shallowwater)

# parallel targets
#ADD_EXECUTABLE(transport_parallel main.cc)
#TARGET_COMPILE_DEFINITIONS(transport_parallel PRIVATE TRANSPORT)
#TARGET_COMPILE_DEFINITIONS(transport_parallel PRIVATE IS_PARALLEL)
#DEAL_II_SETUP_TARGET(transport_parallel)
#
#ADD_EXECUTABLE(burgers_parallel main.cc)
#TARGET_COMPILE_DEFINITIONS(burgers_parallel PRIVATE BURGERS)
#TARGET_COMPILE_DEFINITIONS(burgers_parallel PRIVATE IS_PARALLEL)
#DEAL_II_SETUP_TARGET(burgers_parallel)
#
#ADD_EXECUTABLE(euler_parallel main.cc)
#TARGET_COMPILE_DEFINITIONS(euler_parallel PRIVATE EULER)
#TARGET_COMPILE_DEFINITIONS(euler_parallel PRIVATE IS_PARALLEL)
#DEAL_II_SETUP_TARGET(euler_parallel)
#
#ADD_EXECUTABLE(shallowwater_parallel main.cc)
#TARGET_COMPILE_DEFINITIONS(shallowwater_parallel PRIVATE SHALLOWWATER)
#TARGET_COMPILE_DEFINITIONS(shallowwater_parallel PRIVATE IS_PARALLEL)
#DEAL_II_SETUP_TARGET(shallowwater_parallel)

# add debug target
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

# add release target
ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

# testing
ENABLE_TESTING()
ADD_SUBDIRECTORY(tests/group_fe_values)
ADD_SUBDIRECTORY(tests/quadrature)
ADD_SUBDIRECTORY(tests/transport_galerkin_solution)
ADD_SUBDIRECTORY(tests/transport_low)
ADD_SUBDIRECTORY(tests/transport_ev)
ADD_SUBDIRECTORY(tests/transport_fct_none)
ADD_SUBDIRECTORY(tests/transport_fct_full)
ADD_SUBDIRECTORY(tests/shallowwater_galerkin_solution)
ADD_SUBDIRECTORY(tests/shallowwater_constant_solution)
ADD_SUBDIRECTORY(tests/shallowwater_low_solution)
ADD_SUBDIRECTORY(tests/shallowwater_ev_solution)
ADD_SUBDIRECTORY(tests/shallowwater_low_lakeatrest)

