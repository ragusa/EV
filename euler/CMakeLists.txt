##
#  CMake script:
##

# add option to prevent c++11 standard error
SET(GNUXX_FLAGS "-std=c++0x")
ADD_DEFINITIONS(${GNUXX_FLAGS})

# Usually, you will not need to modify anything beyond this point...

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
           " *** You may want to either pass the -DDEAL_II_DIR=/path/to/deal.II flag to cmake \n"
           " *** or set an environment variable \"DEAL_II_DIR\" that contains this path.")
ENDIF ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(myproject)

# serial targets
ADD_EXECUTABLE(burgers burgers.cc)
DEAL_II_SETUP_TARGET(burgers)

ADD_EXECUTABLE(euler euler.cc)
DEAL_II_SETUP_TARGET(euler)

ADD_EXECUTABLE(shallowwater shallowwater.cc)
DEAL_II_SETUP_TARGET(shallowwater)

# parallel targets
ADD_EXECUTABLE(burgers_parallel burgers.cc)
TARGET_COMPILE_DEFINITIONS(burgers_parallel PRIVATE IS_PARALLEL)
DEAL_II_SETUP_TARGET(burgers_parallel)

ADD_EXECUTABLE(euler_parallel euler.cc)
TARGET_COMPILE_DEFINITIONS(euler_parallel PRIVATE IS_PARALLEL)
DEAL_II_SETUP_TARGET(euler_parallel)

ADD_EXECUTABLE(shallowwater_parallel euler.cc)
TARGET_COMPILE_DEFINITIONS(shallowwater_parallel PRIVATE IS_PARALLEL)
DEAL_II_SETUP_TARGET(shallowwater_parallel)

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
