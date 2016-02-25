##
#  CMake script
##

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# find deal.II package
FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS
    ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ ../../../ $ENV{DEAL_II_DIR}
  )

# report if deal.II was not bound
IF (NOT ${deal.II_FOUND})
   MESSAGE(FATAL_ERROR
           "\n\n"
	   " *** Could not locate deal.II. *** "
	   "\n\n"
           " *** You may want to either pass the -DDEAL_II_DIR=/path/to/deal.II\n"
           " *** flag to cmake or set an environment variable \"DEAL_II_DIR\"\n"
           " *** that contains this path.")
ENDIF ()

# Set the name of the project and target:
SET(TARGET "transport")

# Declare all source files the target consists of:
SET(TARGET_SRC
  ${TARGET}.cc
  # You can specify additional files here!
  )

# This must be invoked before the PROJECT() call. It sets an uninitialized
# CMAKE_BUILD_TYPE to Debug (or Release if Debug is not available). Also it
# sets the C++ compiler to the one used for compiling deal.II. The variables
# CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_DEBUG, and CMAKE_CXX_FLAGS_RELEASE will
# be initialized to empty strings.
DEAL_II_INITIALIZE_CACHED_VARIABLES()

# specify the project
PROJECT(${TARGET})

# get lists of header and source files
FILE(GLOB_RECURSE my_SOURCES ${PROJECT_SOURCE_DIR}/*.cc)
FILE(GLOB_RECURSE my_HEADERS ${PROJECT_SOURCE_DIR}/*.h)

# get list of include directories
SET(my_INCLUDE_DIRS "")
FOREACH(_headerFile ${my_HEADERS})
    GET_FILENAME_COMPONENT(_dir ${_headerFile} PATH)
    LIST(APPEND my_INCLUDE_DIRS ${_dir})
ENDFOREACH()
LIST(REMOVE_DUPLICATES my_INCLUDE_DIRS) # remove duplicates

# get list of source directories
SET(my_SOURCE_DIRS "")
FOREACH(_sourceFile ${my_SOURCES})
    GET_FILENAME_COMPONENT(_dir ${_sourceFile} PATH)
    LIST(APPEND my_SOURCE_DIRS ${_dir})
ENDFOREACH()
LIST(REMOVE_DUPLICATES my_SOURCE_DIRS) # remove duplicates

# specify include directories
INCLUDE_DIRECTORIES(${my_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${my_SOURCE_DIRS})

# specify include directories in other directory
INCLUDE_DIRECTORIES(../conservation_law)

# don't know what this does
DEAL_II_INVOKE_AUTOPILOT()

