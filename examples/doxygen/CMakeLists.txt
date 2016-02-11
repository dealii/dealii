##
#  CMake script for compiling the doxygen examples
##

# Usually, you will not need to modify anything beyond this point...

FILE(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.5.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(doxygen_examples)

FOREACH(ccfile ${sources})
  GET_FILENAME_COMPONENT(file ${ccfile} NAME_WE)
  ADD_EXECUTABLE(${file} ${ccfile})
  DEAL_II_SETUP_TARGET(${file})
ENDFOREACH()
