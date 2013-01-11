##
#  CMake script for the step-42 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "step-42")

# Declare all source files the target consists of:
SET(TARGET_SRC
  ${TARGET}.cc
  # You can specify additional files here!
  )

# A custom command line to run the program
SET(TARGET_RUN mpirun -np 2 ${TARGET})

# Define the output that should be cleaned:
SET(CLEAN_UP_FILES *.vtu *.pvtu *.visit)

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS
    ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ ../../../ ../../../../../ $ENV{DEAL_II_DIR}
  #
  # If the deal.II library cannot be found (because it is not installed at a
  # default location or your project resides at an uncommon place), you
  # can specify additional hints for search paths here, e.g.
  # "$ENV{HOME}/workspace/deal.II"
  )

#
# Are all dependencies fullfilled?
#
IF (NOT ${deal.II_FOUND})
   MESSAGE(FATAL_ERROR
           "\n\n"
	   " *** Could not locate deal.II. *** "
	   "\n\n"
           " *** You may want to either pass the -DDEAL_II_DIR=/path/to/deal.II flag to cmake \n"
           " *** or set an environment variable \"DEAL_II_DIR\" that contains this path.")
ENDIF ()

IF( NOT DEAL_II_WITH_MPI OR
    NOT DEAL_II_WITH_P4EST OR
    NOT DEAL_II_WITH_TRILINOS )
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_MPI = ON
    DEAL_II_WITH_P4EST = ON
    DEAL_II_WITH_TRILINOS = ON
which is required for this tutorial step."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
