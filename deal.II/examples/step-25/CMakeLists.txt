##
#         An autopilot CMakeLists.txt example for a small project
#                        using the deal.II library
##

#
# Set the name of the project and target:
#
SET(TARGET "step-25")

#
# Declare all source files the target consists of:
#
SET(TARGET_SRC
  ${TARGET}.cc
  # You can specify additional files here!
  )

#
# Specify files that will be deleted when calling 'make runclean' or
# 'make distclean':
#
SET(CLEAN_UP_FILES *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.vtk *.ucd *.d2)


# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# Search for the deal.II library:
FIND_PACKAGE(deal.II 8.0 REQUIRED
  HINTS
    ${DEAL_II_DIRS}
  # Additional hints for possible library locations:
    ../ ../../ ../../../ ../../../../../
  # Also look for hints from environment:
    $ENV{DEAL_II_DIR} $ENV{DEAL_DIR}
  #
  # If the deal.II library cannot be found (because it is not installed at a
  # default location or your project resides at an uncommon place), you
  # can specify additional hints for search paths here, e.g.
  # "$ENV{HOME}/workspace/deal.II"
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()

