##
#  CMake script for the step-43 tutorial program:
##

# Set the name of the project and target:
set(TARGET "step-43")

# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#  file(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#  file(GLOB_RECURSE TARGET_INC  "include/*.h")
#  set(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.
set(TARGET_SRC
  ${TARGET}.cc
  )

# Usually, you will not need to modify anything beyond this point...

cmake_minimum_required(VERSION 3.13.4)

find_package(deal.II 9.7.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()

#
# Are all dependencies fulfilled?
#
if(NOT DEAL_II_WITH_TRILINOS)
  message(FATAL_ERROR "
Error! This tutorial requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_TRILINOS = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options:
    DEAL_II_WITH_TRILINOS = ${DEAL_II_WITH_TRILINOS}
This conflicts with the requirements."
    )
endif()

deal_ii_initialize_cached_variables()
project(${TARGET})
deal_ii_invoke_autopilot()
