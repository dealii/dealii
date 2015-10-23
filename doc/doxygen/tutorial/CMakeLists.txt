## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II Authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------


#
# Define target for the tutorial. It depends on the
# file tutorial.h built via the next target below, as well
# as the various files we create from the tutorial 
# directories below that. These dependencies are added
# below the respective targets.
#
# This file uses the DEAL_II_STEPS variable set in
# ../CMakeLists.txt.
#

ADD_CUSTOM_TARGET(tutorial)

#
# Describe how to build tutorial.h:
#

file(GLOB DEAL_II_STEPS_BUILDSON
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/builds-on
  )
file(GLOB DEAL_II_STEPS_KIND
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/kind
  )
file(GLOB DEAL_II_STEPS_TOOLTIP
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/tooltip
  )

ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
  COMMAND ${PERL_EXECUTABLE}
  ARGS
    ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/steps.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.h.in
    ${DEAL_II_STEPS}
    > ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
  DEPENDS
    ${DEAL_II_STEPS}
    ${DEAL_II_STEPS_BUILDSON}
    ${DEAL_II_STEPS_KIND}
    ${DEAL_II_STEPS_TOOLTIP}
    ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.h.in
  )
ADD_CUSTOM_TARGET(build_tutorial_h
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h)
ADD_DEPENDENCIES(tutorial build_tutorial_h)


#
# Prepare the steps for documentation generation
#

FOREACH(_step ${DEAL_II_STEPS})
  GET_FILENAME_COMPONENT(_step "${_step}" NAME)

  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_step}.cc
    COMMAND ${PERL_EXECUTABLE}
    ARGS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2plain
      < ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.cc
      > ${CMAKE_CURRENT_BINARY_DIR}/${_step}.cc
    DEPENDS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2plain
      ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.cc
    VERBATIM
    )

  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
    COMMAND ${PERL_EXECUTABLE}
    ARGS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/make_step.pl
      ${_step} ${CMAKE_SOURCE_DIR}
      > ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
    WORKING_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/make_step.pl
      ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.cc
      ${CMAKE_SOURCE_DIR}/examples/${_step}/doc/intro.dox
      ${CMAKE_SOURCE_DIR}/examples/${_step}/doc/results.dox
    )

  ADD_CUSTOM_TARGET(tutorial_${_step}
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
      ${CMAKE_CURRENT_BINARY_DIR}/${_step}.cc
    )
  ADD_DEPENDENCIES(tutorial tutorial_${_step})
ENDFOREACH()
