## ------------------------------------------------------------------------
##
## Copyright (C) 2012 - 2023 by the deal.II Authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------



# Collect all of the directory names for the tutorial programs
file(GLOB _deal_ii_steps
  ${CMAKE_SOURCE_DIR}/examples/step-*
  )

# Also collect the names of all code gallery projects. To
# do so, find all 'author' files, then strip the last two
# levels of these paths.
#
# For unclear reasons, the glob returns these files as
# "/a/b/c/name//doc//author", so make sure we eat the
# double slashes in the second step
set_if_empty(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
file(GLOB _code_gallery_names
     "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
string(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}")

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

add_custom_target(tutorial)

#
# Describe how to build tutorial.h.
#
# First collect the various files this all depends on. This
# is simpler for the tutorials since the directory names
# are structured. For the code gallery, we have already
# determined the list of names of the relevant subdirectories
# (which may be empty), so we can create the list of files
# the following target depends on without a GLOB and instead
# using just a FOREACH loop.
#
# For tutorial programs, the steps.pl script will read
# tooltip, kind, and builds-on files. For the code gallery
# programs, the 'kind' is obvious and we only need the other
# two categories.
#

file(GLOB _deal_ii_steps_tooltip
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/tooltip
  )
foreach(_gallery ${_code_gallery_names})
  list(APPEND _deal_ii_code_gallery_tooltip
       ${_gallery}/doc/tooltip)
endforeach()
file(GLOB _deal_ii_steps_kind
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/kind
  )
file(GLOB _deal_ii_steps_buildson
  ${CMAKE_SOURCE_DIR}/examples/step-*/doc/builds-on
  )
foreach(_gallery ${_code_gallery_names})
  list(APPEND _deal_ii_code_gallery_buildson
       ${_gallery}/doc/builds-on)
endforeach()


add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
  COMMAND ${PERL_EXECUTABLE}
  ARGS
    ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/steps.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.h.in
    ${_deal_ii_steps}
    ${_code_gallery_names}
    > ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
  DEPENDS
    ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/steps.pl
    ${CMAKE_CURRENT_SOURCE_DIR}/tutorial.h.in
    ${_deal_ii_steps_tooltip}
    ${_deal_ii_steps_kind}
    ${_deal_ii_steps_buildson}
    ${_deal_ii_code_gallery_tooltip}
    ${_deal_ii_code_gallery_buildson}
  )
add_custom_target(build_tutorial_h
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tutorial.h
    COMMENT
      "Building tutorial.h")
add_dependencies(tutorial build_tutorial_h)


#
# Prepare the steps for documentation generation
#

foreach(_step ${_deal_ii_steps})
  get_filename_component(_step "${_step}" NAME)

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_step}.cc
    COMMAND ${PERL_EXECUTABLE}
    ARGS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2plain.pl
      < ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.cc
      > ${CMAKE_CURRENT_BINARY_DIR}/${_step}.cc
    DEPENDS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2plain.pl
      ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.cc
    VERBATIM
    )

  add_custom_command(
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
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/intro2toc.pl
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/create_anchors.pl
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2doxygen.pl
      ${CMAKE_SOURCE_DIR}/examples/${_step}/${_step}.cc
      ${CMAKE_SOURCE_DIR}/examples/${_step}/doc/intro.dox
      ${CMAKE_SOURCE_DIR}/examples/${_step}/doc/results.dox
    )

  add_custom_target(tutorial_${_step}
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
      ${CMAKE_CURRENT_BINARY_DIR}/${_step}.cc
    COMMENT
      "Building doxygen input file for tutorial program <${_step}>"
    )
  add_dependencies(tutorial tutorial_${_step})
endforeach()
