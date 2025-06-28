## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Generate the source documentation via doxygen:
#

find_package(Perl)
find_package(Doxygen)

#
# Do we have all necessary dependencies?
#
if(NOT PERL_FOUND)
  message(FATAL_ERROR
    "Could not find a perl installation which is required for building the documentation"
    )
endif()
if(NOT DOXYGEN_FOUND)
  message(FATAL_ERROR
    "Could not find doxygen which is required for building the documentation"
    )
else()
  message(STATUS "Using doxygen version ${DOXYGEN_VERSION}")
endif()


########################################################################
#
# Process the tutorial and code-gallery files into inputs for doxygen
#
########################################################################

add_subdirectory(tutorial)
add_subdirectory(code-gallery)


########################################################################
#
# Set up all of the other input pieces we want to give to doxygen
#
########################################################################

#
# Prepare auxiliary files for doxygen:
#
if(DEAL_II_DOXYGEN_USE_MATHJAX)
  set(_use_mathjax YES)
  if(NOT DEAL_II_DOXYGEN_USE_ONLINE_MATHJAX)
    deal_ii_find_path(MATHJAX_PATH MathJax.js
      PATHS ${MATHJAX_ROOT}
      $ENV{MATHJAX_ROOT}
      "${CMAKE_PREFIX_PATH}/share/javascript/mathjax"
      "${CMAKE_INSTALL_DATADIR}/javascript/mathjax"
      "/usr/share/javascript/mathjax"
      "/usr/share/yelp/mathjax"
      DOC "Root path of MathJax.js"
      )
    if(MATHJAX_PATH MATCHES "MATHJAX_PATH-NOTFOUND")
      message(FATAL_ERROR "MathJax was not found.")
    endif()
  endif()
  if(MATHJAX_PATH AND NOT DEAL_II_DOXYGEN_USE_ONLINE_MATHJAX)
    set(_mathjax_relpath ${MATHJAX_PATH})
    message(STATUS "Using local MathJax: " ${MATHJAX_PATH})
  else()
    set(_mathjax_relpath "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/")
    message(STATUS "Using web-based MathJax: " ${_mathjax_relpath})
  endif()
else()
  set(_use_mathjax NO)
endif()


# Note where the build system will write the .inst files after
# expansion from the .inst.in files. We need this because we let
# doxygen read the .cc files and they #include the .inst files.
#
# The .cc files #include the .inst files with the directory name
# relative to the /source directory under the build directory, so the
# correct include path is <build dir>/source
#
# This variable is used in options.dox.in where ${_inst_in_dir} is
# added to INCLUDE_PATH.
set(_inst_in_dir "${CMAKE_BINARY_DIR}/source/")


message(STATUS "Additional doxygen include path: ${_inst_in_dir}")

# make sure doxygen sees the extra.sty stylefile:
set(_extra_packages "${CMAKE_CURRENT_SOURCE_DIR}/extra")

# Generate the input files for doxygen using template .in files in which
# we have to substitute some CMake variables
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/options.dox.in
  ${CMAKE_CURRENT_BINARY_DIR}/options.dox
  )

#
# Set up DoxygenLayout.xml. Note that the "Modules" tab has been renamed to
# "Topics" in Doxygen 1.9.8 to avoid confusion with C++ modules.
# Unfortunately, this requires us to play a little bit of a configuration
# dance:
#
if(${DOXYGEN_VERSION} VERSION_LESS 1.9.8)
  set(_topics_layout_flag "<tab type=\"modules\" visible=\"yes\" title=\"\" intro=\"\"/>")
else()
  set(_topics_layout_flag "<tab type=\"topics\" visible=\"yes\" title=\"\" intro=\"\"/>")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml
  )

# Figure out the last copyright date of any of the deal.II source
# files. We will use this then to set the copyright date of the
# doxygen-generated HTML files.
set(_last_year "1997")
file(GLOB _source_files
     "${CMAKE_SOURCE_DIR}/source/*/*.cc"
     "${CMAKE_SOURCE_DIR}/include/deal.II/*/*.h"
     "${CMAKE_SOURCE_DIR}/examples/*/*.cc")
foreach (_source ${_source_files})
  file(READ "${_source}" contents)
  string(REGEX MATCH "Copyright .C. +([0-9]+ - )?([0-9]+)" _copyright "${contents}")
  if (NOT "${_copyright}" STREQUAL "")
    string(REGEX MATCH "[0-9]+\$" _year "${_copyright}")

    if (NOT _year STREQUAL ""   AND   _year GREATER _last_year)
      set (_last_year "${_year}")
    endif()
  endif()
endforeach()
set(DEAL_II_PACKAGE_YEAR ${_last_year})
message(STATUS "doxygen output will have copyright date ${_last_year}")


# Put the copyright year into the correct files
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mod_header.pl.in
  ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl
  )

#
# Generate header, footer and style files for doxygen.
#
# A bug in (at least) doxygen 1.8.12 required that these files already exist
# if they are listed in the options.dox file, even though the -w command is
# specifically intended to create them. See
#    https://bugzilla.gnome.org/show_bug.cgi?id=771606
# To work around this, do a 'touch' operation on them first to ensure they're
# there before we call "doxygen -w".
#
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/header.html
         ${CMAKE_CURRENT_BINARY_DIR}/footer.html
  COMMAND ${CMAKE_COMMAND} -E touch header.html
  COMMAND ${CMAKE_COMMAND} -E touch footer.html
  COMMAND ${DOXYGEN_EXECUTABLE} -w html header.html footer.html
                                   style.css options.dox
  COMMAND ${PERL_EXECUTABLE} -pi ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl header.html
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl
  )

#
# Finalize the doxygen configuration:
#

set(_doxygen_input
  ${CMAKE_CURRENT_SOURCE_DIR}/headers/
)

file(GLOB _changelog_files
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  )

list(APPEND _doxygen_input
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_SOURCE_DIR}/source/
  ${CMAKE_BINARY_DIR}/include/
  ${_changelog_files}
  ${CMAKE_BINARY_DIR}/doc/news/
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial/tutorial.h
  )

# Add other directories in which to find images
set(_doxygen_image_path
  ${CMAKE_CURRENT_SOURCE_DIR}/images
  )

file(GLOB _doxygen_depend
  ${CMAKE_CURRENT_SOURCE_DIR}/headers/*.h
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  ${CMAKE_SOURCE_DIR}/include/deal.II/**/*.h
  )

# Specifically list a few files that are generated by
# other cmake targets as dependencies for doxygen, to ensure
# that they are in fact generated. This would not work if they
# are captured via a GLOB since they may not exist (and so not
# be captured via the GLOB) at the time cmake runs
list(APPEND _doxygen_depend
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial/tutorial.h
  )

# Add the changelog (a CMake custom target) as a dependency for doxygen too
list(APPEND _doxygen_depend changelog)

# find all tutorial programs so we can add dependencies as appropriate
file(GLOB _deal_ii_steps
  ${CMAKE_SOURCE_DIR}/examples/step-*
  )
foreach(_step ${_deal_ii_steps})
  get_filename_component(_step "${_step}" NAME)
  list(APPEND _doxygen_depend
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
  list(APPEND _doxygen_input
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
endforeach()

# Also find all code gallery programs (if available) for the same reason.
# The logic here follows the same as in code-gallery/CMakeLists.txt
set_if_empty(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
if (EXISTS ${DEAL_II_CODE_GALLERY_DIRECTORY}/README.md)
  file(GLOB _code_gallery_names
       "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
  string(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}")

  foreach(_step ${_code_gallery_names})
    get_filename_component(_step "${_step}" NAME)

    list(APPEND _doxygen_depend
      ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/${_step}.h
      )
    list(APPEND _doxygen_input
      ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/${_step}.h
      )
  endforeach()
endif()


# always make the doxygen run depend on the code-gallery.h file
# (whether generated from the code gallery or copied from
# no-code-gallery.h; both happen in code-gallery/CMakeLists.txt)
list(APPEND _doxygen_input
  ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/code-gallery.h
)
list(APPEND _doxygen_depend
  ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/code-gallery.h
)


to_string(_doxygen_image_path_string ${_doxygen_image_path})
to_string(_doxygen_input_string ${_doxygen_input})

file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
  "
  INPUT=${_doxygen_input_string}
  IMAGE_PATH=${_doxygen_image_path_string}
  "
  )

# If we use a reasonably modern doxygen version, make sure it is run in parallel
if(NOT (${DOXYGEN_VERSION} VERSION_LESS 1.9))
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
    "
    NUM_PROC_THREADS       = 0
    DOT_NUM_THREADS        = 0
    "
    )
  message(STATUS "Letting doxygen run with multiple threads")
endif()


# doxygen version 1.9.1 can not correctly parse C++20-style 'requires'
# clauses and just drops certain classes and shows documentation for
# DEAL_II_CXX20_REQUIRES instead. 1.9.7 gets this right, but has other
# bugs. It's unclear for versions in between. Err on the safe side, and
# only make them visible to doxygen if we have at least 1.9.8.
if(${DOXYGEN_VERSION} VERSION_LESS 1.9.8)
  message(STATUS "Suppressing 'requires' clause parsing because doxygen is too old")
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
    "
    PREDEFINED += DEAL_II_DOXYGEN_DO_NOT_PARSE_REQUIRES_CLAUSES=1
    PREDEFINED += DEAL_II_CXX20_REQUIRES(x)=
    "
    )
endif()



########################################################################
#
# And, finally, call doxygen:
#
########################################################################

add_custom_command(
  OUTPUT
    ${CMAKE_BINARY_DIR}/doxygen.log
  COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    > ${CMAKE_BINARY_DIR}/doxygen.log 2>&1 # *pssst*
    || ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/doxygen.log ${CMAKE_BINARY_DIR}/doxygen.err
  COMMAND ${PERL_EXECUTABLE}
  ARGS
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/set_canonical_doxygen.pl
  COMMAND ${CMAKE_COMMAND} -E echo "-- Documentation is available at ${CMAKE_CURRENT_BINARY_DIR}/deal.II/index.html"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    tutorial
    code-gallery
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${CMAKE_CURRENT_BINARY_DIR}/header.html
    ${CMAKE_CURRENT_BINARY_DIR}/footer.html
    ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter.pl
    ${_doxygen_depend}
    expand_all_instantiations
  COMMENT "Generating documentation via doxygen."
  VERBATIM
  )
add_custom_target(doxygen ALL
  DEPENDS ${CMAKE_BINARY_DIR}/doxygen.log
  )
add_dependencies(documentation doxygen)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/deal.tag
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
  COMPONENT documentation
  )

install(FILES
  ${CMAKE_SOURCE_DIR}/doc/deal.ico
  ${CMAKE_SOURCE_DIR}/doc/doxygen/custom.js
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen/deal.II
  COMPONENT documentation
  )

install(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
  COMPONENT documentation
  )
