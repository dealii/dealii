## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Generate the source documentation via doxygen:
#

FIND_PACKAGE(Perl)
FIND_PACKAGE(Doxygen)

#
# Do we have all necessary dependencies?
#
IF(NOT PERL_FOUND)
  MESSAGE(FATAL_ERROR
    "Could not find a perl installation which is required for building the documentation"
    )
ENDIF()
IF(NOT DOXYGEN_FOUND)
  MESSAGE(FATAL_ERROR
    "Could not find doxygen which is required for building the documentation"
    )
ENDIF()

########################################################################
#
# Process the tutorial and code-gallery files into inputs for doxygen
#
########################################################################

ADD_SUBDIRECTORY(tutorial)
ADD_SUBDIRECTORY(code-gallery)


########################################################################
#
# Set up all of the other input pieces we want to give to doxygen
#
########################################################################

#
# Prepare auxiliary files for doxygen:
#
IF(DEAL_II_DOXYGEN_USE_MATHJAX)
  SET(_use_mathjax YES)
  IF(NOT DEAL_II_DOXYGEN_USE_ONLINE_MATHJAX)
    DEAL_II_FIND_PATH(MATHJAX_PATH MathJax.js
      PATHS ${MATHJAX_ROOT}
      $ENV{MATHJAX_ROOT}
      "${CMAKE_PREFIX_PATH}/share/javascript/mathjax"
      "${CMAKE_INSTALL_DATADIR}/javascript/mathjax"
      "/usr/share/javascript/mathjax"
      "/usr/share/yelp/mathjax"
      DOC "Root path of MathJax.js"
      )
    IF(MATHJAX_PATH MATCHES "MATHJAX_PATH-NOTFOUND")
      MESSAGE(FATAL_ERROR "MathJax was not found.")
    ENDIF()
  ENDIF()
  IF(MATHJAX_PATH AND NOT DEAL_II_DOXYGEN_USE_ONLINE_MATHJAX)
    SET(_mathjax_relpath ${MATHJAX_PATH})
    MESSAGE(STATUS "Using local MathJax: " ${MATHJAX_PATH})
  ELSE()
    SET(_mathjax_relpath "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/")
    MESSAGE(STATUS "Using web-based MathJax: " ${_mathjax_relpath})
  ENDIF()
ELSE()
  SET(_use_mathjax NO)
ENDIF()


# Collect where the build system will write the .inst files after
# expansion from the .inst.in files. We need this because we let
# doxygen read the .cc files and they #include the .inst files.
#
# To get at this, find all of the .inst.in files, strip all but the
# directory name under $CMAKE_SOURCE_DIR/source, and then prepend
# the build directory to it
FILE(GLOB _inst_in_files
     "${CMAKE_SOURCE_DIR}/source/*/*.inst.in")
STRING(REPLACE "${CMAKE_SOURCE_DIR}/source/" "" _inst_in_dirs_x "${_inst_in_files}")
STRING(REGEX REPLACE "/[a-zA-Z0-9_\\.]*.inst.in" "/" _inst_in_dirs_y "${_inst_in_dirs_x}")
LIST(REMOVE_DUPLICATES _inst_in_dirs_y)
SET(_inst_in_dirs "")
FOREACH (_dir ${_inst_in_dirs_y})
  SET(_inst_in_dirs "${_inst_in_dirs} ${CMAKE_BINARY_DIR}/source/${_dir}")
ENDFOREACH()

MESSAGE(STATUS "translated this into")
MESSAGE(STATUS "found ${_inst_in_dirs}")

# make sure doxygen sees the extra.sty stylefile:
SET(_extra_packages "${CMAKE_CURRENT_SOURCE_DIR}/extra")

# Generate the input files for doxygen using template .in files in which
# we have to substitute some CMake variables
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/options.dox.in
  ${CMAKE_CURRENT_BINARY_DIR}/options.dox
  )

# Figure out the last copyright date of any of the deal.II source
# files. We will use this then to set the copyright date of the
# doxygen-generated HTML files.
SET(_last_year "1997")
FILE(GLOB _source_files
     "${CMAKE_SOURCE_DIR}/source/*/*.cc"
     "${CMAKE_SOURCE_DIR}/include/deal.II/*/*.h"
     "${CMAKE_SOURCE_DIR}/examples/*/*.cc")
FOREACH (_source ${_source_files})
  FILE(READ "${_source}" contents)
  STRING(REGEX MATCH "Copyright .C. +([0-9]+ - )?([0-9]+)" _copyright "${contents}")
  IF (NOT _copright STREQUAL "")
    STRING(REGEX MATCH "[0-9]+\$" _year "${_copyright}")

    IF (NOT _year STREQUAL ""   AND   _year GREATER _last_year)
      SET (_last_year "${_year}")
    ENDIF()
  ENDIF()
ENDFOREACH()
SET(DEAL_II_PACKAGE_YEAR ${_last_year})
MESSAGE(STATUS "doxygen output will have copyright date ${_last_year}")


# Put the copyright year into the correct files
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mod_header.pl.in
  ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mod_footer.pl.in
  ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_footer.pl
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
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/header.html
         ${CMAKE_CURRENT_BINARY_DIR}/footer.html
  COMMAND ${CMAKE_COMMAND} -E touch header.html
  COMMAND ${CMAKE_COMMAND} -E touch footer.html
  COMMAND ${DOXYGEN_EXECUTABLE} -w html header.html footer.html style.css options.dox
  COMMAND ${PERL_EXECUTABLE} -pi ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl header.html
  COMMAND ${PERL_EXECUTABLE} -pi ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_footer.pl footer.html
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_header.pl
    ${CMAKE_CURRENT_BINARY_DIR}/scripts/mod_footer.pl
  )

#
# Finalize the doxygen configuration:
#

SET(_doxygen_input
  ${CMAKE_CURRENT_SOURCE_DIR}/headers/
)

file(GLOB _changelog_files
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  )

LIST(APPEND _doxygen_input
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_SOURCE_DIR}/source/
  ${CMAKE_BINARY_DIR}/include/
  ${_changelog_files}
  ${CMAKE_BINARY_DIR}/doc/news/
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial/tutorial.h
  )

# Add other directories in which to find images
SET(_doxygen_image_path
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
LIST(APPEND _doxygen_depend
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial/tutorial.h
  )

# Add the changelog (a CMake custom target) as a dependency for doxygen too
LIST(APPEND _doxygen_depend changelog)

# find all tutorial programs so we can add dependencies as appropriate
FILE(GLOB _deal_ii_steps
  ${CMAKE_SOURCE_DIR}/examples/step-*
  )
FOREACH(_step ${_deal_ii_steps})
  GET_FILENAME_COMPONENT(_step "${_step}" NAME)
  LIST(APPEND _doxygen_depend
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
  LIST(APPEND _doxygen_input
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial/${_step}.h
    )
ENDFOREACH()

# Also find all code gallery programs (if available) for the same reason.
# The logic here follows the same as in code-gallery/CMakeLists.txt
SET_IF_EMPTY(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
IF (EXISTS ${DEAL_II_CODE_GALLERY_DIRECTORY}/README.md)
  FILE(GLOB _code_gallery_names
       "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
  STRING(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}")

  FOREACH(_step ${_code_gallery_names})
    GET_FILENAME_COMPONENT(_step "${_step}" NAME)

    LIST(APPEND _doxygen_depend
      ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/${_step}.h
      )
    LIST(APPEND _doxygen_input
      ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/${_step}.h
      )
  ENDFOREACH()
ENDIF()


# always make the doxygen run depend on the code-gallery.h file
# (whether generated from the code gallery or copied from
# no-code-gallery.h; both happen in code-gallery/CMakeLists.txt)
LIST(APPEND _doxygen_input
  ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/code-gallery.h
)
LIST(APPEND _doxygen_depend
  ${CMAKE_CURRENT_BINARY_DIR}/code-gallery/code-gallery.h
)


TO_STRING(_doxygen_image_path_string ${_doxygen_image_path})
TO_STRING(_doxygen_input_string ${_doxygen_input})

FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
  "
  INPUT=${_doxygen_input_string}
  IMAGE_PATH=${_doxygen_image_path_string}
  "
  )

########################################################################
#
# And, finally, call doxygen:
#
########################################################################

ADD_CUSTOM_COMMAND(
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
    ${CMAKE_CURRENT_SOURCE_DIR}/DoxygenLayout.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/filter
    ${_doxygen_depend}
    expand_all_instantiations
  COMMENT "Generating documentation via doxygen."
  VERBATIM
  )
ADD_CUSTOM_TARGET(doxygen ALL
  DEPENDS ${CMAKE_BINARY_DIR}/doxygen.log
  )
ADD_DEPENDENCIES(documentation doxygen)

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/deal.tag
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
  COMPONENT documentation
  )

INSTALL(FILES
  ${CMAKE_SOURCE_DIR}/doc/deal.ico
  ${CMAKE_SOURCE_DIR}/doc/doxygen/custom.js
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen/deal.II
  COMPONENT documentation
  )

INSTALL(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
  COMPONENT documentation
  )
