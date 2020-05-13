## ---------------------------------------------------------------------
##
## Copyright (C) 2015 - 2020 by the deal.II Authors
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
# A target for the preparation of all the stuff happening in here...
#

ADD_CUSTOM_TARGET(code-gallery)


#
# Check whether someone has either specified a global variable that
# points to a code gallery directory, or whether it has been put into
# a top-level code-gallery/ directory alongside the tutorial/
# directory
#
SET_IF_EMPTY(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
IF (EXISTS ${DEAL_II_CODE_GALLERY_DIRECTORY}/README.md)

  MESSAGE(STATUS "Setting up code gallery documentation from ${DEAL_II_CODE_GALLERY_DIRECTORY}.")

  # Collect the names of all code gallery projects. To
  # do so, find all 'author' files, then strip the last two
  # levels of these paths.
  #
  # For unclear reasons, the glob returns these files as
  # "/a/b/c/name//doc//author", so make sure we eat the
  # double slashes in the second step
  FILE(GLOB _code_gallery_names
       "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
  STRING(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}")

  FOREACH(_step ${_code_gallery_names})
    GET_FILENAME_COMPONENT(_step "${_step}" NAME)
    LIST(APPEND _code_gallery_names_sans_dir "${_step}")
  ENDFOREACH()

  # Describe how to build code-gallery.h. Make sure we properly
  # track dependencies on the input files, by assuming that the PERL
  # script is going to read all of the files in the doc/ subdirectories
  FILE(GLOB _code_gallery_h_deps
       "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/*")
  STRING(REPLACE "//" "/" _code_gallery_h_deps "${_code_gallery_h_deps}")

  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
    COMMAND ${PERL_EXECUTABLE}
    ARGS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/code-gallery.pl
      ${CMAKE_CURRENT_SOURCE_DIR}/code-gallery.h.in
      ${DEAL_II_CODE_GALLERY_DIRECTORY}
      ${_code_gallery_names_sans_dir}
      > ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
    DEPENDS
      ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/code-gallery.pl
      ${CMAKE_CURRENT_SOURCE_DIR}/code-gallery.h.in
      ${_code_gallery_h_deps}
    )
  ADD_CUSTOM_TARGET(build_code-gallery_h
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
      ${_code_gallery_h_deps}
    COMMENT
      "Building code-gallery.h")
  ADD_DEPENDENCIES(code-gallery build_code-gallery_h)


  # Now set up targets for each of the code gallery programs
  FOREACH(_step ${_code_gallery_names})
    GET_FILENAME_COMPONENT(_step "${_step}" NAME)
    MESSAGE(STATUS "  Setting up ${_step}")

    # Get all source files so we can let the perl script work on
    # them and so we properly describe the dependencies. exclude
    # meta-files necessary to describe each code gallery project
    FILE(GLOB_RECURSE _src_files
         ${DEAL_II_CODE_GALLERY_DIRECTORY}/${_step}/*)
    STRING(REPLACE "${DEAL_II_CODE_GALLERY_DIRECTORY}/${_step}/" "" _relative_src_files
           "${_src_files}")
    LIST(REMOVE_ITEM _relative_src_files doc/author)
    LIST(REMOVE_ITEM _relative_src_files doc/tooltip)
    LIST(REMOVE_ITEM _relative_src_files doc/dependencies)
    LIST(REMOVE_ITEM _relative_src_files doc/builds-on)


    # Also remove files that were created by running 'cmake' in the
    # code-gallery directory. These greatly confuse doxygen when we
    # run them through doxygen.
    #
    # We would really like to also exclude executables, but there is
    # no option to cmake's 'IF' command to test for that :-(
    FOREACH(_file ${_relative_src_files})
      IF(IS_DIRECTORY ${_file}
         OR
         ${_file} MATCHES "CMakeFiles/.*"
         OR
         ${_file} MATCHES "CMakeCache.txt")
        LIST(REMOVE_ITEM _relative_src_files ${_file})
      ENDIF()
    ENDFOREACH()

    MESSAGE(STATUS "    Setting up ${_relative_src_files}")

    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
      COMMAND ${PERL_EXECUTABLE}
      ARGS
        ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/make_gallery.pl
        ${CMAKE_SOURCE_DIR}
        ${_step}
        ${DEAL_II_CODE_GALLERY_DIRECTORY}/${_step}
        ${_relative_src_files}
        > ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
      WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS
        ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/make_gallery.pl
        ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2doxygen
        ${_src_files}
      )

    # Copy files of interest (non-metadata) to the build directory
    # so we can link to them, and schedule them for installation
    FILE(COPY ${DEAL_II_CODE_GALLERY_DIRECTORY}/${_step}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
      PATTERN REGEX "doc/tooltip|doc/dependencies|doc/builds-on" EXCLUDE
      )
    INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_step}
      DESTINATION ${DEAL_II_EXAMPLES_RELDIR}/code-gallery
      COMPONENT examples
      )

    # Create a target for this program and add it to the top-level
    # target of this directory
    ADD_CUSTOM_TARGET(code-gallery_${_step}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
      COMMENT "Building doxygen input file for code gallery program <${_step}>"
      )
    ADD_DEPENDENCIES(code-gallery code-gallery_${_step})


  ENDFOREACH()

ELSE()

  # no copy of the code gallery is available. say so. but also
  # install a file that creates a doxygen page we can link to
  # nonetheless, so we don't get bad doxygen references
  MESSAGE(STATUS "Setting up code gallery documentation.")
  MESSAGE(STATUS "  Skipping as no code gallery exists in ${DEAL_II_CODE_GALLERY_DIRECTORY}.")

  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/no-code-gallery.h ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/no-code-gallery.h
	)

  # Make the custom command for code-gallery.h visible to the parent CMakeLists.txt by attaching to the code-gallery
  # custom target:
  ADD_CUSTOM_TARGET(build_code-gallery_h
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
    COMMENT
      "Building code-gallery.h")
  ADD_DEPENDENCIES(code-gallery build_code-gallery_h)

ENDIF()

