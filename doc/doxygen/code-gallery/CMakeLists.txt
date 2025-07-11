## ------------------------------------------------------------------------
##
## Copyright (C) 2015 - 2024 by the deal.II Authors
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
# A target for the preparation of all the stuff happening in here...
#

add_custom_target(code-gallery)


#
# Check whether someone has either specified a global variable that
# points to a code gallery directory, or whether it has been put into
# a top-level code-gallery/ directory alongside the tutorial/
# directory
#
set_if_empty(DEAL_II_CODE_GALLERY_DIRECTORY ${CMAKE_SOURCE_DIR}/code-gallery)
if (EXISTS ${DEAL_II_CODE_GALLERY_DIRECTORY}/README.md)

  message(STATUS "Setting up code gallery documentation from ${DEAL_II_CODE_GALLERY_DIRECTORY}.")

  # Collect the names of all code gallery projects. To
  # do so, find all 'author' files, then strip the last two
  # levels of these paths.
  #
  # For unclear reasons, the glob returns these files as
  # "/a/b/c/name//doc//author", so make sure we eat the
  # double slashes in the second step
  file(GLOB _code_gallery_names
       "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/author")
  string(REGEX REPLACE "/+doc/+author" "" _code_gallery_names "${_code_gallery_names}")

  foreach(_step ${_code_gallery_names})
    get_filename_component(_step "${_step}" NAME)
    list(APPEND _code_gallery_names_sans_dir "${_step}")
  endforeach()

  # Describe how to build code-gallery.h. Make sure we properly
  # track dependencies on the input files, by assuming that the PERL
  # script is going to read all of the files in the doc/ subdirectories
  file(GLOB _code_gallery_h_deps
       "${DEAL_II_CODE_GALLERY_DIRECTORY}/*/doc/*")
  string(REPLACE "//" "/" _code_gallery_h_deps "${_code_gallery_h_deps}")

  add_custom_command(
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
  add_custom_target(build_code-gallery_h
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
      ${_code_gallery_h_deps}
    COMMENT
      "Building code-gallery.h")
  add_dependencies(code-gallery build_code-gallery_h)


  # Now set up targets for each of the code gallery programs
  foreach(_step ${_code_gallery_names})
    get_filename_component(_step "${_step}" NAME)
    message(STATUS "  Setting up ${_step}")

    # Get all source files so we can let the perl script work on
    # them and so we properly describe the dependencies. exclude
    # meta-files necessary to describe each code gallery project
    file(GLOB_RECURSE _src_files
         ${DEAL_II_CODE_GALLERY_DIRECTORY}/${_step}/*)
    string(REPLACE "${DEAL_II_CODE_GALLERY_DIRECTORY}/${_step}/" "" _relative_src_files
           "${_src_files}")
    list(REMOVE_ITEM _relative_src_files doc/author)
    list(REMOVE_ITEM _relative_src_files doc/builds-on)
    list(REMOVE_ITEM _relative_src_files doc/dependencies)
    list(REMOVE_ITEM _relative_src_files doc/entry-name)
    list(REMOVE_ITEM _relative_src_files doc/tooltip)


    # Also remove files that were created by running 'cmake' in the
    # code-gallery directory. These greatly confuse doxygen when we
    # run them through doxygen.
    #
    # We would really like to also exclude executables, but there is
    # no option to cmake's 'IF' command to test for that :-(
    foreach(_file ${_relative_src_files})
      if(IS_DIRECTORY ${_file}
         OR
         ${_file} MATCHES "CMakeFiles/.*"
         OR
         ${_file} MATCHES "CMakeCache.txt")
        list(REMOVE_ITEM _relative_src_files ${_file})
      endif()
    endforeach()

    message(STATUS "    Source files: ${_relative_src_files}")

    add_custom_command(
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
        ${CMAKE_SOURCE_DIR}/doc/doxygen/scripts/program2doxygen.pl
        ${_src_files}
      )

    # Copy files of interest (non-metadata) to the build directory
    # so we can link to them, and schedule them for installation
    file(COPY ${DEAL_II_CODE_GALLERY_DIRECTORY}/${_step}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
      PATTERN REGEX "doc/tooltip|doc/dependencies|doc/builds-on" EXCLUDE
      )
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_step}
      DESTINATION ${DEAL_II_EXAMPLES_RELDIR}/code-gallery
      COMPONENT examples
      )

    # Create a target for this program and add it to the top-level
    # target of this directory
    add_custom_target(code-gallery_${_step}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_step}.h
      COMMENT "Building doxygen input file for code gallery program <${_step}>"
      )
    add_dependencies(code-gallery code-gallery_${_step})


  endforeach()

  # Now also install the directory in which the script puts copies of everything
  # that isn't explicitly annotated and sent through doxygen.
  install(DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
    DESTINATION ${DEAL_II_DOCHTML_RELDIR}/doxygen
    COMPONENT documentation
    PATTERN CMakeFiles EXCLUDE
    PATTERN cmake_install.cmake EXCLUDE
    )


else()

  # no copy of the code gallery is available. say so. but also
  # install a file that creates a doxygen page we can link to
  # nonetheless, so we don't get bad doxygen references
  message(STATUS "Setting up code gallery documentation.")
  message(STATUS "  Skipping as no code gallery exists in ${DEAL_II_CODE_GALLERY_DIRECTORY}.")

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/no-code-gallery.h ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/no-code-gallery.h
	)

  # Make the custom command for code-gallery.h visible to the parent CMakeLists.txt by attaching to the code-gallery
  # custom target:
  add_custom_target(build_code-gallery_h
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/code-gallery.h
    COMMENT
      "Building code-gallery.h")
  add_dependencies(code-gallery build_code-gallery_h)

endif()
