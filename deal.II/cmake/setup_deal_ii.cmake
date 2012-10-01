#
# Set up deal.II specific definitions and look for available components
#
# This file defines a long list of variables, used throughout the
# configuration to determine paths, locations and names:
#
# General configuration options:
#
#     DEAL_II_FEATURE_AUTODETECTION   **)
#     DEAL_II_ALLOW_BUNDLED           **)
#     DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS **)
#     DEAL_II_COMPONENT_COMPAT_FILES  **)
#     DEAL_II_COMPONENT_DOCUMENTATION **)
#     DEAL_II_COMPONENT_EXAMPLES      **)
#
# General information about deal.II:
#
#     DEAL_II_PACKAGE_NAME            *)
#     DEAL_II_PACKAGE_VERSION
#     DEAL_II_PACKAGE_STRING
#     DEAL_II_PACKAGE_BUGREPORT       *)
#     DEAL_II_PACKAGE_TARNAME         *)
#     DEAL_II_PACKAGE_URL             *)
#     DEAL_II_VERSION_MAJOR
#     DEAL_II_VERSION_MINOR
#
# Information about paths, install locations and names:
#
#     DEAL_II_PROJECT_CONFIG_NAME     *)
#     DEAL_II_BASE_NAME               *)
#     DEAL_II_DEBUG_SUFFIX            *)
#     DEAL_II_RELEASE_SUFFIX          *)
#     DEAL_II_LIBRARY_NAME_DEBUG
#     DEAL_II_LIBRARY_NAME_RELEASE
#
#     DEAL_II_PATH                    *)
#     DEAL_II_CMAKE_MACROS_RELDIR     *)
#     DEAL_II_DOCUMENTATION_RELDIR    *)
#     DEAL_II_EXAMPLES_RELDIR         *)
#     DEAL_II_EXECUTABLE_RELDIR       *)
#     DEAL_II_INCLUDE_RELDIR          *)
#     DEAL_II_LIBRARY_RELDIR          *)
#     DEAL_II_PROJECT_CONFIG_RELDIR   *)
#
#     DEAL_II_INCLUDE_DIRS
#     DEAL_II_LIBRARIES
#     DEAL_II_LIBRARIES_DEBUG
#     DEAL_II_LIBRARIES_RELEASE
#
#     DEAL_II_BUILD_TYPES
#     DEAL_II_WITH_BUNDLED_DIRECTORY
#     DEAL_II_WITH_DOC_DIRECTORY
#
# *)  Values marked with *) can be overwritten by the command line via
#     -D<...>
# **) Cached Option. Can be set via ccmake or on the command line via
#     -D<...>
#


#
# Check whether the doc or bundled folder is available:
#
IF(EXISTS ${CMAKE_SOURCE_DIR}/bundled/CMakeLists.txt)
  SET(DEAL_II_WITH_BUNDLED_DIRECTORY TRUE)
ENDIF()

IF(EXISTS ${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt)
  SET(DEAL_II_WITH_DOC_DIRECTORY TRUE)
ENDIF()


###########################################################################
#                                                                         #
#                     General configuration options:                      #
#                                                                         #
###########################################################################

OPTION(DEAL_II_FEATURE_AUTODETECTION
  "Enables feature autodetection. This will automatically set DEAL_II_WITH_<...> options."
  ON
  )

OPTION(DEAL_II_ALLOW_BUNDLED
  "Allow the use of libraries bundled with the source tarball. (DEAL_II_FORCE_BUNDLED* will overwrite this option.)"
  ON
  )

#
# Build configuration: configuration options regarding compilation and
# installation of the deal.II library
#
OPTION(DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS
  "configure sensible default CFLAGS and CXXFLAGS depending on platform, compiler and build target."
  ON
  )

#
# Component selection: configuration options regarding the setup of
# different components of the deal.II library
#

OPTION(DEAL_II_COMPONENT_COMPAT_FILES
  "Enable installation of the example steps. This adds a COMPONENT \"compat_files\" to the build system."
  ON
  )

OPTION(DEAL_II_COMPONENT_EXAMPLES
  "Enable configuration and installation of the example steps. This adds a COMPONENT \"examples\" to the build system."
  ON
  )

If(DEAL_II_WITH_DOC_DIRECTORY)
  OPTION(DEAL_II_COMPONENT_DOCUMENTATION
    "Enable configuration, build and installation of the documentation. This adds a COMPONENT \"documentation\" to the build system."
    OFF
    )
ENDIF()


###########################################################################
#                                                                         #
#                   General information about deal.II:                    #
#                                                                         #
###########################################################################

SET_IF_EMPTY(DEAL_II_PACKAGE_NAME "deal.II")

SET(DEAL_II_PACKAGE_VERSION ${VERSION})
SET(DEAL_II_PACKAGE_STRING
  "${DEAL_II_PACKAGE_NAME} ${DEAL_II_PACKAGE_VERSION}"
  )

SET_IF_EMPTY(DEAL_II_PACKAGE_BUGREPORT "dealii@dealii.org")
SET_IF_EMPTY(DEAL_II_PACKAGE_TARNAME ${DEAL_II_PACKAGE_NAME}) #TODO
SET_IF_EMPTY(DEAL_II_PACKAGE_URL "http://www.dealii.org") #TODO

STRING(REGEX REPLACE
  "^([0-9]+)\\..*" "\\1" DEAL_II_VERSION_MAJOR "${VERSION}"
  )
STRING(REGEX REPLACE
  "^[0-9]+\\.([0-9]+).*" "\\1" DEAL_II_VERSION_MINOR "${VERSION}"
  )


###########################################################################
#                                                                         #
#          Information about paths, install locations and names:          #
#                                                                         #
###########################################################################

SET(DEAL_II_PROJECT_CONFIG_NAME "${DEAL_II_PACKAGE_NAME}")

SET_IF_EMPTY(DEAL_II_BASE_NAME "deal_II")
SET_IF_EMPTY(DEAL_II_DEBUG_SUFFIX ".g")
SET_IF_EMPTY(DEAL_II_RELEASE_SUFFIX "")

SET(DEAL_II_PATH ${CMAKE_INSTALL_PREFIX})

IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # The good, old directory structure:
  #
  SET_IF_EMPTY(DEAL_II_CMAKE_MACROS_RELDIR "cmake/macros")
  SET_IF_EMPTY(DEAL_II_DOCUMENTATION_RELDIR "doc")
  SET_IF_EMPTY(DEAL_II_EXAMPLES_RELDIR "examples")
  SET_IF_EMPTY(DEAL_II_EXECUTABLE_RELDIR "bin")
  SET_IF_EMPTY(DEAL_II_INCLUDE_RELDIR "include")
  SET_IF_EMPTY(DEAL_II_LIBRARY_RELDIR "lib")
  SET_IF_EMPTY(DEAL_II_PROJECT_CONFIG_RELDIR ".")
ELSE()
  #
  # IF DEAL_II_COMPONENT_COMPAT_FILES is not set, we assume that we have to
  # obey the FSHS...
  #
  SET_IF_EMPTY(DEAL_II_CMAKE_MACROS_RELDIR "share/${DEAL_II_PACKAGE_NAME}/cmake/Macros")
  SET_IF_EMPTY(DEAL_II_DOCUMENTATION_RELDIR "share/doc/${DEAL_II_PACKAGE_NAME}/html")
  SET_IF_EMPTY(DEAL_II_EXAMPLES_RELDIR "share/doc/${DEAL_II_PACKAGE_NAME}/examples")
  SET_IF_EMPTY(DEAL_II_EXECUTABLE_RELDIR "bin")
  SET_IF_EMPTY(DEAL_II_INCLUDE_RELDIR "include")
  SET_IF_EMPTY(DEAL_II_LIBRARY_RELDIR "lib${LIB_SUFFIX}")
  SET_IF_EMPTY(DEAL_II_PROJECT_CONFIG_RELDIR "${DEAL_II_LIBRARY_RELDIR}/cmake/${DEAL_II_PROJECT_CONFIG_NAME}")
ENDIF()

LIST(APPEND DEAL_II_INCLUDE_DIRS
  "${CMAKE_INSTALL_PREFIX}/${DEAL_II_INCLUDE_RELDIR}"
  "${CMAKE_INSTALL_PREFIX}/${DEAL_II_INCLUDE_RELDIR}/deal.II"
  "${CMAKE_INSTALL_PREFIX}/${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled"
  )


#
# A lot of magic to guess the resulting library names:
#
IF(BUILD_SHARED_LIBS)
  SET(DEAL_II_LIBRARY_NAME_RELEASE ${CMAKE_SHARED_LIBRARY_PREFIX}${DEAL_II_BASE_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX} )
  SET(DEAL_II_LIBRARY_NAME_DEBUG ${CMAKE_SHARED_LIBRARY_PREFIX}${DEAL_II_BASE_NAME}${DEAL_II_DEBUG_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX} )
ELSE()
  SET(DEAL_II_LIBRARY_NAME_RELEASE ${CMAKE_STATIC_LIBRARY_PREFIX}${DEAL_II_BASE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX} )
  SET(DEAL_II_LIBRARY_NAME_DEBUG ${CMAKE_STATIC_LIBRARY_PREFIX}${DEAL_II_BASE_NAME}${DEAL_II_DEBUG_SUFFIX}${CMAKE_STATIC_LIBRARY_SUFFIX} )
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  LIST(APPEND DEAL_II_LIBRARIES_DEBUG
    "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}/${DEAL_II_LIBRARY_NAME_DEBUG}"
    )
  LIST(APPEND DEAL_II_LIBRARIES
    debug
    "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}/${DEAL_II_LIBRARY_NAME_DEBUG}"
    )

  LIST(APPEND DEAL_II_BUILD_TYPES "DEBUG")
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Release")
  LIST(APPEND DEAL_II_LIBRARIES_RELEASE
    "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}/${DEAL_II_LIBRARY_NAME_RELEASE}"
    )
  LIST(APPEND DEAL_II_LIBRARIES
    optimized
    "${CMAKE_INSTALL_PREFIX}/${DEAL_II_LIBRARY_RELDIR}/${DEAL_II_LIBRARY_NAME_RELEASE}"
    )

  LIST(APPEND DEAL_II_BUILD_TYPES "RELEASE")
ENDIF()



#
# Cleanup some files used for storing the names of alle object targets that
# will be bundled to the deal.II library. (Right now, i.e. cmake 2.8.8,
# this is the only reliable way to get information in a global scope...)
#
FOREACH(build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER "${build}" build_lowercase)
  FILE(REMOVE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_objects_${build_lowercase}
    )
ENDFOREACH()

