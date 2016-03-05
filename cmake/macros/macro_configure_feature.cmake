## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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
# This macro is used for the feature configuration in deal.II
#
# Usage:
#     CONFIGURE_FEATURE(feature)
#
#
# This macro uses the following optional variables and macros:
#
# FEATURE_${feature}_DEPENDS    (a variable)
#    a variable which contains an optional list of other features
#    this feature depends on (and which have to be enabled for this feature
#    to work.)
#    Features must be given with short name, i.e. without DEAL_II_WITH_
#
# FEATURE_${feature}_after      (a variable)
#    a variable which contains an optional list of other features
#    that have to be configured prior to this feature
#    Features must be given with short name, i.e. without DEAL_II_WITH_
#
# FEATURE_${feature}_HAVE_BUNDLED   (a variable)
#    which should either be set to TRUE if all necessary libraries of the
#    features comes bundled with deal.II and hence can be supported
#    without external dependencies, or unset.
#
# FEATURE_${feature}_CONFIGURE_BUNDLED()   (a macro)
#    which should setup all necessary configuration for the feature with
#    bundled source dependencies. If something goes wrong this macro must
#    issue a FATAL_ERROR.
#
# FEATURE_${feature}_FIND_EXTERNAL(var)   (a macro)
#    which should set var to TRUE if all dependencies for the feature are
#    fulfilled. In this case all necessary variables for
#    FEATURE_${feature}_CONFIGURE_EXTERNAL must be set. Otherwise
#    var should remain unset.
#    If not defined, FIND_PACKAGE(${feature}) is called.
#
# FEATURE_${feature}_CONFIGURE_EXTERNAL()   (macro)
#    which should setup all necessary configuration for the feature with
#    external dependencies.
#
# FEATURE_${feature}_ERROR_MESSAGE()  (macro)
#    which should print a meaningful error message (with FATAL_ERROR) for
#    the case that no usable library was found.
#    If not defined, a suitable default error message will be printed.
#


########################################################################
#                                                                      #
#                            Helper Macros:                            #
#                                                                      #
########################################################################

#
# Some black magic to have substitution in command names:
#
MACRO(RUN_COMMAND _the_command)
  FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/macro_configure_feature.tmp"
    "${_the_command}")
  INCLUDE("${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/macro_configure_feature.tmp")
ENDMACRO()


#
# A small macro to set the DEAL_II_WITH_${_feature} variables:
#
MACRO(SET_CACHED_OPTION _str _value)
  STRING(TOLOWER "${_str}" _str_lower)
  SET(DEAL_II_WITH_${_str}
    ${_value}
    CACHE BOOL
    "Build deal.II with support for ${_str_lower}."
    FORCE)
ENDMACRO()


#
# A small macro to post a default error message:
#
MACRO(FEATURE_ERROR_MESSAGE _feature)
  STRING(TOLOWER ${_feature} _feature_lowercase)

  IF(DEFINED ${_feature}_DIR)
    SET(_hint_snippet "
    $ ${_feature}_DIR=\"...\" cmake <...>
    $ cmake -D${_feature}_DIR=\"...\" <...>
or set the relevant variables by hand in ccmake."
      )
  ELSE()
    SET(_hint_snippet
      " or set the relevant variables by hand in ccmake."
      )
  ENDIF()

  IF(FEATURE_${_feature}_HAVE_BUNDLED)
    SET(_bundled_snippet
      "\nAlternatively you may choose to compile the bundled library of "
      "${_feature_lowercase} by setting DEAL_II_ALLOW_BUNDLED=on or "
      "DEAL_II_FORCE_BUNDLED_${_feature}=on.\n"
      )
  ELSE()
    SET(_bundled_snippet "\n")
  ENDIF()

  MESSAGE(FATAL_ERROR "\n"
    "Could not find the ${_feature_lowercase} library!\n"
    ${${_feature}_ADDITIONAL_ERROR_STRING}
    "Please ensure that a suitable ${_feature_lowercase} library is installed on your computer.\n"
    "If the library is not at a default location, either provide some hints "
    "for autodetection,${_hint_snippet}${_bundled_snippet}"
    )
ENDMACRO()


#
# Default macro for finding an external library:
#
MACRO(FEATURE_FIND_EXTERNAL _feature _var)
  FIND_PACKAGE(${_feature})
  IF(${_feature}_FOUND)
    SET(${_var} TRUE)
  ENDIF()
ENDMACRO()


########################################################################
#                                                                      #
#                          CONFIGURE_FEATURE:                          #
#                                                                      #
########################################################################

MACRO(CONFIGURE_FEATURE _feature)

  #
  # Register the feature in the DEAL_II_FEATURES list
  #
  LIST(APPEND DEAL_II_FEATURES ${_feature})

  #
  # This script is arcane black magic. But at least for the better good: We
  # don't have to copy the configuration logic to every single
  # configure_<feature>.cmake script...
  #

  #
  # Check for correct include order of the configure_*.cmake files:
  # If feature B explicitly states to come after feature A, or if feature B
  # depends on feature A, configure_A.cmake has to be included before
  # configure_B.cmake:
  #
  FOREACH(_dependency
      ${FEATURE_${_feature}_AFTER}
      ${FEATURE_${_feature}_DEPENDS}
      )
    IF(NOT FEATURE_${_dependency}_PROCESSED)
      MESSAGE(FATAL_ERROR "\n"
        "Internal build system error: The configuration of "
        "DEAL_II_WITH_${_feature} depends on "
        "DEAL_II_WITH_${_dependency}, but CONFIGURE_FEATURE(${_feature}) "
        "was called before CONFIGURE_FEATURE(${_dependency}).\n\n"
        )
    ENDIF()
  ENDFOREACH()

  #
  # Obey the user overrides:
  #
  IF( (NOT DEAL_II_ALLOW_AUTODETECTION) AND
      (NOT DEFINED DEAL_II_WITH_${_feature}) )
    PURGE_FEATURE(${_feature})
    SET_CACHED_OPTION(${_feature} OFF)
  ENDIF()

  #
  # Only try to configure ${_feature} if we have to, i.e.
  # DEAL_II_WITH_${_feature} is set to true or not set at all.
  #
  IF((NOT DEFINED DEAL_II_WITH_${_feature}) OR
     DEAL_II_WITH_${_feature})

    #
    # Are all dependencies fulfilled?
    #
    SET(_dependencies_ok TRUE)
    FOREACH(_dependency ${FEATURE_${_feature}_DEPENDS})
      IF(NOT DEAL_II_WITH_${_dependency})
        IF(DEAL_II_WITH_${_feature})
          MESSAGE(FATAL_ERROR "\n"
            "DEAL_II_WITH_${_feature} has unmet configuration requirements: "
            "DEAL_II_WITH_${_dependency} has to be set to \"ON\".\n\n"
            )
        ELSE()
          MESSAGE(STATUS
            "DEAL_II_WITH_${_feature} has unmet configuration requirements: "
            "DEAL_II_WITH_${_dependency} has to be set to \"ON\"."
            )
          PURGE_FEATURE(${_feature})
          SET_CACHED_OPTION(${_feature} OFF)
        ENDIF()
        SET(_dependencies_ok FALSE)
      ENDIF()
    ENDFOREACH()

    IF(_dependencies_ok)
      IF(DEAL_II_FORCE_BUNDLED_${_feature})
        #
        # First case: DEAL_II_FORCE_BUNDLED_${_feature} is defined:
        #

        PURGE_FEATURE(${_feature})

        IF(FEATURE_${_feature}_HAVE_BUNDLED)
          RUN_COMMAND("FEATURE_${_feature}_CONFIGURE_BUNDLED()")
          MESSAGE(STATUS "DEAL_II_WITH_${_feature} successfully set up with bundled packages.")
          SET(FEATURE_${_feature}_BUNDLED_CONFIGURED TRUE)
          SET_CACHED_OPTION(${_feature} ON)
        ELSE()
          MESSAGE(FATAL_ERROR "\n"
            "Internal build system error: DEAL_II_FORCE_BUNDLED_${_feature} "
            "defined, but FEATURE_${_feature}_HAVE_BUNDLED not present.\n"
            )
        ENDIF()

      ELSE(DEAL_II_FORCE_BUNDLED_${_feature})
        #
        # Second case: We are allowed to search for an external library
        #
        IF(COMMAND FEATURE_${_feature}_FIND_EXTERNAL)
          RUN_COMMAND("FEATURE_${_feature}_FIND_EXTERNAL(FEATURE_${_feature}_EXTERNAL_FOUND)")
        ELSE()
          FEATURE_FIND_EXTERNAL(${_feature} FEATURE_${_feature}_EXTERNAL_FOUND)
        ENDIF()

        IF(FEATURE_${_feature}_EXTERNAL_FOUND)
          IF(COMMAND FEATURE_${_feature}_CONFIGURE_EXTERNAL)
            RUN_COMMAND("FEATURE_${_feature}_CONFIGURE_EXTERNAL()")
          ENDIF()

          MESSAGE(STATUS "DEAL_II_WITH_${_feature} successfully set up with external dependencies.")
          SET(FEATURE_${_feature}_EXTERNAL_CONFIGURED TRUE)
          SET_CACHED_OPTION(${_feature} ON)

        ELSE(FEATURE_${_feature}_EXTERNAL_FOUND)

          PURGE_FEATURE(${_feature})

          MESSAGE(STATUS "DEAL_II_WITH_${_feature} has unmet external dependencies.")

          IF(FEATURE_${_feature}_HAVE_BUNDLED AND DEAL_II_ALLOW_BUNDLED)
            RUN_COMMAND("FEATURE_${_feature}_CONFIGURE_BUNDLED()")

            MESSAGE(STATUS "DEAL_II_WITH_${_feature} successfully set up with bundled packages.")
            SET(FEATURE_${_feature}_BUNDLED_CONFIGURED TRUE)
            SET_CACHED_OPTION(${_feature} ON)

          ELSE()
            IF(DEAL_II_WITH_${_feature})
              IF(COMMAND FEATURE_${_feature}_ERROR_MESSAGE)
                RUN_COMMAND("FEATURE_${_feature}_ERROR_MESSAGE()")
              ELSE()
                FEATURE_ERROR_MESSAGE(${_feature})
              ENDIF()
            ELSE()
              SET_CACHED_OPTION(${_feature} OFF)
            ENDIF()
          ENDIF()

        ENDIF(FEATURE_${_feature}_EXTERNAL_FOUND)

      ENDIF()
    ENDIF()
  ELSE()
    #
    # DEAL_II_WITH_${_feature} is defined and set to OFF, promote it to
    # cache nevertheless:
    #
    MESSAGE(STATUS "DEAL_II_WITH_${_feature} is set to off.")
    PURGE_FEATURE(${_feature})
    SET_CACHED_OPTION(${_feature} OFF)
  ENDIF()

  SET(FEATURE_${_feature}_PROCESSED TRUE)

ENDMACRO()
