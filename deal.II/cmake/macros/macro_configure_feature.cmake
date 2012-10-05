#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# This macro is used for the feature configuration in deal.II
#
# Usage:
#     CONFIGURE_FEATURE(feature)
#
#
# For a feature ${feature} (written in all caps) the following options,
# variables and macros have to be defined (except marked as optional):
#
# DEAL_II_WITH_${feature} (bool, mandatory)
#    If DEAL_II_FEATURE_AUTODETECTION is OFF, this boolean determines
#    whether the feature will be configured.
#    If DEAL_II_FEATURE_AUTODETECTION is ON, this boolean will
#    automatically be set if configuring the feature was successful.
#
# FEATURE_${feature}_DEPENDS (variable, optional)
#    a variable which contains an optional list of other features
#    this feature depends on (and which have to be enbled for this feature
#    to work.) The features must be given with the full option toggle:
#    DEAL_II_WITH_[...]
#
# FEATURE_${feature}_HAVE_BUNDLED  (variable, optional)
#    which should either be set to TRUE if all necessary libraries of the
#    features comes bundled with deal.II and hence can be supported
#    without external dependencies, or unset.
#
# FEATURE_${feature}_CONFIGURE_BUNDLED(var)  (macro, optional)
#    which should setup all necessary configuration for the feature with
#    bundled source dependencies. var set to TRUE indicates success,
#    otherwise this script gives an error.
#
# FEATURE_${feature}_FIND_EXTERNAL(var)  (macro, mandatory)
#    which should set var to TRUE if all dependencies for the feature are
#    fullfilled. In this case all necessary variables for
#    FEATURE_${feature}_CONFIGURE_EXTERNAL must be set. Otherwise
#    var should remain unset.
#    This macro should give an error (FATAL_ERROR or FATAL_ERROR).
#
# FEATURE_${feature}_CONFIGURE_EXTERNAL(var)  (macro, mandatory)
#    which should setup all necessary configuration for the feature with
#    external dependencies. var set to TRUE indicates success,
#    otherwise this script gives an error.
#
# FEATURE_${feature}_CUSTOM_ERROR_MESSAGE() (variable, optional)
#    which should either be set to TRUE if FEATURE_${feature}_ERROR_MESSAGE
#    is set up, or be undefined.
#
# FEATURE_${feature}_ERROR_MESSAGE()  (macro, optional)
#    which should print a meaningfull error message (with FATAL_ERROR) for
#    the case that no external library was found (and bundled is not
#    allowed to be used.) If not defined, a suitable default error message
#    will be printed.
#
#


#
# Some helper macros:
#

#
# Some black magic to have substitution in command names:
#
MACRO(RUN_COMMAND the_command)
  FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/macro_configure_feature.tmp"
    "${the_command}")
  INCLUDE("${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/macro_configure_feature.tmp")
ENDMACRO()


#
# A small macro to set the DEAL_II_WITH_${feature} variables:
#
MACRO(SET_CACHED_OPTION str value)
  STRING(TOLOWER "${str}" str_lower)
  SET(DEAL_II_WITH_${str}
    ${value}
    CACHE BOOL
    "Build deal.II with support for ${str_lower}."
    FORCE)
ENDMACRO()


#
# A small macro to post a default error message:
#
MACRO(FEATURE_ERROR_MESSAGE feature)
  STRING(TOLOWER ${feature} feature_lowercase)
  IF(FEATURE_${feature}_HAVE_BUNDLED)
    MESSAGE(FATAL_ERROR "\n"
      "Could not find the ${feature_lowercase} library!\n\n"
      "Please ensure that the ${feature_lowercase} library is installed on your computer.\n"
      "If the library is not at a default location, either provide some hints\n"
      "for the autodetection, or set the relevant variables by hand in ccmake.\n\n"
      "Alternatively you may choose to compile the bundled library of\n"
      "${feature_lowercase} by setting DEAL_II_ALLOW_BUNDLED=on or\n"
      "DEAL_II_FORCE_BUNDLED_${feature}=on.\n\n"
      )
 ELSE()
    MESSAGE(FATAL_ERROR "\n"
      "Could not find the ${feature_lowercase} library!\n\n"
      "Please ensure that the ${feature_lowercase} library is installed on your computer.\n"
      "If the library is not at a default location, either provide some hints\n"
      "for the autodetection, or set the relevant variables by hand in ccmake.\n\n"
      )
  ENDIF()
ENDMACRO()


MACRO(CONFIGURE_FEATURE feature)
  #
  # This script is arcane black magic. But at least for the better good: We
  # don't have to copy the configuration logic to every single
  # configure_<feature>.cmake script...
  #


  #
  # Check for correct include order of the configure_*.cmake files:
  # If feature B depends on feature A, configure_A.cmake has to be
  # included before configure_B.cmake:
  #
  FOREACH(macro_dependency ${FEATURE_${feature}_DEPENDS})
    STRING(REGEX REPLACE "^DEAL_II_WITH_" "" macro_dependency ${macro_dependency})
    IF(NOT FEATURE_${macro_dependency}_PROCESSED)
      MESSAGE(FATAL_ERROR "\n"
        "Internal build system error:\nDEAL_II_WITH_${feature} depends on "
        "DEAL_II_WITH_${macro_dependency},\nbut CONFIGURE_FEATURE(${feature}) "
        "was called before CONFIGURE_FEATURE(${macro_dependency}).\n\n"
        )
    ENDIF()
  ENDFOREACH()


  #
  # Obey the user overrides:
  #
  IF(FORCE_AUTODETECTION)
    UNSET(DEAL_II_WITH_${feature})
  ENDIF()
  IF(DISABLE_AUTODETECTION AND
    (NOT DEFINED DEAL_II_WITH_${feature}) )
    SET_CACHED_OPTION(${feature} OFF)
  ENDIF()


  #
  # Only try to configure ${feature} if we have to, i.e.
  # DEAL_II_WITH_${feature} is set to true or not set at all.
  #
  IF((NOT DEFINED DEAL_II_WITH_${feature}) OR
     DEAL_II_WITH_${feature})

    #
    # Are all dependencies fullfilled?
    #
    SET(macro_dependencies_ok TRUE)
    FOREACH(macro_dependency ${FEATURE_${feature}_DEPENDS})
      IF(NOT ${macro_dependency})
        IF(DEAL_II_WITH_${feature})
          MESSAGE(FATAL_ERROR "\n"
            "DEAL_II_WITH_${feature} has unmet configuration requirements: "
            "${macro_dependency} has to be set to \"ON\".\n\n"
            )
        ELSE()
          MESSAGE(STATUS
            "DEAL_II_WITH_${feature} has unmet configuration requirements: "
            "${macro_dependency} has to be set to \"ON\"."
            )
          SET_CACHED_OPTION(${feature} OFF)
        ENDIF()
        SET(macro_dependencies_ok FALSE)
      ENDIF()
    ENDFOREACH()

    IF(macro_dependencies_ok)
      IF(DEAL_II_FORCE_BUNDLED_${feature})
        #
        # First case: DEAL_II_FORCE_BUNDLED_${feature} is defined:
        #

        IF(FEATURE_${feature}_HAVE_BUNDLED)
          RUN_COMMAND(
            "FEATURE_${feature}_CONFIGURE_BUNDLED(FEATURE_${feature}_BUNDLED_CONFIGURED)"
            )
          IF(FEATURE_${feature}_BUNDLED_CONFIGURED)
            MESSAGE(STATUS
              "DEAL_II_WITH_${feature} successfully set up with bundled packages."
              )
            SET_CACHED_OPTION(${feature} ON)
          ELSE()
            # This should not happen. So give an error
            MESSAGE(FATAL_ERROR
              "\nInternal build system error: Failed to set up "
              "DEAL_II_WITH_${feature} with bundled packages.\n\n"
              )
          ENDIF()
        ELSE()
          MESSAGE(FATAL_ERROR "\n"
            "Internal build system error: DEAL_II_FORCE_BUNDLED_${feature} "
            "defined, but FEATURE_${feature}_HAVE_BUNDLED not present.\n"
            )
        ENDIF()

      ELSE(DEAL_II_FORCE_BUNDLED_${feature})
        #
        # Second case: We are allowed to search for an external library
        #

        RUN_COMMAND(
          "FEATURE_${feature}_FIND_EXTERNAL(FEATURE_${feature}_EXTERNAL_FOUND)"
          )

        IF(FEATURE_${feature}_EXTERNAL_FOUND)
          MESSAGE(STATUS
            "All external dependencies for DEAL_II_WITH_${feature} are fullfilled."
            )
          RUN_COMMAND(
            "FEATURE_${feature}_CONFIGURE_EXTERNAL(FEATURE_${feature}_EXTERNAL_CONFIGURED)"
            )

          IF(FEATURE_${feature}_EXTERNAL_CONFIGURED)
            MESSAGE(STATUS
              "DEAL_II_WITH_${feature} successfully set up with external dependencies."
              )
            SET_CACHED_OPTION(${feature} ON)
          ELSE()
            # This should not happen. So give an error
            MESSAGE(FATAL_ERROR
              "\nInternal build system error: Failed to set up "
              "DEAL_II_WITH_${feature} with external dependencies.\n\n"
              )
          ENDIF()

        ELSE(FEATURE_${feature}_EXTERNAL_FOUND)

          MESSAGE(STATUS
            "DEAL_II_WITH_${feature} has unmet external dependencies."
            )

          IF(FEATURE_${feature}_HAVE_BUNDLED AND DEAL_II_ALLOW_BUNDLED)
            RUN_COMMAND(
              "FEATURE_${feature}_CONFIGURE_BUNDLED(FEATURE_${feature}_BUNDLED_CONFIGURED)"
              )
            IF(FEATURE_${feature}_BUNDLED_CONFIGURED)
              MESSAGE(STATUS
                "DEAL_II_WITH_${feature} successfully set up with bundled packages."
                )
              SET_CACHED_OPTION(${feature} ON)
            ELSE()
              # This should not happen. So give an error
              MESSAGE(FATAL_ERROR
                "\nInternal build system error: Failed to set up "
                "DEAL_II_WITH_${feature} with bundled packages.\n\n"
                )
            ENDIF()
          ELSE()
            IF(DEAL_II_WITH_${feature})
              IF(FEATURE_${feature}_CUSTOM_ERROR_MESSAGE)
                RUN_COMMAND("FEATURE_${feature}_ERROR_MESSAGE()")
              ELSE()
                FEATURE_ERROR_MESSAGE(${feature})
              ENDIF()
            ELSE()
              SET_CACHED_OPTION(${feature} OFF)
            ENDIF()
          ENDIF()

        ENDIF(FEATURE_${feature}_EXTERNAL_FOUND)
      ENDIF()
    ENDIF()
  ELSE()
    #
    # DEAL_II_WITH_${feature} is defined and set to OFF, promote it to
    # cache nevertheless:
    #
    SET_CACHED_OPTION(${feature} OFF)
  ENDIF()

  SET(FEATURE_${feature}_PROCESSED TRUE)

ENDMACRO()
