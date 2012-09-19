#
# This macro is used in the feature configuration mechanism of deal.II
#
# Usage:
#
# CONFIGURE_FEATURE(feature)
#
#
# This macro assumes the following macros and variables to be defined:
#
# FEATURE_${feature}_DEPENDS (variable, optional)
#
#    a variable which contains an optional list of features this feature
#    depends on (and which have to be enbled for this feature to work.)
#
#
# FEATURE_${feature}_HAVE_CONTRIB  (variable, optional)
#
#    which shall either be set to TRUE if all necessary libraries of the
#    features comes bundled with deal.II and hence can be supported
#    without external dependencies, or unset.
#
# FEATURE_${feature}_CONFIGURE_CONTRIB(var)  (macro, optional)
#
#    which shall setup all necessary configuration for the feature with
#    contrib source dependencies. var set to TRUE shall indicate success.
#
#
# FEATURE_${feature}_FIND_EXTERNAL(var)  (macro)
#
#    which shall set var to TRUE if all dependencies for the feature are
#    fullfilled. In this case all necessary variables for
#    FEATURE_${feature}_CONFIGURE_EXTERNAL must be set. Otherwise
#    var should remain unset.
#
# FEATURE_${feature}_CONFIGURE_EXTERNAL(var)  (macro)
#
#    which shall setup all necessary configuration for the feature with
#    external dependencies. var set to TRUE shall indicate success.
#
#
# FEATURE_${feature}_CUSTOM_ERROR_MESSAGE() (variable, optional)
#
#    which shall either be set to TRUE if FEATURE_${feature}_ERROR_MESSAGE
#    is set up, or be undefined.
#
# FEATURE_${feature}_ERROR_MESSAGE()  (macro, optional)
#
#    which shall print a meaningfull error message if case no external
#    library is found (and contrib is not allowed to be used.)
#    If not defined, a suitable default error message will be printed.
#
#

#
# FAT NOTE:
# This script is arcane black magic. But at least for the better good: We
# don't have to copy the configuration logic to every single
# configure_<feature>.cmake script.
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
  SET(${str}
    ${value}
    CACHE STRING
    "Automatically set due to DEAL_II_FEATURE_AUTODETECTION"
    FORCE)
ENDMACRO()


#
# A small macro to post an error message:
#
MACRO(FEATURE_ERROR_MESSAGE feature)
  STRING(TOLOWER ${feature} feature_lowercase)
  IF(FEATURE_${feature}_HAVE_CONTRIB)
    MESSAGE(SEND_ERROR "\n"
      "Could not find the ${feature_lowercase} library!\n\n"
      "Please ensure that the ${feature_lowercase} library is installed on your computer.\n"
      "If the library is not at a default location, either provide some hints\n"
      "for the autodetection, or set the relevant variables by hand in ccmake.\n\n"
      "Alternatively you may choose to compile the bundled contrib library of\n"
      "${feature_lowercase} by setting DEAL_II_ALLOW_CONTRIB=on or\n"
      "DEAL_II_FORCE_CONTRIB_${feature}=on.\n\n"
      )
 ELSE()
    MESSAGE(SEND_ERROR "\n"
      "Could not find the ${feature_lowercase} library!\n\n"
      "Please ensure that the ${feature_lowercase} library is installed on your computer.\n"
      "If the library is not at a default location, either provide some hints\n"
      "for the autodetection, or set the relevant variables by hand in ccmake.\n\n"
      )
  ENDIF()
ENDMACRO()


MACRO(CONFIGURE_FEATURE feature)

  #
  # Check for correct include order of the configure_*.cmake files:
  # If feature B depends on feature A, configure_A.cmake has to be
  # included before configure_B.cmake:
  #
  FOREACH(macro_dependency ${FEATURE_${feature}_DEPENDS})
    STRING(REGEX REPLACE "^DEAL_II_WITH_" "" macro_dependency ${macro_dependency})
    MESSAGE("${macro_dependency}")
    IF(NOT FEATURE_${macro_dependency}_PROCESSED)
      MESSAGE(FATAL_ERROR "\n"
        "Internal build system error:\nDEAL_II_WITH_${feature} depends on "
        "DEAL_II_WITH_${macro_dependency},\nbut CONFIGURE_FEATURE(${feature}) "
        "was called before CONFIGURE_FEATURE(${macro_dependency}).\n\n"
        )
    ENDIF()
  ENDFOREACH()


  # Only try to configure ${feature} if we have to:
  IF(DEAL_II_FEATURE_AUTODETECTION OR DEAL_II_WITH_${feature})

    #
    # Are all dependencies fullfilled?
    #
    SET(macro_dependencies_ok TRUE)

    FOREACH(macro_dependency ${FEATURE_${feature}_DEPENDS})
      IF(NOT ${macro_dependency})

        IF(DEAL_II_FEATURE_AUTODETECTION)
          MESSAGE(STATUS
            "DEAL_II_WITH_${feature} has unmet configuration requirements: "
            "${macro_dependency} has to be set to \"ON\"."
            )
          SET_CACHED_OPTION(DEAL_II_WITH_${feature} OFF)
        ELSE()
          MESSAGE(SEND_ERROR
            "DEAL_II_WITH_${feature} has unmet configuration requirements: "
            "${macro_dependency} has to be set to \"ON\"."
            )
        ENDIF()

        SET(macro_dependencies_ok FALSE)
      ENDIF()
    ENDFOREACH()


    IF(macro_dependencies_ok)
      #
      # First case: DEAL_II_FORCE_CONTRIB_${feature} is defined:
      #
      IF(DEAL_II_FORCE_CONTRIB_${feature})

        IF(FEATURE_${feature}_HAVE_CONTRIB)
          RUN_COMMAND(
            "
            FEATURE_${feature}_CONFIGURE_CONTRIB(
              FEATURE_${feature}_CONTRIB_CONFIGURED
              )
            "
            )
          IF(FEATURE_${feature}_CONTRIB_CONFIGURED)
            MESSAGE(STATUS
              "DEAL_II_WITH_${feature} successfully set up with contrib packages."
              )
            IF(DEAL_II_FEATURE_AUTODETECTION)
              SET_CACHED_OPTION(DEAL_II_WITH_${feature} ON)
            ENDIF()
          ELSE()
            # This should not happen. So give an error
            MESSAGE(SEND_ERROR
              "Failed to set up DEAL_II_WITH_${feature} with contrib packages."
              )
          ENDIF()
        ELSE()
          MESSAGE(FATAL_ERROR
            "Internal build system error: DEAL_II_FORCE_CONTRIB_${feature} "
            "defined, but FEATURE_${feature}_HAVE_CONTRIB not present."
            )
        ENDIF()

      ELSE(DEAL_II_FORCE_CONTRIB_${feature})

        #
        # Second case: We are allowed to search for an external library:
        #
        RUN_COMMAND(
          "FEATURE_${feature}_FIND_EXTERNAL(FEATURE_${feature}_EXTERNAL_FOUND)"
          )

        IF(FEATURE_${feature}_EXTERNAL_FOUND)

          MESSAGE(STATUS
            "All external dependencies for DEAL_II_WITH_${feature} are fullfilled."
            )

          RUN_COMMAND(
            "
            FEATURE_${feature}_CONFIGURE_EXTERNAL(
              FEATURE_${feature}_EXTERNAL_CONFIGURED
              )
            "
            )

          IF(FEATURE_${feature}_EXTERNAL_CONFIGURED)
            MESSAGE(STATUS
              "DEAL_II_WITH_${feature} successfully set up with external dependencies."
              )
            IF(DEAL_II_FEATURE_AUTODETECTION)
              SET_CACHED_OPTION(DEAL_II_WITH_${feature} ON)
            ENDIF()
          ELSE()
            # This should not happen. So give an error
            MESSAGE(SEND_ERROR
              "Failed to set up DEAL_II_WITH_${feature} with external dependencies."
              )
          ENDIF()

        ELSE()

          MESSAGE(STATUS
            "DEAL_II_WITH_${feature} has unmet external dependencies."
            )

          IF(FEATURE_${feature}_HAVE_CONTRIB AND DEAL_II_ALLOW_CONTRIB)
            RUN_COMMAND(
              "
              FEATURE_${feature}_CONFIGURE_CONTRIB(
                FEATURE_${feature}_CONTRIB_CONFIGURED
                )
              "
              )
            IF(FEATURE_${feature}_CONTRIB_CONFIGURED)
              MESSAGE(STATUS
                "DEAL_II_WITH_${feature} successfully set up with contrib packages."
                )
              IF(DEAL_II_FEATURE_AUTODETECTION)
                SET_CACHED_OPTION(DEAL_II_WITH_${feature} ON)
              ENDIF()
            ELSE()
              # This should not happen. So give an error
              MESSAGE(SEND_ERROR
                "Failed to set up DEAL_II_WITH_${feature} with contrib packages."
                )
            ENDIF()
          ELSE()
              IF(DEAL_II_FEATURE_AUTODETECTION)
                SET_CACHED_OPTION(DEAL_II_WITH_${feature} OFF)
              ELSE()
                IF(FEATURE_${feature}_CUSTOM_ERROR_MESSAGE)
                  RUN_COMMAND(
                    "
                    FEATURE_${feature}_ERROR_MESSAGE()
                    "
                    )
                ELSE()
                  FEATURE_ERROR_MESSAGE(${feature})
                ENDIF()
              ENDIF()
          ENDIF()
        ENDIF()
      ENDIF(DEAL_II_FORCE_CONTRIB_${feature})

    ENDIF()

  ENDIF()

  SET(FEATURE_${feature}_PROCESSED TRUE)

ENDMACRO()
