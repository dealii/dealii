#
# This macro is used in the feature configuration mechanism of deal.II
#
# Usage:
#
# CONFIGURE_FEATURE(feature_name)
#
# This macro assumes the following macros and variables to be defined:
#
#
# HAVE_CONTRIB_FEATURE_${feature_name}  (variable)
#
#   - Which is either set to TRUE if all necessary libraries of the
#     features comes bundled with deal.II and hence can be supported
#     without external dependencies, or unset.
#
# CONFIGURE_FEATURE_${feature_name}_CONTRIB(var)  (macro)
#
#   - which shall setup all necessary configuration for the feature with
#     contrib source dependencies. var set to TRUE shall indicate success.
#
#
# FIND_FEATURE_${feature_name}_EXTERNAL(var)  (macro)
#
#   - which shall set var to TRUE if all dependencies for the feature are
#     fullfilled. In this case all necessary variables for
#     CONFIGURE_FEATURE_${feature_name}_EXTERNAL shall be set. Otherwise
#     var should remain unset.
#
# CONFIGURE_FEATURE_${feature_name}_EXTERNAL(var)  (macro)
#
#   - which shall setup all necessary configuration for the feature with
#     external dependencies. var set to TRUE shall indicate success.
#
# CONFIGURE_FEATURE_${feature_name}_ERROR_MESSAGE()  (macro)
#
#   - which shall print a meaningfull error message if case no external
#     library is found (and contrib is not allowed to be used.)
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


MACRO(CONFIGURE_FEATURE feature)

  # Only (try) to configure ${feature} if we have to:
  IF(DEAL_II_FEATURE_AUTODETECTION OR DEAL_II_WITH_${feature})

    #
    # First case: DEAL_II_FORCE_CONTRIB_${feature} is defined:
    #
    IF(DEAL_II_FORCE_CONTRIB_${feature})

      IF(HAVE_CONTRIB_FEATURE_${feature})
        RUN_COMMAND(
          "
          CONFIGURE_FEATURE_${feature}_CONTRIB(
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
          "Internal build system error: DEAL_II_FORCE_CONTRIB_${feature} defined, but HAVE_CONTRIB_FEATURE_${feature} not present."
          )
      ENDIF()

    ELSE()

      #
      # Second case: We are allowed to search for an external library:
      #
      RUN_COMMAND(
        "FIND_FEATURE_${feature}_EXTERNAL(FEATURE_${feature}_EXTERNAL_FOUND)"
        )

      IF(FEATURE_${feature}_EXTERNAL_FOUND)

        MESSAGE(STATUS
          "All external dependencies for DEAL_II_WITH_${feature} are fullfilled."
          )

        RUN_COMMAND(
          "
          CONFIGURE_FEATURE_${feature}_EXTERNAL(
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

        IF(HAVE_CONTRIB_FEATURE_${feature} AND DEAL_II_ALLOW_CONTRIB)
          RUN_COMMAND(
            "
            CONFIGURE_FEATURE_${feature}_CONTRIB(
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
              SET_CACHED_OPTION(DEAL_II_WITH_${feature} ON)
            ELSE()
              RUN_COMMAND(
                "CONFIGURE_FEATURE_${feature}_ERROR_MESSAGE()"
                )
            ENDIF()
        ENDIF()

      ENDIF()
    ENDIF()

  ENDIF()
ENDMACRO()
