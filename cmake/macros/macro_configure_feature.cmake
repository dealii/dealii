## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
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
# This macro is used for the feature configuration in deal.II
#
# Usage:
#     configure_feature(feature)
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
#    If not defined, find_package(${feature}) is called.
#
# FEATURE_${feature}_CONFIGURE_EXTERNAL()   (macro)
#    which should setup all necessary configuration for the feature with
#    external dependencies.
#
# FEATURE_${feature}_ERROR_message()  (macro)
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
# A small macro to set the DEAL_II_WITH_${_feature} variables:
#
macro(set_cached_option _str _value)
  string(TOLOWER "${_str}" _str_lower)
  set(DEAL_II_WITH_${_str}
    ${_value}
    CACHE BOOL
    "Build deal.II with support for ${_str_lower}."
    FORCE)
endmacro()


#
# A small macro to post a default error message:
#
macro(feature_error_message _feature)
  string(TOLOWER ${_feature} _feature_lowercase)

  if(DEFINED ${_feature}_DIR)
    set(_hint_snippet "
    $ ${_feature}_DIR=\"...\" cmake <...>
    $ cmake -D${_feature}_DIR=\"...\" <...>
or set the relevant variables by hand in ccmake."
      )
  else()
    set(_hint_snippet
      " or set the relevant variables by hand in ccmake."
      )
  endif()

  if(FEATURE_${_feature}_HAVE_BUNDLED)
    set(_bundled_snippet
      "\nAlternatively you may choose to compile the bundled library of "
      "${_feature_lowercase} by setting DEAL_II_ALLOW_BUNDLED=on or "
      "DEAL_II_FORCE_BUNDLED_${_feature}=on.\n"
      )
  else()
    set(_bundled_snippet "\n")
  endif()

  message(FATAL_ERROR "\n"
    "Could not find the ${_feature_lowercase} library!\n"
    ${${_feature}_ADDITIONAL_ERROR_STRING}
    "Please ensure that a suitable ${_feature_lowercase} library is installed on your computer.\n"
    "If the library is not at a default location, either provide some hints "
    "for autodetection,${_hint_snippet}"
    ${_bundled_snippet}
    )
endmacro()


#
# Default macro for finding an external library:
#
macro(feature_find_external _feature _var)
  find_package(DEAL_II_${_feature})
  if(${_feature}_FOUND)
    set(${_var} TRUE)
  endif()
endmacro()


########################################################################
#                                                                      #
#                          CONFIGURE_FEATURE:                          #
#                                                                      #
########################################################################

macro(configure_feature _feature)

  #
  # Register the feature in the DEAL_II_FEATURES list
  #
  list(APPEND DEAL_II_FEATURES ${_feature})

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
  foreach(_dependency
      ${FEATURE_${_feature}_AFTER}
      ${FEATURE_${_feature}_DEPENDS}
      )
    if(NOT FEATURE_${_dependency}_PROCESSED)
      message(FATAL_ERROR "\n"
        "Internal build system error: The configuration of "
        "DEAL_II_WITH_${_feature} depends on "
        "DEAL_II_WITH_${_dependency}, but configure_feature(${_feature}) "
        "was called before configure_feature(${_dependency}).\n\n"
        )
    endif()
  endforeach()

  #
  # Obey the user overrides:
  #
  if( (NOT DEAL_II_ALLOW_AUTODETECTION) AND
      (NOT DEFINED DEAL_II_WITH_${_feature}) )
    purge_feature(${_feature})
    set_cached_option(${_feature} OFF)
  endif()

  #
  # Only try to configure ${_feature} if we have to, i.e.
  # DEAL_II_WITH_${_feature} is set to true or not set at all.
  #
  if((NOT DEFINED DEAL_II_WITH_${_feature}) OR
     DEAL_II_WITH_${_feature})

    #
    # Are all dependencies fulfilled?
    #
    set(_dependencies_ok TRUE)
    foreach(_dependency ${FEATURE_${_feature}_DEPENDS})
      if(NOT DEAL_II_WITH_${_dependency})
        if(DEAL_II_WITH_${_feature})
          message(FATAL_ERROR "\n"
            "DEAL_II_WITH_${_feature} has unmet configuration requirements: "
            "DEAL_II_WITH_${_dependency} has to be set to \"ON\".\n\n"
            )
        else()
          message(STATUS
            "DEAL_II_WITH_${_feature} has unmet configuration requirements: "
            "DEAL_II_WITH_${_dependency} has to be set to \"ON\"."
            )
          purge_feature(${_feature})
          set_cached_option(${_feature} OFF)
        endif()
        set(_dependencies_ok FALSE)
      endif()
    endforeach()

    if(_dependencies_ok)
      if(DEAL_II_FORCE_BUNDLED_${_feature})
        #
        # First case: DEAL_II_FORCE_BUNDLED_${_feature} is defined:
        #

        purge_feature(${_feature})

        if(FEATURE_${_feature}_HAVE_BUNDLED)
          evaluate_expression("feature_${_feature}_configure_bundled()")
          message(STATUS "")
          message(STATUS "DEAL_II_WITH_${_feature} successfully set up with bundled packages.")
          set(DEAL_II_FEATURE_${_feature}_BUNDLED_CONFIGURED TRUE)
          set_cached_option(${_feature} ON)
        else()
          message(FATAL_ERROR "\n"
            "Internal build system error: DEAL_II_FORCE_BUNDLED_${_feature} "
            "defined, but FEATURE_${_feature}_HAVE_BUNDLED not present.\n"
            )
        endif()

      else()
        #
        # Second case: We are allowed to search for an external library
        #
        if(NOT FEATURE_${_feature}_EXTERNAL_FOUND AND NOT DEAL_II_${_feature}_FOUND)
          if(COMMAND FEATURE_${_feature}_FIND_EXTERNAL)
            evaluate_expression("feature_${_feature}_find_external(FEATURE_${_feature}_EXTERNAL_FOUND)")
          else()
            feature_find_external(${_feature} FEATURE_${_feature}_EXTERNAL_FOUND)
          endif()
        else()
          set(FEATURE_${_feature}_EXTERNAL_FOUND TRUE)
        endif()

        if(FEATURE_${_feature}_EXTERNAL_FOUND)
          if(COMMAND feature_${_feature}_configure_external)
            evaluate_expression("feature_${_feature}_configure_external()")
          endif()
          define_interface_target(${_feature})

          message(STATUS "")
          message(STATUS "DEAL_II_WITH_${_feature} successfully set up with external dependencies.")
          set(FEATURE_${_feature}_EXTERNAL_CONFIGURED TRUE)
          set_cached_option(${_feature} ON)

        else()

          purge_feature(${_feature})

          message(STATUS "DEAL_II_WITH_${_feature} has unmet external dependencies.")

          if(FEATURE_${_feature}_HAVE_BUNDLED AND DEAL_II_ALLOW_BUNDLED)
            evaluate_expression("feature_${_feature}_configure_bundled()")

            message(STATUS "")
            message(STATUS "DEAL_II_WITH_${_feature} successfully set up with bundled packages.")
            set(DEAL_II_FEATURE_${_feature}_BUNDLED_CONFIGURED TRUE)
            set_cached_option(${_feature} ON)

          else()
            if(DEAL_II_WITH_${_feature})
              if(COMMAND feature_${_feature}_error_message)
                 evaluate_expression("feature_${_feature}_error_message()")
              else()
                feature_error_message(${_feature})
              endif()
            else()
              set_cached_option(${_feature} OFF)
            endif()
          endif()

        endif()

      endif()
    endif()
  else()
    #
    # DEAL_II_WITH_${_feature} is defined and set to OFF, promote it to
    # cache nevertheless:
    #
    message(STATUS "DEAL_II_WITH_${_feature} is set to off.")
    purge_feature(${_feature})
    set_cached_option(${_feature} OFF)
  endif()

  set(FEATURE_${_feature}_PROCESSED TRUE)

endmacro()
