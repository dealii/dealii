#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
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
# Configuration for the boost library:
#

#
# This configure script has to be included after configure_threads.
# We need some of the variables defined in SETUP_THREADING for
# the setup of the bundled boost library (if used)
#
IF(NOT FEATURE_THREADS_PROCESSED)
  MESSAGE(FATAL_ERROR "\n"
    "Internal build system error:\n"
    "configure_boost.cmake included before configure_tbb.cmake\n\n"
    )
ENDIF()


SET(DEAL_II_WITH_BOOST ON # Always true. We need it :-]
  CACHE BOOL "Build deal.II with support for boost." FORCE
  )


MACRO(FEATURE_BOOST_FIND_EXTERNAL var)
  #
  # We require at least version 1.44
  #
  IF(DEAL_II_WITH_THREADS)
    FIND_PACKAGE(Boost 1.44 COMPONENTS serialization system thread)
  ELSE()
    FIND_PACKAGE(Boost 1.44 COMPONENTS serialization system)
  ENDIF()

  IF( Boost_SERIALIZATION_FOUND AND
      Boost_SYSTEM_FOUND AND
      (NOT DEAL_II_WITH_THREADS OR Boost_THREAD_FOUND)
    )
    SET(${var} TRUE)
    MARK_AS_ADVANCED(Boost_DIR)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_BOOST_CONFIGURE_EXTERNAL)
  INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})

  # The user has to know the location of the boost headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})

  #
  # Remove "pthread" from Boost_LIBRARIES. Threading, if necessary, is
  # already set up via configure_1_threads.cmake.
  #
  LIST(REMOVE_ITEM Boost_LIBRARIES "pthread")

  #
  # Transform  Boost_LIBRARIES into a list of debug and release libraries
  # without keywords:
  #
  SPLIT_DEBUG_RELEASE(_boost_debug _boost_release ${Boost_LIBRARIES})

  IF (CMAKE_BUILD_TYPE MATCHES "Debug")
    LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_DEBUG ${_boost_debug})
  ENDIF()

  IF (CMAKE_BUILD_TYPE MATCHES "Release")
    LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_RELEASE ${_boost_release})
  ENDIF()
ENDMACRO()


MACRO(FEATURE_BOOST_CONFIGURE_BUNDLED)
  #
  # We need to set some definitions to use the headers of the bundled boost
  # library:
  #
  LIST(APPEND DEAL_II_DEFINITIONS
    "BOOST_NO_HASH" "BOOST_NO_SLIST"
    )
  LIST(APPEND DEAL_II_USER_DEFINITIONS
    "BOOST_NO_HASH" "BOOST_NO_SLIST"
    )

  INCLUDE_DIRECTORIES(${BOOST_FOLDER}/include)
ENDMACRO()


CONFIGURE_FEATURE(BOOST)

#
# DEAL_II_WITH_BOOST is always required.
#
IF(NOT DEAL_II_WITH_BOOST)
  IF(DEAL_II_FEATURE_AUTODETECTION)
    FEATURE_ERROR_MESSAGE("BOOST")
  ELSE()
    MESSAGE(FATAL_ERROR "\n"
      "Unmet configuration requirements: "
      "DEAL_II_WITH_BOOST required, but set to OFF!.\n\n"
      )
  ENDIF()
ENDIF()
