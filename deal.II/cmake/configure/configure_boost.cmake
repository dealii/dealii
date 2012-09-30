#
# Configuration for the boost library:
#


#
# This configure script has to be included after configure_tbb.
# We need some of the variables defined in SETUP_THREADING for
# the setup of the bundled boost library (if used)
#
IF(NOT FEATURE_TBB_PROCESSED)
  MESSAGE(FATAL_ERROR "\n"
    "Internal build system error:\n"
    "configure_boost.cmake included before configure_tbb.cmake\n\n"
    )
ENDIF()


MACRO(FEATURE_BOOST_FIND_EXTERNAL var)
  FIND_PACKAGE (Boost COMPONENTS serialization thread)

  IF(Boost_THREAD_FOUND AND Boost_SERIALIZATION_FOUND)
    SET(${var} TRUE)

    # Get rid of this annoying unimportant variable:
    MARK_AS_ADVANCED(Boost_DIR)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_BOOST_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIR})

  IF (CMAKE_BUILD_TYPE MATCHES "Debug")
    LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_DEBUG
      ${Boost_THREAD_LIBRARY_DEBUG} ${Boost_SERIALIZATION_LIBRARY_DEBUG}
      )
  ENDIF()

  IF (CMAKE_BUILD_TYPE MATCHES "Release")
    LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_RELEASE
      ${Boost_THREAD_LIBRARY_RELEASE} ${Boost_SERIALIZATION_LIBRARY_RELEASE}
      )
  ENDIF()

  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_BOOST_HAVE_BUNDLED TRUE)


MACRO(FEATURE_BOOST_CONFIGURE_BUNDLED var)

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

  SET(boost_folder "${CMAKE_SOURCE_DIR}/bundled/boost-1.49.0")
  INCLUDE_DIRECTORIES(${boost_folder}/include)

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(BOOST)


#
# DEAL_II_WITH_BOOST is always required.
#
IF(NOT DEAL_II_WITH_BOOST)
  IF(DEAL_II_FEATURE_AUTODETECTION)
    FEATURE_ERROR_MESSAGE("BOOST")
  ELSE()
    MESSAGE(SEND_ERROR "\n"
      "Unmet configuration requirements: "
      "DEAL_II_WITH_BOOST required, but set to OFF!.\n\n"
      )
  ENDIF()
ENDIF()

