IF(NOT DEAL_II_FORCE_CONTRIB_BOOST)

  IF(DEAL_II_ALLOW_CONTRIB)
    FIND_PACKAGE (Boost COMPONENTS serialization thread)
  ELSE()
    FIND_PACKAGE (Boost COMPONENTS serialization thread REQUIRED)
  ENDIF()

  # Get rid of this annoying unimportant variable:
  MARK_AS_ADVANCED(Boost_DIR)

  IF(Boost_THREAD_FOUND AND Boost_SERIALIZATION_FOUND)
    INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIR})

    LIST(APPEND deal_ii_external_libraries
      ${Boost_THREAD_LIBRARY} ${Boost_SERIALIZATION_LIBRARY}
      )
    LIST(APPEND deal_ii_external_debug_libraries
      ${Boost_THREAD_LIBRARY_DEBUG} ${Boost_SERIALIZATION_LIBRARY_DEBUG}
      )

    SET(DEAL_II_USE_EXTERNAL_BOOST TRUE)
  ENDIF()

ENDIF()

# ELSE() nothing to do. We compile the necessary boost source files
# directly if DEAL_II_USE_EXTERNAL_BOOST is not set.
