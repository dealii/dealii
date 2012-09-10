IF(DEAL_II_ALLOW_CONTRIB)
  FIND_PACKAGE (Boost COMPONENTS serialization thread)
ELSE()
  FIND_PACKAGE (Boost COMPONENTS serialization thread REQUIRED)
ENDIF()

#
# Get rid of this annoying unimportant variable:
#
MARK_AS_ADVANCED(Boost_DIR)

IF(NOT DEAL_II_FORCE_CONTRIB_BOOST)
  IF(Boost_THREAD_FOUND AND Boost_SERIALIZATION_FOUND)
    INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIR})

    LIST(APPEND deal_ii_external_libraries
      ${Boost_THREAD_LIBRARY} ${Boost_SERIALIZATION_LIBRARY}
      )
    LIST(APPEND deal_ii_external_debug_libraries
      ${Boost_THREAD_LIBRARY_DEBUG} ${Boost_SERIALIZATION_LIBRARY_DEBUG}
      )

    # TODO: Renaming!
    SET(DEAL_II_USE_EXTERNAL_BOOST TRUE)
  ENDIF()
ENDIF()
