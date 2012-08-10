IF(NOT DEAL_II_USE_CONTRIB)
  FIND_PACKAGE (Boost COMPONENTS serialization thread REQUIRED)

  INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIR})

  SET(deal_ii_external_libraries
    ${deal_ii_external_libraries}
    ${Boost_THREAD_LIBRARY} ${Boost_SERIALIZATION_LIBRARY}
    )
  SET(deal_ii_external_debug_libraries
    ${deal_ii_external_debug_libraries}
    ${Boost_THREAD_LIBRARY_DEBUG} ${Boost_SERIALIZATION_LIBRARY_DEBUG}
    ) #TODO

  SET(DEAL_II_USE_EXTERNAL_BOOST TRUE)
ENDIF()
