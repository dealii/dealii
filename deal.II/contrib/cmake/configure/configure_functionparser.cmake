IF(NOT DEAL_II_FORCE_CONTRIB_FUNCTIONPARSER)
  IF(DEAL_II_ALLOW_CONTRIB)
    # TODO: Write a module to search for functionparser
  ELSE()
    MESSAGE(FATAL_ERROR "FindFunctionparser.cmake not written, yet. :-[")
  ENDIF()
ENDIF()

IF(DEAL_II_FORCE_CONTRIB_FUNCTIONPARSER OR NOT Functionparser_FOUND)

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/contrib/functionparser/
    )

  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/contrib/functionparser)

  #
  # Add functionparser directly to the object files of deal.II
  #
  LIST(APPEND deal_ii_additional_object_files
    $<TARGET_OBJECTS:obj_functionparser>
    )

ELSE()

  INCLUDE_DIRECTORIES(${Functionparser_INCLUDE_DIR})

  LIST(APPEND deal_ii_external_libraries ${Functionparser_LIBRARY})
  LIST(APPEND deal_ii_external_debug_libraries ${Functionparser_LIBRARY})

ENDIF()

SET(HAVE_FUNCTIONPARSER TRUE)
