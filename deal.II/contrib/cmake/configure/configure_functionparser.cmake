MACRO(FIND_FEATURE_FUNCTIONPARSER_EXTERNAL var)
ENDMACRO()


MACRO(CONFIGURE_FEATURE_FUNCTIONPARSER_EXTERNAL var)
ENDMACRO()


SET(HAVE_CONTRIB_FEATURE_FUNCTIONPARSER TRUE)


MACRO(CONFIGURE_FEATURE_FUNCTIONPARSER_CONTRIB var)

  # compile the necessary parts of functionparser out of ./contrib

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

  SET(HAVE_FUNCTIONPARSER TRUE)

  SET(${var} TRUE)
ENDMACRO()


MACRO(CONFIGURE_FEATURE_FUNCTIONPARSER_ERROR_MESSAGE)
    MESSAGE(SEND_ERROR "
Could not find the functionparser library!

Module not written, yet...

You may want to choose to compile the bundled contrib library of
functionparser by setting DEAL_II_ALLOW_CONTRIB=on or
DEAL_II_FORCE_CONTRIB_FUNCTIONPARSER=on.

")
ENDMACRO()


CONFIGURE_FEATURE(FUNCTIONPARSER)
