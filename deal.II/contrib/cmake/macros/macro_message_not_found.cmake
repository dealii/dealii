FUNCTION(macro_message_not_found library library_uppercase)

  # We have contrib source and our own Find<...> module:
  SET(contrib_and_find amd tbb umfpack)

  # We have contrib source and use an external Find<...> module:
  SET(contrib boost)

  # We use our own Find<...> module:
  SET(find_module netcdf)


  IF(library_uppercase MATCHES AMD)
    SET(library_contrib UMFPACK)
  ELSE()
    SET(library_contrib ${library_uppercase})
  ENDIF()


  LIST_CONTAINS(result library ${contrib_and_find})
  IF(result)
    MESSAGE(SEND_ERROR "
Could not find the ${library} library!

Please ensure that the ${library} library is installed on your computer.
If the library is not at a default location, either provide some hints
via environment variables:
${library_uppercase}_LIBRARY_DIR ${library_uppercase}_INCLUDE_DIR
Or set the relevant variables by hand in ccmake.

Alternatively you may choose to compile the bundled contrib library of
${library} by setting DEAL_II_ALLOW_CONTRIB=on or
DEAL_II_FORCE_CONTRIB_${library_contrib}=on.

")
    RETURN()
  ENDIF()


  LIST_CONTAINS(result library ${contrib})
  IF(result)
    MESSAGE(SEND_ERROR "
Could not find the ${library} library!

Please ensure that the ${library} library is installed on your computer.
If the library is not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

Alternatively you may choose to compile the bundled contrib library of
${library} by setting DEAL_II_ALLOW_CONTRIB=on or
DEAL_II_FORCE_CONTRIB_${library_contrib}=on.

")
    RETURN()
  ENDIF()


  LIST_CONTAINS(result library ${find_module})
  IF(result)
    MESSAGE(SEND_ERROR "
Could not find the ${library} library!

Please ensure that the ${library} library is installed on your computer.
If the library is not at a default location, either provide some hints
via environment variables:
${library_uppercase}_LIBRARY_DIR ${library_uppercase}_INCLUDE_DIR
Or set the relevant variables by hand in ccmake.

")
    RETURN()
  ENDIF()


    MESSAGE(SEND_ERROR "
Could not find the ${library} library!

Please ensure that the ${library} library is installed on your computer.
If the library is not at a default location, either provide some hints
for the autodetection, or set the relevant variables by hand in ccmake.

")

ENDFUNCTION()
