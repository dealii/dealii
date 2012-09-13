MACRO(macro_message_not_found library)
  # TODO: Check for each Find<...> module how to hint for a library
  # location...

  STRING(TOUPPER "${library}" library_uppercase)

  MESSAGE(SEND_ERROR "
Could not find the ${library} library!

Please ensure that the ${library} library is installed on your computer.
If the library is not at a default location, either provide some hints
via environment variables:
${library_uppercase}_LIBRARY_DIR ${library_uppercase}_INCLUDE_DIR
Or set the relevant variables by hand in ccmake.

")

ENDMACRO()
