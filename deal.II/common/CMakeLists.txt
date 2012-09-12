ADD_EXECUTABLE(expand_instantiations expand_instantiations/expand_instantiations.cc)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
  ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
  )
