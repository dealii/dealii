#
# TODO:
#
# A very primitive proof of concept on how to use doxygen, atm.
#
#
#  doxygen/tutorial/toc.html.in
#  doxygen/header.tex.in
#

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/options.dox.in
  ${CMAKE_CURRENT_BINARY_DIR}/options.dox
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/header.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/header.html
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/footer.html
  ${CMAKE_CURRENT_BINARY_DIR}/footer.html
  COPYONLY
  )

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/filter
  ${CMAKE_CURRENT_BINARY_DIR}/filter
  COPYONLY
  )



file(GLOB doxygen_input
  ${CMAKE_SOURCE_DIR}/include/deal.II/*
  ${CMAKE_BINARY_DIR}/include/deal.II/*
  ${CMAKE_SOURCE_DIR}/doc/news/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/*.h
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/headers
  )
TO_STRING(doxygen_input)


file(GLOB doxygen_image_path
  ${CMAKE_CURRENT_SOURCE_DIR}/images
  ${CMAKE_SOURCE_DIR}/examples/*/doc
  ${CMAKE_SOURCE_DIR}/contrib/parameter_gui/images
  )
TO_STRING(doxygen_image_path)


FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/options.dox"
  "
  INPUT=${doxygen_input}
  IMAGE_PATH=${doxygen_image_path}
  "
  )


ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/deal.II
  )


ADD_CUSTOM_TARGET(doxygen ALL # TODO: Really all?
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/options.dox
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/deal.II
    ${CMAKE_CURRENT_BINARY_DIR}/options.dox
    #already generated at configuration time:
    #${CMAKE_CURRENT_BINARY_DIR}/header.html
    #${CMAKE_CURRENT_BINARY_DIR}/footer.html
    #${CMAKE_CURRENT_BINARY_DIR}/filter
    #${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
    #${CMAKE_BINARY_DIR}/include/deal.II/lac/lapack_templates.h
  COMMENT "Generating documentation via doxygen"
  VERBATIM
  )

