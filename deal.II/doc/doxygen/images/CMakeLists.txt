ADD_CUSTOM_COMMAND(
  OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/structure.png
    ${CMAKE_CURRENT_BINARY_DIR}/structure.cmapx
  COMMAND ${DOXYGEN_DOT_EXECUTABLE}
    -Tpng -o ${CMAKE_CURRENT_BINARY_DIR}/structure.png
    -Tcmapx -o ${CMAKE_CURRENT_BINARY_DIR}/structure.cmapx
    ${CMAKE_CURRENT_SOURCE_DIR}/structure.dot
  DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/structure.dot
  )

ADD_CUSTOM_TARGET(images
  DEPENDS
    ${CMAKE_CURRENT_BINARY_DIR}/structure.png
    ${CMAKE_CURRENT_BINARY_DIR}/structure.cmapx
  )
