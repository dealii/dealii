ADD_EXECUTABLE(mesh_conversion
  Main.cc
  MeshConversion.cc
  )

INSTALL(TARGETS mesh_conversion
  RUNTIME
  DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
  COMPONENT contrib
  )

