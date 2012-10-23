#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Configuration for the MUMPS library:
#


SET(FEATURE_MUMPS_DEPENDS DEAL_II_WITH_MPI)


MACRO(FEATURE_MUMPS_FIND_EXTERNAL var)
  FIND_PACKAGE(MUMPS)

  IF(MUMPS_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_MUMPS_CONFIGURE_EXTERNAL var)
  INCLUDE_DIRECTORIES(${MUMPS_INCLUDE_DIRS})
  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES
    ${MUMPS_LIBRARIES}
    ${MPI_CXX_LIBRARIES} # for good measure
    )
  ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${MUMPS_LINKER_FLAGS}")

  SET(DEAL_II_USE_MUMPS TRUE)

  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_MUMPS_CUSTOM_ERROR_MESSAGE TRUE)


MACRO(FEATURE_MUMPS_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find the mumps library!\n"
    "Please ensure that the library is installed on your computer.\n"
    "If the libraries is not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ MUMPS_DIR=\"...\" cmake <...>\n"
    "    $ cmake -DMUMPS_DIR=\"...\" <...>\n"
    "or set the relevant variables by hand in ccmake.\n"
    "Relevant hints for MUMPS are MUMPS_DIR, SCALAPACK_DIR (and BLACS_DIR).\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(MUMPS)
