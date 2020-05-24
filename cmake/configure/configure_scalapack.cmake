## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2019 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Configuration for the SCALAPACK library:
#

SET(FEATURE_SCALAPACK_DEPENDS MPI LAPACK)


MACRO(FEATURE_SCALAPACK_FIND_EXTERNAL var)
  FIND_PACKAGE(SCALAPACK)

  IF(SCALAPACK_FOUND)
    SET(${var} TRUE)
    CHECK_MPI_INTERFACE(SCALAPACK ${var})

    IF (${var})
      CLEAR_CMAKE_REQUIRED()
      SET(CMAKE_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})
      CHECK_C_SOURCE_COMPILES("
        void pdsyevr_();
        void pssyevr_();
        int main(){
          pdsyevr_();
          pssyevr_();
          return 0;
        }"
        DEAL_II_SCALAPACK_HAS_PDSYEVR_PSSYEVR)
        RESET_CMAKE_REQUIRED()

      IF(NOT DEAL_II_SCALAPACK_HAS_PDSYEVR_PSSYEVR)
        MESSAGE(STATUS "Could not find a sufficient SCALAPACK installation: "
          "The required symbols pdsyevr_ and pssyevr_ were not found."
          )
        SET(SCALAPACK_ADDITIONAL_ERROR_STRING
          ${SCALAPACK_ADDITIONAL_ERROR_STRING}
          "Could not find a sufficient SCALAPACK installation: \n"
          "SCALAPACK symbol check for pdsyevr_ and pssyevr_ failed! "
          "This usually means that your SCALAPACK installation is incomplete "
          "or the link line is broken. Consult\n"
          "  CMakeFiles/CMakeError.log\n"
          "for further information.\n"
          )
        SET(${var} FALSE)
      ENDIF()
    ENDIF()
  ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(SCALAPACK)
