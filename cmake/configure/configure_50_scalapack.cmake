## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2020 by the deal.II authors
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

set(FEATURE_SCALAPACK_DEPENDS MPI LAPACK)


macro(feature_scalapack_find_external var)
  find_package(DEAL_II_SCALAPACK)

  if(SCALAPACK_FOUND)
    set(${var} TRUE)
    check_mpi_interface(SCALAPACK ${var})

    if (${var})
      clear_cmake_required()
      set(CMAKE_REQUIRED_LIBRARIES ${SCALAPACK_LIBRARIES} ${LAPACK_LIBRARIES})
      CHECK_C_SOURCE_COMPILES("
        void pdsyevr_();
        void pssyevr_();
        int main(){
          pdsyevr_();
          pssyevr_();
          return 0;
        }"
        DEAL_II_SCALAPACK_HAS_PDSYEVR_PSSYEVR)
        reset_cmake_required()

      if(NOT DEAL_II_SCALAPACK_HAS_PDSYEVR_PSSYEVR)
        message(STATUS "Could not find a sufficient SCALAPACK installation: "
          "The required symbols pdsyevr_ and pssyevr_ were not found."
          )
        set(SCALAPACK_ADDITIONAL_ERROR_STRING
          ${SCALAPACK_ADDITIONAL_ERROR_STRING}
          "Could not find a sufficient SCALAPACK installation: \n"
          "SCALAPACK symbol check for pdsyevr_ and pssyevr_ failed! "
          "This usually means that your SCALAPACK installation is incomplete "
          "or the link line is broken. Consult\n"
          "  CMakeFiles/CMakeError.log\n"
          "for further information.\n"
          )
        set(${var} FALSE)
      endif()
    endif()
  endif()
endmacro()

configure_feature(SCALAPACK)
