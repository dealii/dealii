## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
      CHECK_CXX_SOURCE_COMPILES("
        extern \"C\" void pdsyevr_();
        extern \"C\" void pssyevr_();
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
