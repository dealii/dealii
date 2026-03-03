## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2026 by the deal.II authors
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
# Configuration for the t8code and sc libraries:
#

set(FEATURE_T8CODE_DEPENDS MPI)


macro(feature_t8code_find_external var)
  message(STATUS ${T8CODE_DIR})
  find_package(DEAL_II_T8CODE)
  if(T8CODE_FOUND)
    set(${var} TRUE)

#
# Check whether t8code supports mpi:
#
if(NOT T8CODE_ENABLE_MPI)
  message(STATUS "Insufficient t8code installation found: "
    "t8code has to be configured with MPI enabled."
    )
  set(T8CODE_ADDITIONAL_ERROR_STRING
    ${T8CODE_ADDITIONAL_ERROR_STRING}
    "Insufficient t8code installation found!\n"
    "t8code has to be configured with MPI enabled.\n"
    )
  set(${var} FALSE)
endif()

    check_mpi_interface(T8CODE ${var})
  endif()
endmacro()


configure_feature(T8CODE)
