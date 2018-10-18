## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
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

FOREACH(_dependency MPI CUDA)
  IF(NOT DEAL_II_WITH_${_dependency})
    IF(DEAL_II_WITH_CUDA_AWARE_MPI)
      MESSAGE(FATAL_ERROR "\n"
        "DEAL_II_WITH_CUDA_AWARE_MPI has unmet configuration requirements: "
        "DEAL_II_WITH_${_dependency} has to be set to \"ON\".\n\n"
        )
    ELSE()
      MESSAGE(STATUS
        "DEAL_II_WITH_CUDA_AWARE_MPI has unmet configuration requirements: "
        "DEAL_II_WITH_${_dependency} has to be set to \"ON\"."
        )
      SET(DEAL_II_WITH_CUDA_AWARE_MPI OFF)
    ENDIF()
  ENDIF()
ENDFOREACH()
