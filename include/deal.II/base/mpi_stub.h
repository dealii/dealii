// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_mpi_stub_h
#define dealii_mpi_stub_h

#include <deal.II/base/config.h>

// If we have mpi.h then include it. Otherwise, define some common MPI data
// types and global constants for the no-MPI case. This way we can still use,
// e.g., MPI_Comm in the API.

#if defined(DEAL_II_WITH_MPI) || defined(DEAL_II_WITH_PETSC)
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <mpi.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#else
// without MPI, we would still like to use
// some constructs with MPI data
// types. Therefore, create some dummies
using MPI_Comm     = int;
using MPI_Request  = int;
using MPI_Datatype = int;
using MPI_Op       = int;
#  ifndef MPI_COMM_WORLD
#    define MPI_COMM_WORLD 0
#  endif
#  ifndef MPI_COMM_SELF
#    define MPI_COMM_SELF 0
#  endif
#  ifndef MPI_COMM_NULL
#    define MPI_COMM_NULL 0
#  endif
#  ifndef MPI_REQUEST_NULL
#    define MPI_REQUEST_NULL 0
#  endif
#  ifndef MPI_MIN
#    define MPI_MIN 0
#  endif
#  ifndef MPI_MAX
#    define MPI_MAX 0
#  endif
#  ifndef MPI_SUM
#    define MPI_SUM 0
#  endif
#  ifndef MPI_LOR
#    define MPI_LOR 0
#  endif
#endif

#endif
