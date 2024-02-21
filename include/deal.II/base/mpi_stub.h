// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mpi_stub_h
#define dealii_mpi_stub_h

#include <deal.II/base/config.h>

// If we have mpi.h then include it. Otherwise, define some common MPI data
// types and global constants for the no-MPI case. This way we can still use,
// e.g., MPI_Comm in the API.

#if defined(DEAL_II_WITH_MPI)
#  include <mpi.h>
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
