// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
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

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <mpi.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#else

// Without MPI, we would still like to use some constructs with MPI
// data types. Therefore, create some dummies. Since we only ever use
// them inside our own constructs, the right thing to do is to put
// them into namespace dealii:
DEAL_II_NAMESPACE_OPEN

using MPI_Comm     = int;
using MPI_Request  = int;
using MPI_Datatype = int;
using MPI_Op       = int;

constexpr MPI_Comm    MPI_COMM_WORLD   = 0;
constexpr MPI_Comm    MPI_COMM_SELF    = 0;
constexpr MPI_Comm    MPI_COMM_NULL    = 0;
constexpr MPI_Request MPI_REQUEST_NULL = 0;
constexpr MPI_Op      MPI_MIN          = 0;
constexpr MPI_Op      MPI_MAX          = 0;
constexpr MPI_Op      MPI_SUM          = 0;
constexpr MPI_Op      MPI_LOR          = 0;

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
