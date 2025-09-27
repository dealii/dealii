// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test the testsuite framework. Verify that we correctly fail when there is a
// PETSc object leak

#include <petscvec.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_context(argc, argv, 1);
  MPILogInitAll                    mpi_log;

  Vec x;
  VecCreateSeq(PETSC_COMM_SELF, 42, &x);

  deallog << "OK" << std::endl;
}
