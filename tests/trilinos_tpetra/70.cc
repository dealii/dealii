// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check TrilinosScalar

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());


  if (typeid(TrilinosScalar) == typeid(double))
    deallog << "double" << std::endl;
  else if (typeid(TrilinosScalar) == typeid(float))
    deallog << "float" << std::endl;
  else
    DEAL_II_NOT_IMPLEMENTED();
}
