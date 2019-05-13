// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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



// check TrilinosScalar

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

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
    Assert(false, ExcNotImplemented());
}
