// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2022 by the deal.II authors
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


#include <deal.II/lac/trilinos_tpetra_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include "../tests.h"

#define PRINTME(name, var)                        \
  deallog << "Vector: " name << ':' << std::endl; \
  for (unsigned int i = 0; i < var.size(); ++i)   \
    deallog << var[i] << ' ';                     \
  deallog << std::endl;

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  {
    IndexSet local_owned(10);
    local_owned.add_range(0, 10);

    LinearAlgebra::TpetraWrappers::Vector<double> temp(local_owned,
                                                       MPI_COMM_WORLD);
    for (unsigned int i = 0; i < temp.size(); ++i)
      temp[i] = (double)(i + 1);

    PRINTME("Vector", temp);

    LinearAlgebra::TpetraWrappers::Vector<double> u(std::move(temp));
    PRINTME("move constructor", u);

    LinearAlgebra::TpetraWrappers::Vector<double> v;
    v = u;
    PRINTME("copy assignment", v);
    PRINTME("old object", u);

    v.clear();
    v = std::move(u);
    PRINTME("move assignment", v);
    deallog << "old object size: " << u.size() << std::endl;
  }
}
