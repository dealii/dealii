// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*
  test ::Vector::swap() (fixed in r 25668)
 */


#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

void
print(LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> &v)
{
  deallog << "size= " << v.size() << " el(0)= " << v(0)
          << " l2norm()= " << v.l2_norm() << std::endl;
}


void
test()
{
  LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> v;
  v.reinit(complete_index_set(5), MPI_COMM_WORLD);
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = 1;
  LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default> w;
  w.reinit(complete_index_set(9), MPI_COMM_WORLD);
  for (unsigned int i = 0; i < w.size(); ++i)
    w(i) = 2;


  deallog << "v: ";
  print(v);
  deallog << "w: ";
  print(w);

  deallog << "**swap**" << std::endl;

  swap(v, w);

  deallog << "v: ";
  print(v);
  deallog << "w: ";
  print(w);

  Assert(v.size() == 9, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());


  try
    {
      test();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
