// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
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
  test that swapping two TrilinosWrapper::MPI::Vector objects
  also swaps has_ghosts and locally_owned_elements.
 */


#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

void
print(TrilinosWrappers::MPI::Vector &v, unsigned int first_element)
{
  deallog << "size= " << v.size() << " el(" << first_element
          << ")= " << v(first_element)
          << " has_ghost_elements: " << v.has_ghost_elements() << std::endl
          << " locally_owned_elements: " << std::endl;
  v.locally_owned_elements().print(deallog);
}


void
test()
{
  unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector v;
  IndexSet                      owned_indices_1(6);
  IndexSet                      ghosted_indices_1(6);
  owned_indices_1.add_range(my_id * 3, (my_id + 1) * 3);
  ghosted_indices_1.add_range((1 - my_id) * 3, (2 - my_id) * 3);
  v.reinit(owned_indices_1, ghosted_indices_1, MPI_COMM_WORLD);

  TrilinosWrappers::MPI::Vector w;
  IndexSet                      owned_indices_2(10);
  owned_indices_2.add_range(my_id * 5, (my_id + 1) * 5);
  w.reinit(owned_indices_2, MPI_COMM_WORLD);
  for (unsigned int i = my_id * 5; i < (my_id + 1) * 5; ++i)
    w(i) = 2;


  deallog << "v: ";
  print(v, my_id * 3);
  deallog << "w: ";
  print(w, my_id * 5);

  deallog << "**swap**" << std::endl;

  swap(v, w);

  deallog << "v: ";
  print(v, my_id * 5);
  deallog << "w: ";
  print(w, my_id * 3);

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    mpi_log;

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
