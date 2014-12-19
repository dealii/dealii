// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Tests basic stuff of Trilinos sparsity patterns

#include "../tests.h"
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <fstream>
#include <iomanip>


void test ()
{
  TrilinosWrappers::SparsityPattern sp;
  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  deallog << "Creating entries..." << std::endl;

  sp.reinit(5,7,3);
  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<7; ++j)
      if ((i+2*j+1) % 3 == 0)
        sp.add (i,j);

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  sp.compress ();

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;
  deallog << "Number of entries: " << sp.n_nonzero_elements() << std::endl;
  deallog << "Number of rows: " << sp.n_rows() << std::endl;
  deallog << "Number of colums: " << sp.n_cols() << std::endl;
  deallog << "Local size: " << sp.local_size() << std::endl;
  deallog << "Max row length: " << sp.max_entries_per_row() << std::endl;
  deallog << "SP::row_length(0): " << sp.row_length(0) << std::endl;
  deallog << "Bandwidth: " << sp.bandwidth() << std::endl;
  deallog << "SP::empty(): " << sp.empty() << std::endl;

  sp.compress ();
  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;

  deallog << "Clearing..." << std::endl;

  sp.clear();

  deallog << "SP::is_compressed(): " << sp.is_compressed() << std::endl;
  deallog << "Bandwidth: " << sp.bandwidth() << std::endl;
  deallog << "SP::empty(): " << sp.empty() << std::endl;
  deallog << "Number of rows: " << sp.n_rows() << std::endl;

  deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
