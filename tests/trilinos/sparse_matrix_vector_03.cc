// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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



// check SparseMatrix::vmult_add

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <iostream>
#include <vector>


void
test (TrilinosWrappers::MPI::Vector &v,
      TrilinosWrappers::MPI::Vector &w)
{
  TrilinosWrappers::SparseMatrix m(v.size(),v.size(),v.size());
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      m.set (i,j, i+2*j);

  for (unsigned int i=0; i<v.size(); ++i)
    {
      v(i) = i;
      w(i) = i;
    }

  m.compress (VectorOperation::insert);
  v.compress (VectorOperation::insert);
  w.compress (VectorOperation::insert);

  // w:=Mv
  m.vmult_add (w,v);

  // make sure we get the expected result
  for (unsigned int i=0; i<v.size(); ++i)
    {
      AssertThrow (v(i) == i, ExcInternalError());

      double result = 0;
      for (unsigned int j=0; j<m.m(); ++j)
        result += (i+2*j)*j;
      AssertThrow (w(i) == i+result, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int
main (int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());


  try
    {
      {
        TrilinosWrappers::MPI::Vector v;
        v.reinit(complete_index_set(100), MPI_COMM_WORLD);
        TrilinosWrappers::MPI::Vector w;
        w.reinit(complete_index_set(100), MPI_COMM_WORLD);
        test (v,w);
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
