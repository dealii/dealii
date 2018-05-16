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



// check TrilinosWrappers::MPI::Vector::operator = (const ::Vector &)

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <vector>


void
test (TrilinosWrappers::MPI::Vector &v)
{
  // set only certain elements of the
  // vector.
  std::vector<bool> pattern (v.size(), false);
  for (unsigned int i=0; i<v.size(); i+=1+i)
    {
      v(i) += i;
      pattern[i] = true;
    }

  v.compress (VectorOperation::add);

  Vector<double> w(v.size());
  w=v;
  Vector<float>  x(v.size());
  x=v;

  TrilinosWrappers::MPI::Vector w1;
  w1.reinit(complete_index_set(v.size()), MPI_COMM_WORLD);
  w1=w;
  TrilinosWrappers::MPI::Vector x1;
  x1.reinit(complete_index_set(v.size()), MPI_COMM_WORLD);
  x1=x;

  for (unsigned int i=0; i<v.size(); ++i)
    {
      AssertThrow (w1(i) == w(i), ExcInternalError());
      AssertThrow (x1(i) == x(i), ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int
main (int argc,char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());


  try
    {
      {
        TrilinosWrappers::MPI::Vector v;
        v.reinit(complete_index_set(100), MPI_COMM_WORLD);
        test (v);
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
