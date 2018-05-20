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

// check TrilinosWrappers::MPI::Vector::operator = (Vector)

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_vector.h>
#include <iostream>
#include <vector>

void
test(TrilinosWrappers::MPI::Vector& v)
{
  // set some entries of the vector
  for(unsigned int i = 0; i < v.size(); ++i)
    if(i % 3 == 0)
      v(i) = i + 1.;
  v.compress(VectorOperation::insert);

  // then copy it
  TrilinosWrappers::MPI::Vector w;
  w.reinit(complete_index_set(v.size()), MPI_COMM_WORLD);
  w = v;

  // make sure they're equal
  deallog << v * w << ' ' << v.l2_norm() * w.l2_norm() << ' '
          << v * w - v.l2_norm() * w.l2_norm() << std::endl;
  const double eps = typeid(TrilinosScalar) == typeid(double) ? 1e-14 : 1e-5;
  Assert(std::fabs(v * w - v.l2_norm() * w.l2_norm()) < eps * (v * w),
         ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main(int argc, char** argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      {
        TrilinosWrappers::MPI::Vector v;
        v.reinit(complete_index_set(100), MPI_COMM_WORLD);
        test(v);
      }
    }
  catch(std::exception& exc)
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
  catch(...)
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
