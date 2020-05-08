// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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



// check VectorTools::subtract_mean_value() for Trilinos vectors

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <vector>

#include "../tests.h"

template <class VectorType>
void
test(VectorType &v)
{
  // set some elements of the vector
  unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  for (unsigned int i = 5 * my_id; i < 5 * (my_id + 1); ++i)
    {
      v(i) = i;
    }
  v.compress(VectorOperation::insert);

  // then check the norm
  VectorTools::subtract_mean_value(v);
  AssertThrow(std::fabs(v.mean_value()) < 1e-10 * v.l2_norm(),
              ExcInternalError());

  deallog << "OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  try
    {
      IndexSet local_range(10);
      local_range.add_range(5 * my_id, 5 * (my_id + 1));
      {
        TrilinosWrappers::MPI::Vector v1(local_range, MPI_COMM_WORLD);
        test(v1);

        LinearAlgebra::ReadWriteVector<double> v_tmp(local_range);
        LinearAlgebra::EpetraWrappers::Vector  v2(local_range, MPI_COMM_WORLD);
        v_tmp.import(v1, VectorOperation::insert);
        v2.import(v_tmp, VectorOperation::insert);
        VectorTools::subtract_mean_value(v2);
        AssertThrow(std::fabs(v2.mean_value()) < 1e-10 * v2.l2_norm(),
                    ExcInternalError());
        deallog << "OK" << std::endl;
      }
      {
        TrilinosWrappers::MPI::BlockVector v(std::vector<IndexSet>(1,
                                                                   local_range),
                                             MPI_COMM_WORLD);
        test(v);
      }
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
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
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
