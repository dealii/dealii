// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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



// check VectorTools::subtract_mean_value() for deal.II parallel vector

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

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
      IndexSet ghost_indices(10);
      ghost_indices.add_range(3, 8);
      {
        LinearAlgebra::distributed::Vector<double> v(local_range,
                                                     ghost_indices,
                                                     MPI_COMM_WORLD);
        test(v);
      }

      {
        LinearAlgebra::distributed::Vector<float> v(local_range,
                                                    ghost_indices,
                                                    MPI_COMM_WORLD);
        test(v);
      }

      {
        LinearAlgebra::distributed::BlockVector<double> v(
          std::vector<IndexSet>(1, local_range),
          std::vector<IndexSet>(1, ghost_indices),
          MPI_COMM_WORLD);
        test(v);
      }

      {
        LinearAlgebra::distributed::BlockVector<float> v(
          std::vector<IndexSet>(1, local_range),
          std::vector<IndexSet>(1, ghost_indices),
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
