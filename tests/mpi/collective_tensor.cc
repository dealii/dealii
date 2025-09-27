// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check Utilities::MPI::sum() for tensor objects

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"

void
test()
{
  Assert(Utilities::MPI::job_supports_mpi(), ExcInternalError());

  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  Tensor<1, 1> tensor_1;

  tensor_1[0] = 1.0;

  Tensor<1, 1> result_1 = Utilities::MPI::sum(tensor_1, MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "1D tensor: " << result_1 << std::endl;

  Tensor<1, 2> tensor_2;

  tensor_2[0] = 1.0;
  tensor_2[1] = 2.0;

  Tensor<1, 2> result_2 = Utilities::MPI::sum(tensor_2, MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "2D tensor: " << result_2 << std::endl;

  Tensor<1, 3> tensor_3;

  tensor_3[0] = 1.0;
  tensor_3[1] = 2.0;
  tensor_3[2] = 3.0;

  Tensor<1, 3> result_3 = Utilities::MPI::sum(tensor_3, MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "3D tensor: " << result_3 << std::endl;

  Tensor<2, 3> rank_2_tensor;

  rank_2_tensor[0][0] = 1.0;
  rank_2_tensor[0][1] = 2.0;
  rank_2_tensor[0][2] = 3.0;
  rank_2_tensor[1][0] = 4.0;
  rank_2_tensor[1][1] = 5.0;
  rank_2_tensor[1][2] = 6.0;
  rank_2_tensor[2][0] = 7.0;
  rank_2_tensor[2][1] = 8.0;
  rank_2_tensor[2][2] = 9.0;

  Tensor<2, 3> result_rank_2 =
    Utilities::MPI::sum(rank_2_tensor, MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "Rank 2 tensor: " << result_rank_2 << std::endl;

  SymmetricTensor<2, 2> symmetric_tensor;

  symmetric_tensor[0][0] = 1.0;
  symmetric_tensor[0][1] = 2.0;
  symmetric_tensor[1][1] = 3.0;

  SymmetricTensor<2, 2> result_symmetric =
    Utilities::MPI::sum(symmetric_tensor, MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "Symmetric tensor: " << result_symmetric << std::endl;
}


int
main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("mpi");
      test();
      deallog.pop();
    }
  else
    test();
}
