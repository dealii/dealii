// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#include "../tests.h"

// test copying FullMatrix into distributed ScaLAPACK matrix and back.

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/multithread_info.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#include <deal.II/lac/vector.h>

#include <deal.II/lac/scalapack.h>

#include <fstream>
#include <iostream>

template <typename NumberType>
void test(const unsigned int size, const unsigned int block_size)
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout (std::cout, (this_mpi_process ==0));

  // test multiplication with random vectors
  boost::random::mt19937 gen;
  boost::random::uniform_01<> dist;

  // Create SPD matrices of requested size:
  FullMatrix<NumberType> full_in(size), full_out(size), diff(size);

  ScaLAPACKMatrix<NumberType> scalapack_matrix(size,
                                               mpi_communicator,
                                               block_size);

  pcout << size << " " << block_size << " " << scalapack_matrix.get_process_grid_rows() << " " << scalapack_matrix.get_process_grid_columns() << std::endl;

  {
    unsigned int index = 0;
    for (unsigned int i = 0; i < size; ++i)
      for (unsigned int j = 0; j < size; ++j)
        full_in(i,j) = index++;
  }

  scalapack_matrix = full_in;
  scalapack_matrix.copy_to(full_out);

  diff = 0;
  diff.add( 1., full_in);
  diff.add(-1., full_out);

  const double error = diff.frobenius_norm();

  if (error > 1e-10 && this_mpi_process == 0)
    {
      std::cout << "Error!" << std::endl
                << "expected to have:" << std::endl;
      full_in.print_formatted (std::cout);
      std::cout << "but got:" << std::endl;
      full_out.print_formatted(std::cout);
      AssertThrow (false, dealii::ExcInternalError());
    }
  else
    pcout << "Ok" << std::endl;
}



int main (int argc,char **argv)
{

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

      const std::vector<unsigned int> sizes = {{64,120,320,640}};
      const std::vector<unsigned int> blocks = {{32,64}};

      for (const auto &s : sizes)
        for (const auto &b : blocks)
          test<double>(s,b);

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
