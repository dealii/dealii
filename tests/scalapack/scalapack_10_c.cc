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

#include "../lapack/create_matrix.h"
#include "../tests.h"

// test saving and loading of the State and Property of distributed ScaLAPACKMatrices

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/scalapack.h>

#include <cstdio>
#include <fstream>
#include <iostream>


template <typename NumberType>
void
test()
{
  const std::string filename("scalapack_10_test.h5");

  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));
  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  pcout << "Saving and restoring the state and property of ScaLAPACKMatrix"
        << std::endl;

  const unsigned int size = 100, block_size = 8;

  //create FullMatrix and fill it
  FullMatrix<NumberType> full(100);
  create_spd(full);

  //create 2d process grid
  std::shared_ptr<Utilities::MPI::ProcessGrid> grid =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size, block_size);

  ScaLAPACKMatrix<NumberType> scalapack_matrix(
    size, size, grid, block_size, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_matrix_copy(
    size, size, grid, block_size, block_size);

  scalapack_matrix.set_property(LAPACKSupport::Property::diagonal);
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_property() ==
                scalapack_matrix_copy.get_property(),
              ExcInternalError());

  scalapack_matrix.set_property(LAPACKSupport::Property::general);
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_property() ==
                scalapack_matrix_copy.get_property(),
              ExcInternalError());

  scalapack_matrix.set_property(LAPACKSupport::Property::hessenberg);
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_property() ==
                scalapack_matrix_copy.get_property(),
              ExcInternalError());

  scalapack_matrix.set_property(LAPACKSupport::Property::lower_triangular);
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_property() ==
                scalapack_matrix_copy.get_property(),
              ExcInternalError());

  scalapack_matrix.set_property(LAPACKSupport::Property::symmetric);
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_property() ==
                scalapack_matrix_copy.get_property(),
              ExcInternalError());

  scalapack_matrix.set_property(LAPACKSupport::Property::upper_triangular);
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_property() ==
                scalapack_matrix_copy.get_property(),
              ExcInternalError());

  // after construction the matrix state is LAPACKSupport::State::unusable
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_state() == scalapack_matrix_copy.get_state(),
              ExcInternalError());

  // the assignment operator changes the state to LAPACKSupport::State::matrix
  scalapack_matrix = full;
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_state() == scalapack_matrix_copy.get_state(),
              ExcInternalError());

  // calling invert changes the state to LAPACKSupport::inverse_matrix
  scalapack_matrix.set_property(LAPACKSupport::Property::symmetric);
  scalapack_matrix.invert();
  scalapack_matrix.save(filename.c_str());
  scalapack_matrix_copy.load(filename.c_str());
  std::remove(filename.c_str());
  AssertThrow(scalapack_matrix.get_state() == scalapack_matrix_copy.get_state(),
              ExcInternalError());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  test<double>();
}
