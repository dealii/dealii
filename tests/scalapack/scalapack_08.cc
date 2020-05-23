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

#include "../tests.h"

#include "../lapack/create_matrix.h"

// test compute_SVD(ScaLAPACKMatrix<NumberType>*,ScaLAPACKMatrix<NumberType>*)

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/process_grid.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/scalapack.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>


template <typename NumberType>
void
test(const unsigned int size,
     const unsigned int block_size,
     const NumberType   tol)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(mpi_communicator));
  const unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  ConditionalOStream pcout(std::cout, (this_mpi_process == 0));

  std::shared_ptr<Utilities::MPI::ProcessGrid> grid_2d =
    std::make_shared<Utilities::MPI::ProcessGrid>(
      mpi_communicator, size, size, block_size, block_size);

  pcout << size << " " << block_size << std::endl;

  // Create s.p.d matrices of requested size:
  FullMatrix<NumberType> full_A(size);
  create_spd(full_A);

  // compute eigenpairs of s.p.d matrix
  ScaLAPACKMatrix<NumberType> scalapack_A_ev(size, grid_2d, block_size);
  scalapack_A_ev.set_property(LAPACKSupport::Property::symmetric);
  scalapack_A_ev = full_A;
  std::vector<NumberType> eigenvalues =
    scalapack_A_ev.eigenpairs_symmetric_by_index(std::make_pair(0, size - 1),
                                                 true);
  FullMatrix<NumberType> eigenvectors(size, size);
  scalapack_A_ev.copy_to(eigenvectors);

  // compute SVD of s.p.d matrix A = U * SIGMA * VT
  ScaLAPACKMatrix<NumberType> scalapack_A_sv(size, grid_2d, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_U(size, grid_2d, block_size);
  ScaLAPACKMatrix<NumberType> scalapack_VT(size, grid_2d, block_size);
  scalapack_A_sv.set_property(LAPACKSupport::Property::symmetric);
  scalapack_A_sv = full_A;
  std::vector<NumberType> singular_values =
    scalapack_A_sv.compute_SVD(&scalapack_U, &scalapack_VT);
  FullMatrix<NumberType> l_singular_vectors(size, size);
  FullMatrix<NumberType> r_singular_vectors(size, size);
  scalapack_U.copy_to(l_singular_vectors);
  scalapack_VT.copy_to(r_singular_vectors);

  const unsigned int max_num_values = 5;
  pcout << "comparing the SVD and Eigendecomposition of a s.p.d matrix"
        << std::endl;
  for (unsigned i = 0; i < max_num_values; ++i)
    AssertThrow(std::abs(eigenvalues[size - 1 - i] - singular_values[i]) < tol,
                ExcMessage("singular and eigenvalues do not match"));
  pcout
    << "   with respect to the given tolerance the singular and eigenvalues coincide"
    << std::endl;

  Vector<NumberType> eigenvector(size), l_singular_vector(size),
    r_singular_vector(size);
  for (unsigned int i = 0; i < max_num_values; ++i)
    {
      for (unsigned int j = 0; j < size; ++j)
        {
          eigenvector[j]       = eigenvectors(j, size - 1 - i);
          l_singular_vector[j] = l_singular_vectors(j, i);
          r_singular_vector[j] = r_singular_vectors(i, j);
        }
      NumberType product_1 = eigenvector * l_singular_vector;
      NumberType product_2 = eigenvector * r_singular_vector;
      // the tolerance is reduced for the singular vectors
      AssertThrow((std::abs(product_1) - 1) < tol * 10,
                  ExcMessage(
                    "left singular vectors and eigenvectors do not coincide"));
      AssertThrow((std::abs(product_2) - 1) < tol * 10,
                  ExcMessage(
                    "right singular vectors and eigenvectors do not coincide"));
    }
  pcout
    << "   with respect to the given tolerance the right and left singular vectors coincide with the eigenvectors"
    << std::endl;
  pcout << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  const std::vector<unsigned int> sizes      = {{200, 400, 600}};
  const std::vector<unsigned int> blocks     = {{32, 64}};
  const double                    tol_double = 1e-10;

  for (const auto &s : sizes)
    for (const auto &b : blocks)
      if (b <= s)
        test<double>(s, b, tol_double);
}
