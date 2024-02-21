// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test the various VectorType::extract_subvector_to functions

#include <deal.II/base/index_set.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <typename VectorType>
void
test(VectorType &vector)
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  for (unsigned int i = 0; i < vector.size(); ++i)
    vector(i) = i;

  // select every other element
  std::vector<typename VectorType::size_type> indices;
  for (unsigned int j = 0; j < vector.size() / 2; ++j)
    indices.push_back(2 * j);

  // do the extraction with the function that takes indices, then
  // assert correctness
  std::vector<typename VectorType::value_type> values1(indices.size());
  vector.extract_subvector_to(indices, values1);
  for (unsigned int j = 0; j < vector.size() / 2; ++j)
    AssertThrow(get_real_assert_zero_imag(values1[j]) == 2 * j,
                ExcInternalError());

  // do the same with the version of the function that takes iterators
  std::vector<typename VectorType::value_type> values2(indices.size());
  vector.extract_subvector_to(indices.begin(), indices.end(), values2.begin());
  for (unsigned int j = 0; j < vector.size() / 2; ++j)
    AssertThrow(get_real_assert_zero_imag(values2[j]) == 2 * j,
                ExcInternalError());

  // done
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  {
    {
      deallog.push("deal.II");
      Vector<double> v(17);
      test(v);
      deallog.pop();
    }

    {
      deallog.push("PETSc");
      PETScWrappers::MPI::Vector v(MPI_COMM_SELF, 17, 17);
      test(v);
      deallog.pop();
    }

    {
      deallog.push("Trilinos");
      TrilinosWrappers::MPI::Vector v;
      v.reinit(complete_index_set(17), MPI_COMM_WORLD);
      test(v);
      deallog.pop();
    }


    {
      deallog.push("deal.II");
      BlockVector<double> v(3);
      v.block(0).reinit(7);
      v.block(1).reinit(5);
      v.block(2).reinit(3);
      v.collect_sizes();
      test(v);
      deallog.pop();
    }

    {
      deallog.push("PETSc");
      std::vector<PETScWrappers::MPI::BlockVector::size_type> sizes(3);
      sizes[0] = 7;
      sizes[1] = 5;
      sizes[2] = 3;
      PETScWrappers::MPI::BlockVector v(sizes, MPI_COMM_SELF, sizes);
      test(v);
      deallog.pop();
    }

    {
      deallog.push("Trilinos");
      TrilinosWrappers::MPI::BlockVector v(3);
      v.block(0).reinit(complete_index_set(7));
      v.block(1).reinit(complete_index_set(5));
      v.block(2).reinit(complete_index_set(3));
      v.collect_sizes();
      test(v);
      deallog.pop();
    }
  }
}
