// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check AffineConstraints<double>.distribute() for a petsc vector
//
// like _01, but with an inhomogeneity

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"



void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  // create a vector that consists of elements indexed from 0 to n
  PETScWrappers::MPI::Vector vec(MPI_COMM_WORLD, 100 * n_processes, 100);
  AssertThrow(vec.locally_owned_size() == 100, ExcInternalError());
  AssertThrow(vec.local_range().first == 100 * myid, ExcInternalError());
  AssertThrow(vec.local_range().second == 100 * myid + 100, ExcInternalError());
  for (unsigned int i = vec.local_range().first; i < vec.local_range().second;
       ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  // verify correctness so far
  {
    double exact_l1 = 0;
    for (unsigned int i = 0; i < vec.size(); ++i)
      exact_l1 += i;
    AssertThrow(vec.l1_norm() == exact_l1, ExcInternalError());
  }


  // create a AffineConstraints<double> with a range that exceeds the locally
  // owned range by 50 on each side
  IndexSet locally_relevant_range(vec.size());
  locally_relevant_range.add_range(
    std::max<int>(100 * myid - 50, 0),
    std::min(static_cast<types::global_dof_index>(100 * myid + 150),
             vec.size()));
  AffineConstraints<PetscScalar> cm(vec.locally_owned_elements(),
                                    locally_relevant_range);

  // add constraints that constrain an element in the middle of the
  // local range of each processor against an element outside, both in
  // the ghost range before and after
  //
  // note that we tell each processor about all constraints, but most
  // of them will throw away this information since it is not for a
  // DoF inside the locally relevant range
  for (unsigned int p = 0; p < n_processes; ++p)
    {
      if ((p != 0) && locally_relevant_range.is_element(p * 100 + 10))
        {
          cm.add_line(p * 100 + 10);
          cm.add_entry(p * 100 + 10, p * 100 - 25, 1);
          cm.set_inhomogeneity(p * 100 + 10, 1);
        }

      if ((p != n_processes - 1) &&
          locally_relevant_range.is_element(p * 100 + 90))
        {
          cm.add_line(p * 100 + 90);
          cm.add_entry(p * 100 + 90, p * 100 + 105, 1);
          cm.set_inhomogeneity(p * 100 + 90, -1);
        }
    }
  cm.close();

  // now distribute these constraints
  cm.distribute(vec);

  // verify correctness
  if (myid != 0)
    AssertThrow(get_real_assert_zero_imag(vec(vec.local_range().first + 10)) ==
                  vec.local_range().first - 25 + 1,
                ExcInternalError());

  if (myid != n_processes - 1)
    AssertThrow(get_real_assert_zero_imag(vec(vec.local_range().first + 90)) ==
                  vec.local_range().first + 105 - 1,
                ExcInternalError());

  for (unsigned int i = vec.local_range().first; i < vec.local_range().second;
       ++i)
    {
      if ((i != vec.local_range().first + 10) &&
          (i != vec.local_range().first + 90))
        {
          PetscScalar val = vec(i);
          AssertThrow(std::fabs(get_real_assert_zero_imag(val) - i) <= 1e-6,
                      ExcInternalError());
        }
    }

  {
    double exact_l1 = 0;

    // add up original values of vector entries
    for (unsigned int i = 0; i < vec.size(); ++i)
      exact_l1 += i;

    // but then correct for the constrained values
    for (unsigned int p = 0; p < n_processes; ++p)
      {
        if (p != 0)
          exact_l1 = exact_l1 - (p * 100 + 10) + (p * 100 - 25);
        if (p != n_processes - 1)
          exact_l1 = exact_l1 - (p * 100 + 90) + (p * 100 + 105);
      }

    const double l1_norm = vec.l1_norm();
    AssertThrow(l1_norm == exact_l1, ExcInternalError());

    if (myid == 0)
      deallog << "Norm = " << l1_norm << std::endl;
  }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      test();
    }
  else
    test();
}
