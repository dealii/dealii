// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/psblas_vector.h>

#include <psb_c_dbase.h>

#include <iostream>

#include "../../tests/tests.h"

using namespace dealii;

// Test for swap(), add(), sadd(), and add_and_dot(), equ(), and so on for
// PSBLAS vectors

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPI_Comm                         mpi_communicator = MPI_COMM_WORLD;
  AssertThrow(Utilities::MPI::n_mpi_processes(mpi_communicator) == 2,
              ExcMessage("This test needs to be run with 2 MPI processes."));

  MPILogInitAll log;

  int id;
  MPI_Comm_rank(mpi_communicator, &id);
  IndexSet locally_owned_dofs(25);
  if (id == 0)
    locally_owned_dofs.add_range(0, 15);
  else if (id == 1)
    locally_owned_dofs.add_range(15, 25);


  IndexSet locally_relevant_dofs(25);
  locally_relevant_dofs = locally_owned_dofs;
  if (id == 0)
    locally_relevant_dofs.add_range(15, 17);
  else if (id == 1)
    locally_relevant_dofs.add_range(12, 15);

  PSCToolkitWrappers::Vector x(locally_owned_dofs, mpi_communicator);

  for (const types::global_dof_index idx : locally_owned_dofs)
    x(idx) += idx;
  x.compress(VectorOperation::add);

  PSCToolkitWrappers::Vector y;
  y.reinit(locally_owned_dofs, mpi_communicator);
  for (const types::global_dof_index idx : locally_owned_dofs)
    y(idx) += idx + 0.5;
  y.compress(VectorOperation::add);

  // Test swap function
  double old_value = x(*locally_owned_dofs.begin());
  x.swap(y);
  AssertThrow(y(*locally_owned_dofs.begin()) == old_value,
              ExcMessage("Test for swap() failed."));

  // Test add()
  x.add(2.0, y);
  for (types::global_dof_index idx : locally_owned_dofs)
    {
      deallog << "x(" << idx << ")=" << x(idx) << std::endl;
      AssertThrow(x(idx) == 3.0 * idx + 0.5,
                  ExcMessage("Test for add() failed."));
    }

  // Test sadd()
  x.reinit(locally_owned_dofs, mpi_communicator);
  y.reinit(locally_owned_dofs, mpi_communicator);
  for (const types::global_dof_index idx : locally_owned_dofs)
    {
      x(idx) += idx;
      y(idx) += idx + 0.5;
    }
  x.compress(VectorOperation::add);
  y.compress(VectorOperation::add);

  x.sadd(3.0, y);
  for (types::global_dof_index idx : locally_owned_dofs)
    {
      deallog << "x(" << idx << ")=" << x(idx) << std::endl;
      AssertThrow(x(idx) == 3.0 * idx + y(idx),
                  ExcMessage("Test for sadd() failed."));
    }

  // Test add_and_dot()
  PSCToolkitWrappers::Vector w;
  x.reinit(locally_owned_dofs, mpi_communicator);
  y.reinit(locally_owned_dofs, mpi_communicator);
  w.reinit(locally_owned_dofs, mpi_communicator);
  for (const types::global_dof_index idx : locally_owned_dofs)
    {
      x(idx) += idx;
      y(idx) += 1.;
      w(idx) += idx + 1;
    }
  x.compress(VectorOperation::add);
  y.compress(VectorOperation::add);
  w.compress(VectorOperation::add);

  const double dot_product = x.add_and_dot(1.0, y, w);
  Assert(dot_product - std::pow(w.l2_norm(), 2) < 1e-11,
         ExcMessage("Test for add_and_dot() failed."));

  //  Test equ(), while there we also use operator[]
  PSCToolkitWrappers::Vector v;
  v.reinit(locally_owned_dofs, mpi_communicator);
  v.equ(1.0, x);
  for (const types::global_dof_index idx : locally_owned_dofs)
    AssertThrow(v[idx] == x[idx], ExcMessage("Test for equ() failed."));

  deallog << "Ok" << std::endl;

  // Test reinit(const Vector &v, const bool omit_zeroing_entries=false)
  PSCToolkitWrappers::Vector z;
  z.reinit(x, true);
  for (const types::global_dof_index idx : locally_owned_dofs)
    AssertThrow(z[idx] == 0.0, ExcMessage("Entry not zero."));

  // Test all_zeros()
  AssertThrow(z.all_zero(), ExcMessage("Vector must be all zeros."));

  // Test add(scalar)
  z.add(2.0);
  for (const types::global_dof_index idx : locally_owned_dofs)
    AssertThrow(z[idx] == 2.0, ExcMessage("Entry not correct after add()."));

  // Test scale(V)
  PSCToolkitWrappers::Vector scale_vector;
  scale_vector.reinit(locally_owned_dofs, mpi_communicator);
  scale_vector = 2.0;
  z.scale(scale_vector);
  for (const types::global_dof_index idx : locally_owned_dofs)
    AssertThrow(z[idx] == 4.0, ExcMessage("Entry not correct after scale()."));

  // Test mean_value()
  AssertThrow(z.mean_value() == 4.0, ExcMessage("Mean value not correct."));

  // Test add(V)
  z = 0.;
  z += v;
  for (const types::global_dof_index idx : locally_owned_dofs)
    AssertThrow(z[idx] == v[idx], ExcMessage("Add vector failed."));

  return 0;
}
