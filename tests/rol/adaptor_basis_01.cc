// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check ROLAdaptor::basis() on a distributed vector.

#include <deal.II/base/index_set.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/trilinos/rol_adaptor.h>

#include "../tests.h"


template <typename VectorType>
void
test()
{
  // set up distributed vector
  const unsigned int nprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int myid   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  IndexSet owned(nprocs);
  owned.add_index(myid);

  VectorType v(owned, MPI_COMM_WORLD);
  v[myid] = myid;

  // wrap for ROL
  TrilinosWrappers::ROLAdaptor<VectorType> rol_v(ROL::makePtrFromRef(v));

  // get all basis vectors
  for (unsigned int b = 0; b < nprocs; ++b)
    {
      const auto basis = rol_v.basis(b);
      basis->print(deallog.get_file_stream());
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<TrilinosWrappers::MPI::Vector>();
  test<LinearAlgebra::distributed::Vector<double>>();
}
