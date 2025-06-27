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



// document bug in the DerivativeApproximation in parallel
/*
 --------------------------------------------------------
 An error occurred in line <3349> of file
 </ssd/deal-trunk/deal.II/include/deal.II/dofs/dof_accessor.templates.h> in
 function void dealii::DoFCellAccessor<DoFHandlerType,
 lda>::get_dof_values(const InputVector&, ForwardIterator, ForwardIterator)
 const [with InputVector = dealii::TrilinosWrappers::MPI::Vector,
 ForwardIterator = double*, DoFHandlerType = dealii::DoFHandler<2>, bool
 level_dof_access = false] The violated condition was: this->is_artificial() ==
 false The name and call sequence of the exception was: ExcMessage ("Can't ask
 for DoF indices on artificial cells.") Additional Information: Can't ask for
 DoF indices on artificial cells.
 --------------------------------------------------------

 */

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/derivative_approximation.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  const IndexSet locally_relevant_set =
    DoFTools::extract_locally_relevant_dofs(dofh);

  TrilinosWrappers::MPI::Vector vec(dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  for (unsigned int i = vec.local_range().first; i < vec.local_range().second;
       ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  TrilinosWrappers::MPI::Vector vec_rel(locally_relevant_set);
  vec_rel = vec;

  MappingQ<dim> mapping(1);
  Vector<float> indicators(tr.n_active_cells());
  DerivativeApproximation::approximate_gradient(mapping,
                                                dofh,
                                                vec_rel,
                                                indicators);

  // we got here, so no exception.
  if (myid == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
