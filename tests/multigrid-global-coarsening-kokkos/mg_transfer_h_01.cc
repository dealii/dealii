// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


/**
 * Test transfer operator for global coarsening on a uniformly refined meshes.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/portable_mg_transfer_global_coarsening.h>

#include "mg_transfer_util.h"

using namespace dealii;

template <int dim, typename Number>
void
do_test(const FiniteElement<dim> &fe_fine, const FiniteElement<dim> &fe_coarse)
{
  // create coarse grid
  parallel::distributed::Triangulation<dim> tria_coarse(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_coarse);
  tria_coarse.refine_global(1);

  // create_fine_grid(tria_coarse);
  parallel::distributed::Triangulation<dim> tria_fine(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_fine);
  tria_fine.refine_global(2);

  // setup dof-handlers
  DoFHandler<dim> dof_handler_fine(tria_fine);
  dof_handler_fine.distribute_dofs(fe_fine);

  DoFHandler<dim> dof_handler_coarse(tria_coarse);
  dof_handler_coarse.distribute_dofs(fe_coarse);

  AffineConstraints<Number> constraint_coarse;
  AffineConstraints<Number> constraint_fine;

  // setup transfer operator
  Portable::MGTwoLevelTransfer<
    dim,
    LinearAlgebra::distributed::Vector<Number, MemorySpace::Default>>
    transfer;
  transfer.reinit_geometric_transfer(dof_handler_fine,
                                     dof_handler_coarse,
                                     constraint_fine,
                                     constraint_coarse);

  deallog << "test first time" << std::endl;
  test_transfer_operator(transfer, dof_handler_fine, dof_handler_coarse);

  deallog << "test second time" << std::endl;
  test_transfer_operator(transfer, dof_handler_fine, dof_handler_coarse);
}

template <int dim, typename Number>
void
test(int fe_degree)
{
  const auto str_fine   = std::to_string(fe_degree);
  const auto str_coarse = std::to_string(fe_degree);

  if (fe_degree > 0)
    {
      deallog.push("CG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
      do_test<dim, double>(FE_Q<dim>(fe_degree), FE_Q<dim>(fe_degree));
      deallog.pop();
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  for (unsigned int i = 0; i < 5; ++i)
    test<2, double>(i);
}
