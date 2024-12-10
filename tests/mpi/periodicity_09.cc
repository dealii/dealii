// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// The program attempts to create periodicity constraints when DoFHandler
// object carries FE_RaviartThomas element which does not provide
// get_subface_interpolation() method. See also pull-request
// Periodicity constraints: skip artificial cell face dofs #17918

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparsity_tools.h>

#include "../tests.h"

using namespace dealii;

// define boundaries for convenience
namespace Boundaries
{
  constexpr types::boundary_id lower_x = 0;
  constexpr types::boundary_id upper_x = 1;
  constexpr types::boundary_id lower_y = 2;
  constexpr types::boundary_id upper_y = 3;
  constexpr types::boundary_id lower_z = 4;
  constexpr types::boundary_id upper_z = 5;
} // namespace Boundaries


int
main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  mpi_initlog();

  constexpr unsigned int dim           = 2;
  constexpr unsigned int degree        = 1;
  unsigned int           n_refinements = 3;

  using namespace dealii;

  using Number = double;

  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  ConditionalOStream pcout(
    std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0));

  parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);

  // setup grid, in this case rectangle as example
  {
    Point<dim, Number> lower;
    Point<dim, Number> upper;
    for (unsigned int i = 0; i < dim; ++i)
      {
        lower[i] = -1.0;
        upper[i] = 1.0;
      }
    GridGenerator::hyper_rectangle(triangulation, lower, upper, true);
  }

  // add periodic boundary information to the triangulation
  {
    std::vector<GridTools::PeriodicFacePair<
      typename parallel::distributed::Triangulation<dim>::cell_iterator>>
      periodicity_vector;

    // identify face pairs in x-direction
    GridTools::collect_periodic_faces(triangulation,
                                      Boundaries::lower_x,
                                      Boundaries::upper_x,
                                      0,
                                      periodicity_vector);

    // identify face pairs in y-direction
    GridTools::collect_periodic_faces(triangulation,
                                      Boundaries::lower_y,
                                      Boundaries::upper_y,
                                      1,
                                      periodicity_vector);

    if (dim == 3)
      {
        // identify face pairs in z-direction
        GridTools::collect_periodic_faces(triangulation,
                                          Boundaries::lower_z,
                                          Boundaries::upper_z,
                                          2,
                                          periodicity_vector);
      }
    triangulation.add_periodicity(periodicity_vector);

    triangulation.refine_global(n_refinements);
  }

  DoFHandler<dim> dof_handler(triangulation);

  FESystem<dim, dim> fe(FE_RaviartThomas<dim>(degree), 1);

  dof_handler.distribute_dofs(fe);

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  AffineConstraints<Number> constraints;

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);


  // setup periodic boundary conditions
  {
    std::vector<
      GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
      matched_pairs;

    // identify face pairs in x-direction
    GridTools::collect_periodic_faces(
      dof_handler, Boundaries::lower_x, Boundaries::upper_x, 0, matched_pairs);

    // identify face pairs in y-direction
    GridTools::collect_periodic_faces(
      dof_handler, Boundaries::lower_y, Boundaries::upper_y, 1, matched_pairs);


    // identify face pairs in z-direction
    if (dim == 3)
      {
        GridTools::collect_periodic_faces(dof_handler,
                                          Boundaries::lower_z,
                                          Boundaries::upper_z,
                                          2,
                                          matched_pairs);
      }

    deallog << "Attempting to make constraints" << std::endl;

    DoFTools::make_periodicity_constraints<dim, dim>(matched_pairs,
                                                     constraints);
    constraints.close();
  }

  deallog << "The test completed successfully!" << std::endl;
}
