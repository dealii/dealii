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


// Test that a shared triangulation with artificial cells has the neighbors
// across periodic faces in its ghost layer immediately after calling
// add_periodicity(), i.e., without any subsequent mesh refinement.

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"

template <int dim>
void
test(const typename parallel::shared::Triangulation<dim>::Settings settings)
{
  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD,
                                            Triangulation<dim>::none,
                                            true,
                                            settings);

  std::vector<unsigned int> repetitions(dim, 8);
  Point<dim>                top_right;
  for (unsigned int d = 0; d < dim; ++d)
    top_right[d] = 1.;
  GridGenerator::subdivided_hyper_rectangle(
    tria, repetitions, Point<dim>(), top_right, true);

  const unsigned int n_owned = tria.n_locally_owned_active_cells();

  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    face_pairs;
  for (unsigned int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(tria, 2 * d, 2 * d + 1, d, face_pairs);
  tria.add_periodicity(face_pairs);

  deallog << "partition unchanged: "
          << (n_owned == tria.n_locally_owned_active_cells()) << std::endl;

  unsigned int n_artificial_partners = 0;
  for (const auto &[left, value] : tria.get_periodic_face_map())
    if (left.first->is_locally_owned() && value.first.first->is_artificial())
      ++n_artificial_partners;
  deallog << "artificial periodic partners: " << n_artificial_partners
          << std::endl;

  unsigned int n_ghosts = 0;
  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_ghost())
      ++n_ghosts;
  deallog << "ghost cells: " << n_ghosts << std::endl;

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FE_Q<dim>(1));

  AffineConstraints<double> constraints;
  constraints.reinit(dof_handler.locally_owned_dofs(),
                     DoFTools::extract_locally_relevant_dofs(dof_handler));
  for (unsigned int d = 0; d < dim; ++d)
    DoFTools::make_periodicity_constraints(
      dof_handler, 2 * d, 2 * d + 1, d, constraints);
  constraints.close();

  const auto   owned               = dof_handler.locally_owned_dofs();
  unsigned int n_owned_constraints = 0;
  for (const auto dof : owned)
    if (constraints.is_constrained(dof))
      ++n_owned_constraints;

  deallog << "constrained dofs: "
          << Utilities::MPI::sum(n_owned_constraints, MPI_COMM_WORLD)
          << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  using Tria1 = parallel::shared::Triangulation<1>;
  using Tria2 = parallel::shared::Triangulation<2>;
  using Tria3 = parallel::shared::Triangulation<3>;

  test<1>(Tria1::partition_zorder);
  test<2>(Tria2::partition_zorder);
  test<3>(Tria3::partition_zorder);

  test<1>(static_cast<Tria1::Settings>(Tria1::partition_zorder |
                                       Tria1::construct_multigrid_hierarchy));
  test<2>(static_cast<Tria2::Settings>(Tria2::partition_zorder |
                                       Tria2::construct_multigrid_hierarchy));
  test<3>(static_cast<Tria3::Settings>(Tria3::partition_zorder |
                                       Tria3::construct_multigrid_hierarchy));
}
