// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Test
// TriangulationDescription::Utilities::create_description_from_triangulation()
// so that it also works for p:d:T set up on a subcommunicator.

#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out.h>

#include "../grid/tests.h"

using namespace dealii;

namespace dealii
{
  namespace parallel
  {
    template <int dim, int spacedim = dim>
    std::function<
      unsigned int(const typename Triangulation<dim, spacedim>::cell_iterator &,
                   const CellStatus)>
    hanging_nodes_weighting(const double weight)
    {
      return [weight](const auto &cell, const auto &) -> unsigned int {
        if (cell->is_locally_owned() == false)
          return 0;

        bool flag = false;

        for (const auto f : cell->face_indices())
          if (!cell->at_boundary(f) &&
              (cell->neighbor(f)->has_children() ||
               cell->level() != cell->neighbor(f)->level()))
            flag = true;

        if (flag)
          return 10000 * weight;
        else
          return 10000;
      };
    }
  } // namespace parallel
} // namespace dealii


template <int dim, int spacedim>
void
create_triangulation(Triangulation<dim, spacedim> &tria)
{
  const unsigned int n_ref_global = 3;
  const unsigned int n_ref_local  = 3;
  GridGenerator::hyper_cube(tria, -1.0, +1.0);
  tria.refine_global(n_ref_global);

  for (unsigned int i = 0; i < n_ref_local; ++i)
    {
      for (auto cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool flag = true;
            for (unsigned int d = 0; d < dim; ++d)
              if (cell->center()[d] > 0.0)
                flag = false;
            if (flag)
              cell->set_refine_flag();
          }
      tria.execute_coarsening_and_refinement();
    }
}


template <int dim>
void
test(const MPI_Comm comm)
{
  parallel::distributed::Triangulation<dim> tria(
    comm,
    Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  create_triangulation(tria);

  const RepartitioningPolicyTools::FirstChildPolicy<dim> policy(tria);

  // const RepartitioningPolicyTools::CellWeightPolicy<dim>
  // policy(parallel::hanging_nodes_weighting<dim>(2.0));

  const auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      tria,
      policy.partition(tria),
      TriangulationDescription::Settings::construct_multigrid_hierarchy);

  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);
  tria_pft.create_triangulation(construction_data);

  if (false)
    {
      GridOut grid_out;
      grid_out.write_mesh_per_processor_as_vtu(tria_pft, "test", true, true);
    }

  FE_Q<dim> fe(2);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  DoFHandler<dim> dof_handler_pft(tria_pft);
  dof_handler_pft.distribute_dofs(fe);
  dof_handler_pft.distribute_mg_dofs();

  // print statistics
  print_statistics(tria, true);
  print_statistics(tria_pft, true);
  print_statistics(dof_handler, true);
  print_statistics(dof_handler_pft, true);
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  MPI_Comm comm = MPI_COMM_WORLD;

  deallog.push("all");
  test<2>(comm);
  deallog.pop();
}
