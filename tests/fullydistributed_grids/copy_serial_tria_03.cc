// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a serial triangulation with manifold and hanging nodes and copy it.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "../grid/tests.h"


template <int dim>
void
test(int n_refinements, const int n_subdivisions, MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria(
    Triangulation<dim>::limit_level_difference_at_vertices);


  const Point<dim> center(1, 0);
  const double     inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    basetria, center, inner_radius, outer_radius, n_subdivisions);
  // basetria.reset_all_manifolds ();
  for (int step = 0; step < n_refinements; ++step)
    {
      for (auto &cell : basetria.active_cell_iterators())
        {
          // if(cell->is_locally_owned ())
          for (const unsigned int v : GeometryInfo<2>::vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));
              if (std::fabs(distance_from_center - inner_radius) < 1e-10)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }
      basetria.execute_coarsening_and_refinement();
    }

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);
  // if(n_refinements!=0)
  GridTools::partition_multigrid_levels(basetria);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  tria_pft.set_manifold(0, SphericalManifold<dim>(center));

  // extract relevant information form pdt
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria,
      comm,
      TriangulationDescription::Settings::construct_multigrid_hierarchy);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);


  // test triangulation
  FE_Q<dim>       fe(2);
  DoFHandler<dim> dof_handler(tria_pft);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();

  // print statistics
  print_statistics(tria_pft, true);
  print_statistics(dof_handler, true);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const MPI_Comm comm = MPI_COMM_WORLD;

  {
    deallog.push("2d");
    const int n_refinements  = 3;
    const int n_subdivisions = 8;
    test<2>(n_refinements, n_subdivisions, comm);
    deallog.pop();
  }
}
