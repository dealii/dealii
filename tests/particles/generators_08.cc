// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Insert particles within an hyper_cube triangulation using a Q1 quadrature
// defined on a non-matching hyper_cube triangulation and then check if the
// particles are correctly positioned


#include <deal.II/base/function_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  parallel::distributed::Triangulation<dim, spacedim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  DoFHandler<dim, spacedim>       dof_handler(tr);
  const FE_Nothing<dim, spacedim> fe_nothing;
  dof_handler.distribute_dofs(fe_nothing);

  const MappingQ<dim, spacedim> mapping(1);

  Particles::ParticleHandler<dim, spacedim> particle_handler(tr, mapping);

  parallel::distributed::Triangulation<dim, spacedim> particles_tr(
    MPI_COMM_WORLD);
  GridGenerator::hyper_cube(particles_tr, 0.1, 0.9);

  const QGauss<dim> quadrature(2);

  // Generate the necessary bounding boxes for the generator
  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    tr, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

  Particles::Generators::quadrature_points(particles_tr,
                                           quadrature,
                                           global_bounding_boxes,
                                           particle_handler);

  {
    deallog << "Locally owned active cells: "
            << tr.n_locally_owned_active_cells() << std::endl;

    deallog << "Global particles: " << particle_handler.n_global_particles()
            << std::endl;

    for (const auto &cell : tr.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            deallog << "Cell " << cell << " has "
                    << particle_handler.n_particles_in_cell(cell)
                    << " particles." << std::endl;
          }
      }

    for (const auto &particle : particle_handler)
      {
        deallog << "Particle index " << particle.get_id() << " is in cell "
                << particle.get_surrounding_cell() << std::endl;
        deallog << "Particle location: " << particle.get_location()
                << std::endl;
      }
  }
  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll init;

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("2d/3d");
  test<2, 3>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
