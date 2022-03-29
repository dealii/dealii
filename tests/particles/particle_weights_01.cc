// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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



// verify particle weighting mechanism


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // initialize grid
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  {
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
  }

  // initialize particle handler
  const MappingQ1<dim>            mapping;
  Particles::ParticleHandler<dim> particle_handler(triangulation, mapping);
  {
    const auto local_bounding_box =
      GridTools::compute_mesh_predicate_bounding_box(
        triangulation, IteratorFilters::LocallyOwnedCell());
    const auto global_bounding_boxes =
      Utilities::MPI::all_gather(triangulation.get_communicator(),
                                 local_bounding_box);

    Particles::Generators::quadrature_points(triangulation,
                                             QMidpoint<dim>(),
                                             global_bounding_boxes,
                                             particle_handler);

    Assert(particle_handler.n_locally_owned_particles() ==
             triangulation.n_locally_owned_active_cells(),
           ExcInternalError());
    Assert(particle_handler.n_global_particles() ==
             triangulation.n_global_active_cells(),
           ExcInternalError());
  }

  // register weighting function that returns the number of particles per cell
  triangulation.signals.weight.connect(
    [&particle_handler](
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
        &cell,
      const typename parallel::distributed::Triangulation<dim>::CellStatus
        status) -> unsigned int {
      unsigned int n_particles_in_cell = 0;
      switch (status)
        {
          case parallel::distributed::Triangulation<dim>::CELL_PERSIST:
          case parallel::distributed::Triangulation<dim>::CELL_REFINE:
            n_particles_in_cell = particle_handler.n_particles_in_cell(cell);
            break;

          case parallel::distributed::Triangulation<dim>::CELL_INVALID:
            break;

          case parallel::distributed::Triangulation<dim>::CELL_COARSEN:
            for (const auto &child : cell->child_iterators())
              n_particles_in_cell +=
                particle_handler.n_particles_in_cell(child);
            break;

          default:
            Assert(false, ExcInternalError());
            break;
        }
      deallog << ' ' << n_particles_in_cell << std::endl;
      return n_particles_in_cell;
    });

  // refine one cell, and coarsen one group of siblings
  for (const auto &cell : triangulation.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    {
      if (cell->id().to_string() == "0_2:00")
        cell->set_refine_flag();
      else if (cell->parent()->id().to_string() == "0_1:3")
        cell->set_coarsen_flag();
    }

  deallog << "weights before repartitioning:" << std::endl;
  triangulation.execute_coarsening_and_refinement();
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  test<2>();
}
