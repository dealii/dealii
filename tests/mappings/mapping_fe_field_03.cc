// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2018 by the deal.II authors
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

// tests VectorTools::get_position_vector() used in MappingFEField
// with parallel distributed triangulation

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  {
    Point<dim> center;
    GridGenerator::hyper_ball(tria, center, 12);

    const types::manifold_id            sphere_id = 0;
    static const SphericalManifold<dim> boundary_ball(center);
    tria.set_all_manifold_ids_on_boundary(sphere_id);
    tria.set_manifold(sphere_id, boundary_ball);

    tria.refine_global(dim == 2 ? 1 : 3);
  }

  FESystem<dim>   fe(FE_Q<dim>(4), dim);
  DoFHandler<dim> dh(tria);

  dh.distribute_dofs(fe);

  deallog << "dim: " << dim << std::endl
          << "cells: " << tria.n_global_active_cells()
          << ", dofs: " << dh.n_dofs() << std::endl;

  // Create a Mapping
  LinearAlgebra::distributed::Vector<double> map_vector;

  IndexSet locally_relevant_euler;
  DoFTools::extract_locally_relevant_dofs(dh, locally_relevant_euler);
  map_vector.reinit(dh.locally_owned_dofs(),
                    locally_relevant_euler,
                    MPI_COMM_WORLD);

  VectorTools::get_position_vector(dh, map_vector);

  deallog << "L2=" << map_vector.l2_norm() << std::endl
          << "L1=" << map_vector.l1_norm() << std::endl
          << "Linfty=" << map_vector.linfty_norm() << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    initlog();

  test<2>();
  test<3>();
}
