// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// check MappingFEField when initialized manually on levels on a distributed
// triangulation

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  parallel::distributed::Triangulation<dim, spacedim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim, spacedim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim, spacedim>::
      construct_multigrid_hierarchy);
  GridGenerator::hyper_ball(tria);

  tria.refine_global(2);
  if (tria.begin_active()->is_locally_owned())
    tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  FESystem<dim, spacedim>   fe(FE_Q<dim, spacedim>(2), spacedim);
  DoFHandler<dim, spacedim> dh(tria);

  dh.distribute_dofs(fe);
  dh.distribute_mg_dofs();

  deallog << "dim, spacedim: " << dim << ", " << spacedim << std::endl
          << "cells: " << tria.n_active_cells() << ", dofs: " << dh.n_dofs()
          << std::endl;

  // Create a Mapping
  std::vector<LinearAlgebra::distributed::Vector<double>> level_vectors(
    tria.n_global_levels());
  MappingQGeneric<dim, spacedim> mapping_ref(fe.degree);
  FEValues<dim>                  fe_values_setup(mapping_ref,
                                dh.get_fe(),
                                Quadrature<dim>(
                                  dh.get_fe().get_unit_support_points()),
                                update_quadrature_points);
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    {
      IndexSet relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(dh, level, relevant_dofs);
      level_vectors[level].reinit(dh.locally_owned_mg_dofs(level),
                                  relevant_dofs,
                                  tria.get_communicator());
      std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
      for (const auto &cell : dh.mg_cell_iterators_on_level(level))
        if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
          {
            fe_values_setup.reinit(cell);
            cell->get_active_or_mg_dof_indices(dof_indices);
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
              {
                const unsigned int coordinate_direction =
                  fe.system_to_component_index(i).first;
                const Point<dim> point = fe_values_setup.quadrature_point(i);
                level_vectors[level](dof_indices[i]) =
                  point[coordinate_direction];
              }
          }
      level_vectors[level].update_ghost_values();
    }

  MappingFEField<dim,
                 spacedim,
                 LinearAlgebra::distributed::Vector<double>,
                 DoFHandler<dim>>
    mapping(dh, level_vectors);

  QGauss<dim>   quad(1);
  FEValues<dim> fe_values_ref(mapping_ref,
                              fe,
                              quad,
                              update_jacobians | update_quadrature_points);
  FEValues<dim> fe_values(mapping,
                          fe,
                          quad,
                          update_jacobians | update_quadrature_points);

  for (const auto &cell : tria.cell_iterators())
    if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
      {
        fe_values_ref.reinit(cell);
        fe_values.reinit(cell);

        if (fe_values_ref.quadrature_point(0).distance(
              fe_values.quadrature_point(0)) > 1e-12)
          deallog << "Mapped point should be "
                  << fe_values_ref.quadrature_point(0) << " and is "
                  << fe_values.quadrature_point(0) << std::endl;
        Tensor<2, dim> jac_ref = fe_values_ref.jacobian(0),
                       jac     = fe_values.jacobian(0);
        if ((jac_ref - jac).norm() > 1e-12)
          deallog << "Jacobian should be " << jac_ref << " and is " << jac
                  << std::endl;
      }

  // shift the Euler vectors and check whether the result is still correct
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    level_vectors[level].add(1.1);

  for (const auto &cell : tria.cell_iterators())
    if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
      {
        fe_values_ref.reinit(cell);
        fe_values.reinit(cell);

        Point<dim> shift;
        for (unsigned int d = 0; d < dim; ++d)
          shift[d] = 1.1;
        if ((fe_values_ref.quadrature_point(0) + shift)
              .distance(fe_values.quadrature_point(0)) > 1e-12)
          deallog << "Mapped point should be "
                  << fe_values_ref.quadrature_point(0) + shift << " and is "
                  << fe_values.quadrature_point(0) << std::endl;
        Tensor<2, dim> jac_ref = fe_values_ref.jacobian(0),
                       jac     = fe_values.jacobian(0);
        if ((jac_ref - jac).norm() > 1e-12)
          deallog << "Jacobian should be " << jac_ref << " and is " << jac
                  << std::endl;
      }
  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc,
                                            argv,
                                            testing_max_num_threads());
  MPILogInitAll                    log;
  test<2, 2>();
  test<3, 3>();
}
