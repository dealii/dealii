// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test VectorTools::get_position_vector() with Mapping as argument.

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
test(const bool         vector_describes_relative_displacement,
     const unsigned int fe_degree,
     const unsigned int fe_degree_mapping)
{
  // surface mesh
  const unsigned int spacedim      = dim + 1;
  const unsigned int n_refinements = 1;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_sphere(tria, Point<spacedim>(), 0.5);
  tria.refine_global(n_refinements);

  FE_Q<dim, spacedim>       fe(fe_degree);
  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  FESystem<dim, spacedim>   fe_dim(fe, spacedim);
  DoFHandler<dim, spacedim> dof_handler_dim(tria);
  dof_handler_dim.distribute_dofs(fe_dim);

  // set up base high-order mapping
  Vector<double> euler_vector_base(dof_handler_dim.n_dofs());
  VectorTools::get_position_vector(dof_handler_dim, euler_vector_base);
  MappingFEField<dim, spacedim> mapping_base(dof_handler_dim,
                                             euler_vector_base);

  // clear manifold
  tria.reset_all_manifolds();

  auto euler_vector = euler_vector_base;

  if (vector_describes_relative_displacement)
    {
      Vector<double> absolute_vector(dof_handler_dim.n_dofs());
      VectorTools::get_position_vector(dof_handler_dim, absolute_vector);

      euler_vector -= absolute_vector;
    }

  // output mesh with with MappingFEField based on a vector constructed
  // with the triangulation without manifolds and the high-order base mapping
  {
    MappingQCache<dim, spacedim> mapping(fe_degree_mapping);
    mapping.initialize(MappingQ1<dim, spacedim>(),
                       dof_handler_dim,
                       euler_vector,
                       vector_describes_relative_displacement);
    DataOutBase::VtkFlags flags;

    DataOut<dim, spacedim> data_out;
    data_out.set_flags(flags);
    data_out.attach_dof_handler(dof_handler);

    data_out.build_patches(
      mapping,
      fe_degree + 1,
      DataOut<dim, spacedim>::CurvedCellRegion::curved_inner_cells);

    static unsigned int counter = 0;

#if false
    std::ofstream output("test." + std::to_string(counter++) + ".vtk");
    data_out.write_vtk(output);
#else
    data_out.write_vtk(deallog.get_file_stream());
#endif
  }
}



int
main(int argc, char **argv)
{
  initlog();

  test<1>(false, 4, 4);
  test<1>(true, 4, 4);

  test<1>(false, 4, 5);
  test<1>(true, 4, 5);
}
