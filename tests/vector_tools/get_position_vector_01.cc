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

// Test VectorTools::get_position_vector() with Mapping as argument.

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
test()
{
  // surface mesh
  const unsigned int spacedim      = dim + 1;
  const unsigned int fe_degree     = 3;
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
  VectorTools::get_position_vector(MappingQGeneric<dim, spacedim>(4),
                                   dof_handler_dim,
                                   euler_vector_base);
  MappingFEField<dim, spacedim> mapping_base(dof_handler_dim,
                                             euler_vector_base);

  // clear manifold
  tria.reset_all_manifolds();

  // output mesh with with MappingQGeneric(degree=4)
  {
    DataOutBase::VtkFlags flags;

    DataOut<dim, DoFHandler<dim, spacedim>> data_out;
    data_out.set_flags(flags);
    data_out.attach_dof_handler(dof_handler);

    data_out.build_patches(
      MappingQGeneric<dim, spacedim>(4),
      fe_degree + 1,
      DataOut<dim,
              DoFHandler<dim, spacedim>>::CurvedCellRegion::curved_inner_cells);

#if false
    std::ofstream output("test.0.vtk");
    data_out.write_vtk(output);
#else
    data_out.write_vtk(deallog.get_file_stream());
#endif
  }

  // output mesh with with MappingFEField based on a vector constructed
  // with the triangulation without manifolds and the high-order base mapping
  {
    Vector<double> euler_vector(dof_handler_dim.n_dofs());
    VectorTools::get_position_vector(mapping_base,
                                     dof_handler_dim,
                                     euler_vector);
    MappingFEField<dim, spacedim> mapping(dof_handler_dim, euler_vector);
    DataOutBase::VtkFlags         flags;

    DataOut<dim, DoFHandler<dim, spacedim>> data_out;
    data_out.set_flags(flags);
    data_out.attach_dof_handler(dof_handler);

    data_out.build_patches(
      mapping,
      fe_degree + 1,
      DataOut<dim,
              DoFHandler<dim, spacedim>>::CurvedCellRegion::curved_inner_cells);

#if false
    std::ofstream output("test.1.vtk");
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

  test<1>();
}
