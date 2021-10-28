// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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



// Test DataOutFaces::write_vtk() for simplex meshes.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction(const unsigned int n_components)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    return p[component % dim] * p[component % dim];
  }
};

template <int dim, int spacedim = dim>
void
test(const FiniteElement<dim, spacedim> &fe, const unsigned int n_components)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, dim == 2 ? 4 : 2);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());

  VectorTools::interpolate(dof_handler,
                           RightHandSideFunction<dim>(n_components),
                           solution);

  DataOutFaces<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

#if true
  static unsigned int counter = 0;
  std::ofstream       output("test." + std::to_string(dim) +
                       std::to_string(counter++) + ".vtk");
  data_out.write_vtk(output);
#endif

  data_out.write_vtk(deallog.get_file_stream());
}

int
main()
{
  initlog();

  {
    const int dim = 2;
    test<dim>(FE_SimplexP<dim>(2) /*=degree*/, 1);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/), dim), dim);
  }

  {
    const int dim = 3;
    test<dim>(FE_SimplexP<dim>(2) /*=degree*/, 1);
    test<dim>(FESystem<dim>(FE_SimplexP<dim>(2 /*=degree*/), dim), dim);
  }
}
