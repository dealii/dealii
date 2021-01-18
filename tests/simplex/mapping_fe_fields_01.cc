// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Test MappingFEField for simplex meshes.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/simplex/fe_lib.h>
#include <deal.II/simplex/grid_generator.h>
#include <deal.II/simplex/quadrature_lib.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
class Solution : public Function<dim>
{
public:
  Solution()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &point, const unsigned int compontent) const
  {
    return point[compontent];
  }
};

void
test()
{
  const int dim = 2;

  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 1);

  Simplex::FE_P<dim> fe(2);
  FESystem<dim>      euler_fe(fe, dim);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  DoFHandler<dim> euler_dof_handler(tria);
  euler_dof_handler.distribute_dofs(euler_fe);

  Vector<double> euler_vector(euler_dof_handler.n_dofs());

  // TODO: not working (missing mapping)
  // VectorTools::get_position_vector(euler_dof_handler, euler_vector);

  MappingFE<dim> mapping_interpolation(Simplex::FE_P<dim>(1));
  VectorTools::interpolate(mapping_interpolation,
                           euler_dof_handler,
                           Solution<dim>(),
                           euler_vector);

  MappingFEField<dim> mapping(euler_dof_handler, euler_vector);

  Simplex::QGauss<dim> quadrature_formula(1);

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
    }

  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    Vector<double> solution(dof_handler.n_dofs());
    data_out.add_data_vector(solution, "solution");

    data_out.build_patches(mapping, 2);

#if false
    std::ofstream output("test.vtk");
    data_out.write_vtk(output);
#else
    data_out.write_vtk(deallog.get_file_stream());
#endif
  }
}

int
main()
{
  initlog();

  test();
}