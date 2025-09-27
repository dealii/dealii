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


// Like mapping_fe_fields_01 but for deformed meshes tested for linear and
// quadratic mapping.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
class Solution : public Function<dim>
{
public:
  Solution()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &point, const unsigned int component) const
  {
    return std::sin(point[component] * 0.5 * numbers::PI);
  }
};

void
test(const unsigned int mapping_degree)
{
  const int dim = 2;

  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 4);

  FE_SimplexP<dim> fe(mapping_degree);
  FESystem<dim>    euler_fe(fe, dim);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  DoFHandler<dim> euler_dof_handler(tria);
  euler_dof_handler.distribute_dofs(euler_fe);

  Vector<double> euler_vector(euler_dof_handler.n_dofs());

  // TODO: not working (missing mapping)
  // VectorTools::get_position_vector(euler_dof_handler, euler_vector);

  MappingFE<dim> mapping_interpolation(FE_SimplexP<dim>(1));
  VectorTools::interpolate(mapping_interpolation,
                           euler_dof_handler,
                           Solution<dim>(),
                           euler_vector);

  MappingFEField<dim> mapping(euler_dof_handler, euler_vector);

  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    Vector<double> solution(dof_handler.n_dofs());
    data_out.add_data_vector(solution, "solution");

    data_out.build_patches(mapping, 2);

#if false
    std::ofstream output("test." + std::to_string(mapping_degree) +  ".vtk");
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

  test(1); // linear mapping
  test(2); // quadratic mapping
}
