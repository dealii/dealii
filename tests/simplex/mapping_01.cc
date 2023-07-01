// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
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

// Test mapping centers with simplices

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

using namespace dealii;

template <int dim>
class Position : public Function<dim>
{
public:
  Position()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &point, const unsigned int component) const
  {
    if (component == 0)
      return point[component] + std::sin(point[dim - 1] * numbers::PI);
    return point[component];
  }
};

void
test(const unsigned int mapping_degree)
{
  const int dim = 3;

  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 2);

  FE_SimplexP<dim> fe(mapping_degree);
  FESystem<dim>    position_fe(fe, dim);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  DoFHandler<dim> position_dof_handler(tria);
  position_dof_handler.distribute_dofs(position_fe);

  Vector<double> position_vector(position_dof_handler.n_dofs());
  MappingFE<dim> mapping_interpolation(FE_SimplexP<dim>(1));
  VectorTools::interpolate(mapping_interpolation,
                           position_dof_handler,
                           Position<dim>(),
                           position_vector);

  MappingFEField<dim> mapping(position_dof_handler, position_vector);

  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << "cell = " << cell << std::endl;
      deallog << "mapped barycenter       = " << mapping.get_center(cell, true)
              << std::endl;
      Point<dim> manually_mapped_center;
      for (unsigned int d = 0; d < dim; ++d)
        manually_mapped_center[d] =
          Position<dim>().value(cell->barycenter(), d);
      deallog << "exact mapped barycenter = " << manually_mapped_center
              << std::endl;
      deallog << "mapped vertex mean      = " << mapping.get_center(cell, false)
              << std::endl;
    }
}

int
main()
{
  initlog();

  test(1); // linear mapping

  deallog << std::endl << "quadratic mapping" << std::endl << std::endl;

  test(2); // quadratic mapping
}
