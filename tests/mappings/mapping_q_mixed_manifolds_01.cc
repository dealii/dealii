// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Compute the area of a square where one boundary is deformed by a
// cylindrical manifold.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <sstream>
#include <string>
#include <vector>

#include "../tests.h"

template <int dim>
struct ManifoldWrapper
{
  Manifold<dim> *
  operator()(const Tensor<1, dim> &direction, const Point<dim> &center) const;
};

template <>
Manifold<2> *
ManifoldWrapper<2>::operator()(const Tensor<1, 2> & /*direction*/,
                               const Point<2> &center) const
{
  return new SphericalManifold<2>(center);
}

template <>
Manifold<3> *
ManifoldWrapper<3>::operator()(const Tensor<1, 3> &direction,
                               const Point<3>     &center) const
{
  return new CylindricalManifold<3>(direction, center);
}

template <int dim>
void
test()
{
  Tensor<1, dim> direction;
  direction[dim - 1] = 1.;

  std::shared_ptr<Manifold<dim>> cylinder_manifold(
    ManifoldWrapper<dim>()(direction, Point<dim>()));

  // create cube and shift it to position 0.7
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -0.5, 0.5);
  Tensor<1, dim> shift;
  shift[0] = 1.;
  GridTools::shift(shift, tria);
  tria.begin()->face(0)->set_all_manifold_ids(1);
  tria.set_manifold(1, *cylinder_manifold);

  FE_Nothing<dim> fe;
  for (unsigned int degree = 6; degree < 7; ++degree)
    {
      MappingQ<dim> mapping(degree);

      QGauss<dim>   quad(degree + 1);
      FEValues<dim> fe_values(mapping, fe, quad, update_JxW_values);
      double        sum = 0.;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tria.begin_active();
           cell != tria.end();
           ++cell)
        {
          fe_values.reinit(cell);
          double local_sum = 0;
          for (unsigned int q = 0; q < quad.size(); ++q)
            local_sum += fe_values.JxW(q);
          sum += local_sum;
        }
      // reference area/volume is the area of the cube minus the area of the
      // circular segment that is missing due to the manifold, whose area is
      // R^2/2(theta-sin(theta)) with theta 90 degrees and R=sqrt(0.5)
      const double reference = 1. - 0.25 * (0.5 * numbers::PI - 1.);
      deallog << "Volume " << dim << "D mapping degree " << degree << ": "
              << sum << " error: " << (sum - reference) / reference
              << std::endl;
    }
}


int
main()
{
  initlog();

  test<2>();
  test<3>();
}
