// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Show the shape functions of the Nedelec element on a grid with only
// one cell. This cell is rotated, stretched, scaled, etc, and on each
// of these cells each time we evaluate the shape functions.

// Shape functions and their gradients are output in gnuplot format,
// prefixed by an identifier for the transformation. To produce a
// gnuplot file for the values of functions of degree K with
// transformation N, do

// perl -n -e 'print if s/DEAL:NedelecK-TransformN::value//' output

#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 2


Point<2>
stretch_coordinates(const Point<2> p)
{
  return Point<2>(2 * p[0], p[1]);
}



Point<2>
tilt_coordinates(const Point<2> p)
{
  return Point<2>(p[0] + p[1], p[1]);
}


void
transform_grid(Triangulation<2> &tria, const unsigned int transform)
{
  switch (transform)
    {
      // first round: take
      // original grid
      case 0:
        break;

      // second round: rotate
      // triangulation
      case 1:
        GridTools::rotate(numbers::PI_2, tria);
        break;

      // third round: inflate
      // by a factor of 2
      case 2:
        GridTools::scale(2, tria);
        break;

      // third round: scale
      // back, rotate back,
      // stretch
      case 3:
        GridTools::scale(.5, tria);
        GridTools::rotate(-numbers::PI_2, tria);
        GridTools::transform(&stretch_coordinates, tria);

        break;

      default:
        DEAL_II_NOT_IMPLEMENTED();
    };
}



template <int dim>
inline void
plot_shape_functions(const unsigned int degree)
{
  FE_Nedelec<dim>    element(degree);
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

  // check the following with a
  // number of transformed
  // triangulations
  for (unsigned int transform = 0; transform < 4; ++transform)
    {
      std::ostringstream ost;
      ost << "Nedelec" << degree << "-Transform" << transform;
      deallog.push(ost.str());

      transform_grid(tr, transform);

      DoFHandler<dim>                         dof(tr);
      typename DoFHandler<dim>::cell_iterator c = dof.begin();
      dof.distribute_dofs(element);

      QTrapezoid<1>      q_trapez;
      const unsigned int div = 2;
      QIterated<dim>     q(q_trapez, div);
      FEValues<dim>      fe(element,
                       q,
                       update_values | update_gradients |
                         update_quadrature_points);
      fe.reinit(c);

      for (unsigned int q_point = 0; q_point < q.size(); ++q_point)
        {
          // Output function in
          // gnuplot readable format,
          // namely x y z u0x u0y u0z u1x...
          deallog << "value    " << q_point << '\t'
                  << fe.quadrature_point(q_point);

          for (unsigned int i = 0; i < element.dofs_per_cell; ++i)
            {
              for (unsigned int c = 0; c < dim; ++c)
                deallog << '\t' << fe.shape_value_component(i, q_point, c);
            }

          // Output the gradients in
          // similar fashion
          deallog << std::endl
                  << "gradient " << q_point << '\t'
                  << fe.quadrature_point(q_point);

          for (unsigned int i = 0; i < element.dofs_per_cell; ++i)
            {
              for (unsigned int c = 0; c < dim; ++c)
                {
                  for (unsigned int d = 0; d < dim; ++d)
                    deallog << '\t'
                            << fe.shape_grad_component(i, q_point, c)[d];
                }
            }
          deallog << std::endl;

          if ((q_point + 1) % (2 * div - 1) == 0)
            {
              deallog << "value    " << std::endl;
              deallog << "gradient " << std::endl;
            }
        }

      deallog << std::endl;
      deallog.pop();
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  plot_shape_functions<2>(0);
  plot_shape_functions<2>(1);

  return 0;
}
