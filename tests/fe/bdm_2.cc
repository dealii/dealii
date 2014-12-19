// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Show the shape functions of the Raviart-Thomas element on a grid
// with only one cell. This cell is rotated, stretched, scaled, etc,
// and on each of these cells each time we evaluate the shape
// functions.

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#define PRECISION 3


Point<2> stretch_coordinates (const Point<2> p)
{
  return Point<2>(2*p(0), p(1));
}



Point<2> tilt_coordinates (const Point<2> p)
{
  return Point<2>(p(0)+p(1), p(1));
}


void
transform_grid (Triangulation<2> &tria,
                const unsigned int  transform)
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
      GridTools::rotate (3.14159265358/2, tria);
      break;

      // third round: inflate
      // by a factor of 2
    case 2:
      GridTools::scale (2, tria);
      break;

      // third round: scale
      // back, rotate back,
      // stretch
    case 3:
      GridTools::scale (.5, tria);
      GridTools::rotate (-3.14159265358/2, tria);
      GridTools::transform (&stretch_coordinates, tria);

      break;

    default:
      Assert (false, ExcNotImplemented());
    };
}




template<int dim>
void
plot_shape_functions(const unsigned int degree)
{
  FE_BDM<dim> element(degree);
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, 0., 1.);

  // check the following with a
  // number of transformed
  // triangulations
  for (unsigned int transform=0; transform<4; ++transform)
    {
      std::ostringstream ost;
      ost << "BDM" << degree << "-Transform" << transform;
      deallog.push(ost.str());

      transform_grid (tr, transform);

      DoFHandler<dim> dof(tr);
      typename DoFHandler<dim>::cell_iterator c = dof.begin();
      dof.distribute_dofs(element);

      QTrapez<1> q_trapez;
      const unsigned int div=2;
      QIterated<dim> q(q_trapez, div);
      FEValues<dim> fe(element, q, update_values|update_gradients|update_q_points);
      fe.reinit(c);

      for (unsigned int q_point=0; q_point< q.size(); ++q_point)
        {
          // Output function in
          // gnuplot readable format,
          // namely x y z u0x u0y u0z u1x...
          deallog << "value    " << q_point << '\t' << fe.quadrature_point(q_point);

          for (unsigned int i=0; i<element.dofs_per_cell; ++i)
            {
              for (unsigned int c=0; c<dim; ++c)
                deallog << '\t' << fe.shape_value_component(i,q_point,c);
            }

          // Output the gradients in
          // similar fashion
          deallog << std::endl << "gradient " << q_point << '\t' << fe.quadrature_point(q_point);

          for (unsigned int i=0; i<element.dofs_per_cell; ++i)
            {
              for (unsigned int c=0; c<dim; ++c)
                {
                  for (unsigned int d=0; d<dim; ++d)
                    deallog << '\t' << fe.shape_grad_component(i,q_point,c)[d];
                }
            }
          deallog << std::endl;

          if ((q_point+1) % (2*div-1) == 0)
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
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int degree=1; degree<4; ++degree)
    plot_shape_functions<2>(degree);

  return 0;
}
