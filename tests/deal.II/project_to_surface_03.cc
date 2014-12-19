// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// like the _02 test, but project to the faces of a distorted cell



#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> tria;
  StraightBoundary<dim> boundary;

  GridGenerator::hyper_cube(tria, 0, 1);

  const typename Triangulation<dim>::active_cell_iterator
  cell=tria.begin_active();

  // distort the cell a bit. all
  // faces but face 0 stay planar;
  // face 0 becomes a saddle
  cell->vertex(0)[0] -= 0.25;
  cell->vertex(6)[0] -= 0.25;
  cell->vertex(2)[0] += 0.25;
  cell->vertex(4)[0] += 0.25;


  for (unsigned int point=0; point<9; ++point)
    {
      // choose the 8 vertices of the
      // original unit cell as well
      // as the center point
      const Point<dim> trial_point = (point < 8
                                      ?
                                      GeometryInfo<dim>::unit_cell_vertex(point)
                                      :
                                      Point<dim>(.5,.5,.5));

      deallog << "Trial point = " << trial_point << std::endl;

      for (unsigned int e=0; e<GeometryInfo<dim>::quads_per_cell; ++e)
        {
          const typename Triangulation<dim>::quad_iterator
          quad = (dim > 2 ? cell->quad(e) :
                  *reinterpret_cast<const typename Triangulation<dim>::quad_iterator *>(&cell));

          deallog << "    Quad " << e << ", projected point=";

          const Point<dim> p = boundary.project_to_surface (quad, trial_point);
          deallog << p;
          deallog << "  (quad is from ";
          deallog << quad->vertex(0);
          deallog << " to ";
          deallog << quad->vertex(1);
          deallog << " to ";
          deallog << quad->vertex(2);
          deallog << " to ";
          deallog << quad->vertex(3);
          deallog << ")" << std::endl;

          // now make sure that p is
          // indeed closer to
          // trial_point than any of
          // the vertices of the quad
          for (unsigned int v=0; v<4; ++v)
            Assert (p.distance (trial_point) <=
                    quad->vertex(v).distance (trial_point),
                    ExcInternalError());
        }
      deallog << std::endl;
    }
}




int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  test<3>();
}
