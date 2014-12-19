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



// test StraightBoundary::project_to_surface for quads



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


class Rotate2d
{
public:
  Rotate2d (const double angle)
    :
    angle(angle)
  {}

  template <int spacedim>
  Point<spacedim> operator() (const Point<spacedim> p) const
  {
    Point<spacedim> q;
    q[0] = std::cos(angle)*p(0) - std::sin(angle) * p(1);
    q[1] = std::sin(angle)*p(0) + std::cos(angle) * p(1);
    for (unsigned d=2; d<spacedim; ++d)
      q[d] = p[d];
    return q;
  }
private:
  const double angle;
};


template <int dim>
void do_rotate (Triangulation<dim> &tria)
{
  GridTools::transform (Rotate2d(numbers::PI/4), tria);
}


void do_rotate (Triangulation<1> &)
{}



template <int dim>
void create_triangulation(const bool rotate,
                          Triangulation<dim> &tria)
{
  GridGenerator::hyper_cube(tria, 1., 3.);

  if (rotate)
    do_rotate (tria);
}


template <int dim>
void test ()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> tria;
  StraightBoundary<dim> boundary;

  for (unsigned int case_no=0; case_no<2; ++case_no)
    {
      deallog << "  Case " << case_no << std::endl;
      create_triangulation((case_no==1), tria);

      const typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active();
      Point<dim> trial_point;
      for (unsigned int i=0; i<dim; ++i)
        trial_point[i] = 1.5;

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
      tria.clear();
    }
}




int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  test<2>();
  test<3>();
}
