/* ---------------------------------------------------------------------
 * $Id$
 *
 * Copyright (C) 2014 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2014
 */

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>


using namespace dealii;


class SineBoundary : public Boundary<3>
{
public:
  virtual
  Point<3>
  get_new_point_on_line (const typename Triangulation<3>::line_iterator &line) const;

  virtual
  Point<3>
  get_new_point_on_quad (const typename Triangulation<3>::quad_iterator &quad) const;
};
  

Point<3>
SineBoundary::
get_new_point_on_line (const typename Triangulation<3>::line_iterator &line) const
{
  const double new_x = (line->vertex(0)[0] + line->vertex(1)[0]) / 2;
  const double new_y = (line->vertex(0)[1] + line->vertex(1)[1]) / 2;
  const double new_z = 1 +
		       0.05 * std::sin (2 * numbers::PI * new_x) +
		       0.05 * std::sin (2 * numbers::PI * new_y);
  return Point<3> (new_x, new_y, new_z);
}


Point<3>
SineBoundary::
get_new_point_on_quad (const typename Triangulation<3>::quad_iterator &quad) const
{
  const double new_x = (quad->vertex(0)[0] + quad->vertex(1)[0] +
			quad->vertex(2)[0] + quad->vertex(3)[0]) / 4;
  const double new_y = (quad->vertex(0)[1] + quad->vertex(1)[1] +
			quad->vertex(2)[1] + quad->vertex(3)[1]) / 4;
  const double new_z = 1 +
		       0.05 * std::sin (2 * numbers::PI * new_x) +
		       0.05 * std::sin (2 * numbers::PI * new_y);
  return Point<3> (new_x, new_y, new_z);
}


void make_grid ()
{
  SineBoundary sine_boundary;
  
  Triangulation<3> triangulation;
  GridGenerator::hyper_cube (triangulation);

  triangulation.begin_active()->face(5)->set_manifold_id(1);
  triangulation.set_boundary (1, sine_boundary);

  triangulation.refine_global (2);

  std::ofstream out ("mesh.mgl");
  GridOut grid_out;
  grid_out.write_mathgl (triangulation, out);
}


  

// @sect3{The main function}

// Finally, the main function. There isn't much to do here, only to call the
// subfunctions.
int main ()
{
  make_grid ();
}
