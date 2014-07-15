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
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>


using namespace dealii;



template <int dim>
void make_grid ()
{
  Point<dim> center;
  for (unsigned int i=0; i<dim; ++i)
    center[i] = .25;

  double radius=center.norm();

  HyperBallBoundary<dim,dim> boundary(center, .25*std::sqrt((double)dim));
  Triangulation<dim,dim> triangulation;
  GridGenerator::hyper_cube (triangulation);

  triangulation.refine_global(1);

  for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active();
      cell!=triangulation.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->center().distance(center)< radius)
        cell->face(f)->set_all_manifold_ids(1);

  triangulation.set_manifold(1,boundary);
  triangulation.refine_global(3);

  const std::string filename = "mesh-" + Utilities::int_to_string(dim) + "d.vtk";
  std::ofstream out (filename.c_str());
  GridOut grid_out;
  grid_out.write_vtk (triangulation, out);
}


  

// @sect3{The main function}

// Finally, the main function. There isn't much to do here, only to call the
// subfunctions.
int main ()
{
  make_grid<2> ();
  make_grid<3> ();
}
