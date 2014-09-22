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



// check that indeed Triangulation::execute_coarsening_and_refinement
// throws an exception if a distorted child cell is created
//
// the 2d case is as described in the example in the Glossary on
// distorted cells resulting from refinement

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>


template <int dim>
class MyBoundary : public Boundary<dim>
{
  virtual Point<dim>
  get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
  {
    deallog << "Finding point between "
            << line->vertex(0) << " and "
            << line->vertex(1) << std::endl;

    // in 2d, find a point that
    // lies on the opposite side
    // of the quad. in 3d, choose
    // the midpoint of the edge
    if (dim == 2)
      return Point<dim>(0,0.75);
    else
      return (line->vertex(0) + line->vertex(1)) / 2;
  }

  virtual Point<dim>
  get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
  {
    deallog << "Finding point between "
            << quad->vertex(0) << " and "
            << quad->vertex(1) << " and "
            << quad->vertex(2) << " and "
            << quad->vertex(3) << std::endl;

    return Point<dim>(0,0,.75);
  }
};



template <int dim>
void check ()
{
  MyBoundary<dim> my_boundary;

  // create a single square/cube
  Triangulation<dim> coarse_grid (Triangulation<dim>::none, true);
  GridGenerator::hyper_cube (coarse_grid, -1, 1);

  // set bottom face to use MyBoundary
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    if (coarse_grid.begin_active()->face(f)->center()[dim-1] == -1)
      coarse_grid.begin_active()->face(f)->set_boundary_indicator (1);
  coarse_grid.set_boundary (1, my_boundary);

  // now try to refine this one
  // cell. we should get an exception
  bool flag = false;
  try
    {
      coarse_grid.begin_active()->set_refine_flag ();
      coarse_grid.execute_coarsening_and_refinement ();
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      flag = true;

      deallog << dcv.distorted_cells.size() << " distorted cells"
              << std::endl;
      Assert (dcv.distorted_cells.front() == coarse_grid.begin(0),
              ExcInternalError());

      // ignore the exception, i.e. do not
      // even attempt to fix up the problem
      // (but see the _04 testcase for this)
    }

  Assert (flag == true, ExcInternalError());
  Assert (coarse_grid.n_levels() == 2, ExcInternalError());
  Assert (coarse_grid.n_active_cells() == 1<<dim, ExcInternalError());

  // output the coordinates of the
  // child cells
  GridOut().write_gnuplot (coarse_grid, deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2> ();
  check<3> ();
}



