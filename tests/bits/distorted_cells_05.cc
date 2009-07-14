//----------------------------  distorted_cells_05.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  distorted_cells_05.cc  ---------------------------


// like _04, except that this time the cell really can't be fixed up
// because we move the node from the bottom face all the way above the
// top face

#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_tools.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>

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
	  return Point<dim>(0,1.25);
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

	return Point<dim>(0,0,1.25);
      }
};



template <int dim>
void check ()
{
  MyBoundary<dim> my_boundary;
  
				   // create a single square/cube
  Triangulation<dim> coarse_grid;
  GridGenerator::hyper_cube (coarse_grid, -1, 1);

				   // set bottom face to use MyBoundary
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    if (coarse_grid.begin_active()->face(f)->center()[dim-1] == -1)
      coarse_grid.begin_active()->face(f)->set_boundary_indicator (1);
  coarse_grid.set_boundary (1, my_boundary);

				   // now try to refine this one
				   // cell. we should get an exception
  try
    {
      coarse_grid.begin_active()->set_refine_flag ();
      coarse_grid.execute_coarsening_and_refinement ();
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      typename Triangulation<dim>::DistortedCellList
	subset = GridTools::fix_up_distorted_child_cells (dcv,
							  coarse_grid);
      Assert (subset.distorted_cells.size() == 1,
	      ExcInternalError());
    }

  Assert (coarse_grid.n_levels() == 2, ExcInternalError());
  Assert (coarse_grid.n_active_cells() == 1<<dim, ExcInternalError());

				   // output the coordinates of the
				   // child cells
  GridOut().write_gnuplot (coarse_grid, deallog.get_file_stream());
}


int main () 
{
  std::ofstream logfile("distorted_cells_05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2> ();
  check<3> ();
}

  
  
