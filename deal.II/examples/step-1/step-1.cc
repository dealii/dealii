/* $Id$ */

				 // Needed for the Triangulation class:
#include <grid/tria.h>
				 // We need the following two includes
				 // for loops over cells and/or faces:
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
				 // Here are some functions to
				 // generate standard grids:
#include <grid/grid_generator.h>
				 // We would like to use boundaries
				 // which are not straight lines, so
				 // we import some classes which
				 // predefine some boundary
				 // descriptions:
#include <grid/tria_boundary_lib.h>
				 // Output of grids in various
				 // graphics formats:
#include <basic/grid_out.h>

				 // This is needed for C++ output:
#include <fstream>



				 // In the following function, we
				 // simply use the unit square as
				 // domain and produce a globally
				 // refined grid from it.
void first_grid ()
{
				   // Define an object for a
				   // triangulation of a
				   // two-dimensional domain
  Triangulation<2> triangulation;
  
				   // Fill it with a square
  GridGenerator::hyper_cube (triangulation);
  
				   // Refine all cells four times, to
				   // yield 4^4=256 cells in total
  triangulation.refine_global (4);

				   // Now we want to write it to some
				   // output, here in postscript
				   // format
  ofstream out ("grid-1.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);
};



				 // The grid in the following function
				 // is slightly more complicated in
				 // that we use a ring domain and
				 // refine the result once globally
void second_grid ()
{
				   // Define an object for a
				   // triangulation of a
				   // two-dimensional domain
  Triangulation<2> triangulation;
  
				   // Fill it with a ring domain. The
				   // center of the ring shall be the
				   // point (1,0), and inner and outer
				   // radius shall be 0.5 and 1. The
				   // number of circumferentical cells
				   // will be adjusted automatically
				   // by this function (in this case,
				   // there will be 10)
  const Point<2> center (1,0);
  const double inner_radius = 0.5,
	       outer_radius = 1.0;
  GridGenerator::hyper_shell (triangulation,
			      center, inner_radius, outer_radius);
				   // By default, the triangulation
				   // assumes that all boundaries are
				   // straight and given by the cells
				   // of the coarse grid (which we
				   // just created). Here, however, we
				   // would like to have a curved
				   // boundary. Furtunately, some good
				   // soul implemented an object which
				   // describes the boundary of a ring
				   // domain; it only needs the center
				   // of the ring and automatically
				   // figures out the inner and outer
				   // radius when needed. Note that we
				   // associate this boundary object
				   // with that part of the boundary
				   // that has the "boundary number"
				   // zero. By default, all boundary
				   // parts have this number, but you
				   // might want to change this number
				   // for some parts, and then the
				   // curved boundary thus associated
				   // with number zero will not apply
				   // there.
  const HyperShellBoundary<2> boundary_description(center);
  triangulation.set_boundary (0, boundary_description);
  
				   // Now, just for the purpose of
				   // demonstration and for no
				   // particular reason, we will
				   // refine the grid in five steps
				   // towards the inner circle of the
				   // domain:
  for (unsigned int step=0; step<5; ++step)
    {
				       // Get an iterator which points
				       // to a cell and which we will
				       // move over all active cells
				       // one by one. Active cells are
				       // those that are not further
				       // refined
      Triangulation<2>::active_cell_iterator cell, endc;
      cell = triangulation.begin_active();
      endc = triangulation.end();

				       // Now loop over all cells...
      for (; cell!=endc; ++cell)
					 // ...and over all vertices
					 // of the cells. Note the
					 // dimension-independent way
					 // by which we find out about
					 // the number of faces of a
					 // cell
	for (unsigned int vertex=0;
	     vertex < GeometryInfo<2>::vertices_per_cell;
	     ++vertex)
	  {
					     // If this cell is at the
					     // inner boundary, then
					     // at least one of its vertices
					     // must have a radial
					     // distance from the center
					     // of 0.5
	    const Point<2> vector_to_center
	      = (cell->vertex(vertex) - center);
	    const double distance_from_center
	      = sqrt(vector_to_center.square());
	    
	    if (fabs(distance_from_center - inner_radius) < 1e-10)
	      {
						 // Ok, this is one of
						 // the cells we were
						 // looking for. Flag
						 // it for refinement
						 // and go to the next
						 // cell by breaking
						 // the loop over all
						 // vertices
		cell->set_refine_flag ();
		break;
	      };
	  };

				       // Refine the cells which we
				       // have marked
      triangulation.execute_coarsening_and_refinement ();
    };
  
  
				   // Now we want to write it to some
				   // output, here in postscript
				   // format
  ofstream out ("grid-2.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);


				   // At this point, all objects
				   // created in this function will be
				   // destroyed in reverse
				   // order. Unfortunately, we defined
				   // the boundary object after the
				   // triangulation, which still has a
				   // pointer to it and the library
				   // will produce an error if the
				   // boundary object is destroyed
				   // before the triangulation. We
				   // therefore have to release it,
				   // which can be done as
				   // follows. Note that this sets the
				   // boundary object used for part
				   // "0" of the boundary back to a
				   // default object, over which the
				   // triangulation has full control.
  triangulation.set_boundary (0);
};



int main () 
{
  first_grid ();
  second_grid ();
};
