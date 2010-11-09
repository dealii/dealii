//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include <grid/tria.h>
#include <grid/grid_generator.h>

/**
 * A set of test meshes for the deal.II test suite.
 *
 * These meshes exhibit certain key features for writing tests. If you
 * want to test certain properties of algorithms, the following table
 * might be of help.
 *
 * <table border=1>
 * <tr><th>Mesh</th><th>Feature tested</th></tr>
 * <tr><td>#hypercube(tr)</td><td>works at all on a single
 * cell</td></tr>
 * <tr><td>#hypercube(tr,2)</td><td>works on uniform meshes</td></tr>
 * <tr><td>#hypercube(tr,3,true)</td><td>works with local
 * refinement</td></tr>
 * <tr><td>#star_shaped(tr,1)</td><td>method is robust if more than
 * usual cells meet in one vertex</td></tr>
 * <tr><td>#star_shaped(tr,2,true)</td><td>method is robust if more than
 * usual cells meet in one vertex and local refinement exceeds one
 * level</td></tr>
 * </table>
 *
 * @author Guido Kanschat
 * @date 2006, 2010
 */
namespace TestGrids
{
				   /**
				    * Generate grids based on
				    * hypercube. These meshes have a
				    * regular geometry and topology.
				    *
				    * @param <tt>refinement</tt>
				    * denotes the number of refinement
				    * steps of the root cell.
				    *
				    * @param if <tt>local</tt> is
				    * <tt>true</tt>, refine only the
				    * cell containing the corner with
				    * only negative coordinates.
				    */
  template <int dim>
  void hypercube(Triangulation<dim>& tr,
		 unsigned int refinement = 0,
		 bool local = false)
  {
    GridGenerator::hyper_cube(tr);
    if (refinement && !local)
      tr.refine_global(refinement);
    if (refinement && local)
      {
	tr.refine_global(1);
	for (unsigned int i=1;i<refinement;++i)
	  {
	    for (typename Triangulation<dim>::active_cell_iterator
		   cell = tr.begin_active(); cell != tr.end(); ++cell)
	      {
		const Point<dim>& p = cell->center();
		bool negative = true;
		for (unsigned int d=0;d<dim;++d)
		  if (p(d) >= 0.)negative = false;
	      }
	    tr.execute_coarsening_and_refinement();
	  }
      }
    deallog << "Triangulation hypercube " << dim << "D refinement " << refinement;
    if (local)
      deallog << " local ";
    deallog << " steps " << tr.n_active_cells() << " active cells "
	    << tr.n_cells() << " total cells " << std::endl;
  }
  
				   /**
				    * Create a star-shaped mesh,
				    * having more than the average
				    * <tt>2<sup>dim</sup></tt> cells
				    * in the central vertex.
				    *
				    * @param <tt>refinement</tt>
				    * denotes the number of refinement
				    * steps of the root mesh.
				    *
				    * @param if <tt>local</tt> is
				    * <tt>true</tt>, refine only one
				    * of the coarse cells.
				    */
  template <int dim>
  void star_shaped(Triangulation<dim>& tr,
		   unsigned int refinement = 0,
		   bool local = false);
				   /**
				    * Local refinement of every other
				    * cell in a checkerboard fashion.
				    */
  template <int dim>
  void checkers(Triangulation<dim>& tr);
				   /**
				    * Islands of local refinement
				    */
  template <int dim>
  void islands(Triangulation<dim>& tr);
				   /**
				    * Local refinement with an
				    * unrefined hole.
				    */
  template <int dim>
  void laguna(Triangulation<dim>& tr);
}
