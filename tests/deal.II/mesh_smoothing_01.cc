//----------------------------  mesh_smoothing_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_smoothing_01.cc  ---------------------------


// something went wrong with smoothing a mesh when a particular set of
// cells were set to be refined. found in step-31
//
// a redux is in mesh_smoothing_02


char logname[] = "mesh_smoothing_01/output";


#include "../tests.h"


#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>

#include <fstream>
#include <iostream>
#include <cstring>


				 // Finally, this is as in previous
				 // programs:
using namespace dealii;



template <int dim>
bool cell_is_patch_level_1 (const typename Triangulation<dim>::cell_iterator &cell)
{
  Assert (cell->active() == false, ExcInternalError());

  unsigned int n_active_children = 0;
  for (unsigned int i=0; i<cell->n_children(); ++i)
    if (cell->child(i)->active())
      ++n_active_children;

  return (n_active_children == 0) || (n_active_children == cell->n_children());
}


void test ()
{
  Triangulation<2> triangulation (Triangulation<2>::maximum_smoothing);
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (3);

  const char *refinement_flags[]
    = 
    {
	  "1010222211110211222222220211211111111111111111111111111111111111",
	  "1212111111111111111111111111111111111111222222222222222211111111111111111101102201110210012211112001111111111111111111111111111111111111",
	  "1111111111111111111111111111111111111111111111111121111111111111111111111111111111111111111111111111111211111111111111111111111111021121222222221112112022222222021222011111111102102021111111111111110210102222111102121111111111211111222201012021111111111111",
	  "1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110112211111111112111110222222222222022202211212222210111111111111111121111111111221101111202221210222222222022222222220011222211111111222111111111111111112212111111112222000111111111111111111111112210022222021222221011220211111111111111111122111122220120222220212220011111111111",
	  "1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111121111111111111111111111111111111111211111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111100110011111111111111111111100010100000000000000001101010100000101000001101111100001000000011010000000011111111111111111010000000000001111111111111111111111111111111111111111110001100111111111111111111100000000000001111111111111111000001010010000000000000010100100100111100000101111111100100000000000101000001010000000011111111000101111111111110110010111111110000000011111111111111111111110011111111111110000000000000000001111100001111101000001111000010000010000011111111111111111100111100000000001000001111111101001111111100000100000000001111010111110000000111111111"
    };
  const unsigned int n_refinement_steps = sizeof(refinement_flags) /
					  sizeof(refinement_flags[0]);

  for (unsigned int i=0; i<n_refinement_steps; ++i)
    {
      Assert (triangulation.n_active_cells() ==
	      std::strlen(refinement_flags[i]),
	      ExcInternalError());

				       // set refinement flags from
				       // the string above
      unsigned int index = 0;
      for (Triangulation<2>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end(); ++cell, ++index)
	if (refinement_flags[i][index] == '2')
	  cell->set_refine_flag();
	else if (refinement_flags[i][index] == '1')
	  cell->set_coarsen_flag();

      triangulation.prepare_coarsening_and_refinement ();

				       // at this point, of course the
				       // mesh needs to be
				       // patch_level_1, but the flags
				       // also need to be set in such
				       // a way that it preserves
				       // patch_level_1. check
				       // this. the reason that the
				       // test failed is because flags
				       // were incorrectly set,
				       // resulting in a mesh that
				       // wasn't patch_level_1, and we
				       // got an assertion in the next
				       // call to prepare_c_and_r
				       // after setting flags next
				       // time around
      for (Triangulation<2>::cell_iterator
	     cell = triangulation.begin();
	   cell != triangulation.end(); ++cell)
	if (!cell->active())
	  {
	    unsigned int n_active_children = 0;
	    for (unsigned int i=0; i<cell->n_children(); ++i)
	      if (cell->child(i)->active())
		++n_active_children;
	    
	    Assert ((n_active_children == 0) || (n_active_children == cell->n_children()),
		    ExcInternalError());

					     // if there are no active
					     // children, then nothing
					     // to do here
	    if (n_active_children == 0)
	      continue;

					     // otherwise: all
					     // children are
					     // active. patch_level_1
					     // then means that either
					     // all of them need to be
					     // flagged for
					     // refinement, or
					     // none. check that this
					     // is so by making sure
					     // that they all have the
					     // same flags as child(0)
	    for (unsigned int i=1; i<cell->n_children(); ++i)
	      Assert(cell->child(i)->refine_flag_set()
		     ==
		     cell->child(0)->refine_flag_set(),
		     ExcInternalError());
	  }

      for (Triangulation<2>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end(); ++cell)
	if (cell->coarsen_flag_set ())
	  deallog << "Coarsen flag set: " << cell
		  << std::endl;
      
      
      triangulation.execute_coarsening_and_refinement ();
      
				       // and at this point, the
				       // mesh still needs to be
				       // patch_level_1, but no flags
				       // should remain set
      for (Triangulation<2>::cell_iterator
	     cell = triangulation.begin();
	   cell != triangulation.end(); ++cell)
	{
	  Assert ((cell->refine_flag_set() == false)
		  &&
		  (cell->coarsen_flag_set() == false),
		  ExcInternalError());
	  if (!cell->active())
	    Assert (cell_is_patch_level_1<2>(cell),
		    ExcInternalError());
	}

      
      deallog << triangulation.n_active_cells()
	      << " (on "
	      << triangulation.n_levels()
	      << " levels)"
	      << std::endl;
    }

  deallog << "OK" << std::endl;
}



int main () 
{
  std::ofstream logfile(logname);
  logfile.precision (3);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  try
    {
      test();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;

      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }

  return 0;
}
