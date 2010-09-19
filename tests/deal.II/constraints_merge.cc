//----------------------------  constraints.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2003, 2004, 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  constraints.cc  ---------------------------


// merge and print a bunch of ConstrainMatrices

#include "../tests.h"
#include <dofs/dof_handler.h>
#include <grid/tria.h>
#include <fe/mapping_q1.h>
#include <fe/fe_q.h>
#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <lac/sparse_matrix.h>
#include <base/parameter_handler.h>
#include <dofs/dof_accessor.h>
#include <lac/constraint_matrix.h>
#include <dofs/dof_tools.h>
#include <grid/grid_out.h>
#include <base/logstream.h>

#include <fstream>
#include <cmath>
#include <cstdlib>


std::ofstream logfile("constraints_merge/output");


void merge_check ()
{
  deallog << "Checking ConstraintMatrix::merge" << std::endl;

				   // check twice, once with closed
				   // objects, once with open ones
  for (unsigned int run=0; run<2; ++run)
    {
      deallog << "Checking with " << (run == 0 ? "open" : "closed")
	      << " objects" << std::endl;
      
				       // check that the `merge' function
				       // works correctly
      ConstraintMatrix c1, c2;

				       // enter simple line
      c1.add_line (0);
      c1.add_entry (0, 11, 1.);
				       // add more complex line
      c1.add_line (1);
      c1.add_entry (1, 3, 0.5);
      c1.add_entry (1, 4, 0.5);
  
				       // fill second constraints
				       // object with one trivial line
				       // and one which further
				       // constrains one of the
				       // entries in the first object
      c2.add_line (10);
      c2.add_entry (10, 11, 1.);
      c2.add_line (3);
      c2.add_entry (3, 12, 0.25);
      c2.add_entry (3, 13, 0.75);

				       // in one of the two runs,
				       // close the objects
      if (run == 1)
	{
	  c1.close ();
	  c2.close ();
	};

				       // now merge the two and print the
				       // results
      c1.merge (c2);
      c1.print (logfile);
    };
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  merge_check ();
}

