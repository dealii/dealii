//----------------------------  sparsity_tools_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparsity_tools_01.cc  ---------------------------


// apply SparsityTools::reorder_Cuthill_McKee to a graph that consists
// of two or more non-connected parts. the reordering algorithm used
// to trip over that

#include "../tests.h"
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <iomanip>
#include <fstream>

int main ()
{
  std::ofstream logfile("sparsity_tools_01/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  SparsityPattern sp (4,4,4);
  for (unsigned int i=0; i<4; ++i)
    sp.add (i,i);

				   // create a graph with components
				   // 0,2 and 1,3 that are
				   // disconnected
  sp.add (0,2);
  sp.add (2,0);

  sp.add (1,3);
  sp.add (3,1);

  sp.compress ();

				   // now find permutation
  std::vector<unsigned int> permutation(4);
  SparsityTools::reorder_Cuthill_McKee (sp, permutation);

  for (unsigned int i=0; i<permutation.size(); ++i)
    deallog << permutation[i] << std::endl;
}



