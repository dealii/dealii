//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test IndexSet::fill_index_vector

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include <deal.II/base/index_set.h>


void test ()
{
  IndexSet is1 (100);

				   // randomly add 90 elements to each
				   // set, some of which may be
				   // repetitions of previous ones
  for (unsigned int i=0; i<9*is1.size()/10; ++i)
    is1.add_index (rand() % is1.size());

  std::vector<unsigned int> indices;
  is1.fill_index_vector (indices);

  deallog << "Original index set: " << std::endl;
  is1.print(deallog);

  deallog << "List of indices: " << std::endl;
  for (unsigned int i=0; i<indices.size(); i++)
    deallog << indices[i] << ' ';
  deallog << std::endl;

  for (unsigned int i=0; i<indices.size(); i++)
    Assert(is1.index_within_set(indices[i])==i, ExcInternalError());

  deallog << "OK" << std::endl;
}




int main()
{
  std::ofstream logfile("index_set_19/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
