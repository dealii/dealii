//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test IndexSet::add_indices(IndexSet)

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include <deal.II/base/index_set.h>

void testor(IndexSet & a,
	    IndexSet & other,
	    unsigned int offset,
	    bool verbose)
{
  IndexSet merged(a);

  merged.add_indices(other, offset);

  if (verbose)
    {
      deallog << "Original index set: " << std::endl;
      a.print(deallog);
      deallog << "other index set: " << std::endl;
      other.print(deallog);
      deallog << "merged index set: " << std::endl;
      merged.print(deallog);
    }

  for (unsigned int i=0;i<merged.size();++i)
    {
      Assert(
	merged.is_element(i)
	==
	(a.is_element(i) || (i>=offset && other.is_element(i-offset))),
	ExcInternalError());
    }
}




void test()
{
  const int size = 10;
  
  IndexSet empty(size);
  IndexSet id(size);

  id.add_index(3);
  id.add_index(4);
  id.add_index(7);
  
  deallog << "* add empty: " << std::endl;
  testor(id, empty, 2, true);

  deallog << "* add self: " << std::endl;
  testor(id, id, 2, true);

  deallog << "* add id2: " << std::endl;
  IndexSet id2(size);
  id2.add_index(0);
  id2.add_index(2);
  id2.add_index(3);
  testor(id, id2, 3, true);
}






int main()
{
  std::ofstream logfile("index_set_20_offset/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
