//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include "../tests.h"
#include <base/path_search.h>
#include <base/logstream.h>

int main()
{
  std::ofstream logfile("path_search.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  PathSearch::add_class("CC");

				   // Open with full debugging
  PathSearch cc("CC", 3);
  cc.add_path(DEAL_II_PATH "/lac/source/");
  cc.add_path(DEAL_II_PATH "/lac/include/lac/", PathSearch::front);
  cc.add_suffix(".c");
  cc.add_suffix(".cc");

  cc.open("block_vector.cc");
  cc.open("block_vector");
  cc.add_suffix(".h", PathSearch::after_none);
  cc.open("block_vector");

  cc.show(deallog);
  
  PathSearch mesh("MESH", 3);
  mesh.show(deallog);
  std::istream& in = mesh.open("backstep");
  std::string line;
  for (unsigned int i=0;i<4;++i)
    {
      in >> line;
      deallog << ' ' << line;
    }
  deallog << std::endl;
}
