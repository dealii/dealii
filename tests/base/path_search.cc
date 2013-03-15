//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/base/path_search.h>
#include <deal.II/base/logstream.h>

int main()
{
  std::ofstream logfile("path_search/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  PathSearch::add_class("CC");

				   // Open with full debugging
  PathSearch cc("CC", 3);
  cc.add_path("../");
  cc.add_path("../scripts/", PathSearch::front);
  cc.add_suffix(".c");
  cc.add_suffix(".cc");

  deallog << cc.find("path_search.cc") << std::endl;
  deallog << cc.find("path_search") << std::endl;
  cc.add_suffix(".h", PathSearch::after_none);
  deallog << cc.find("tests") << std::endl;
  cc.show(deallog);
  
  PathSearch mesh("MESH", 3);
  mesh.add_path("../bits/");
  mesh.show(deallog);
  std::ifstream in(mesh.find("grid_in_msh_01.2d").c_str());
  std::string line;
  for (unsigned int i=0;i<4;++i)
    {
      in >> line;
      deallog << ' ' << line;
    }
  deallog << std::endl;
}
