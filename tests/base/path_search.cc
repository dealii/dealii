// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
  for (unsigned int i=0; i<4; ++i)
    {
      in >> line;
      deallog << ' ' << line;
    }
  deallog << std::endl;
}
