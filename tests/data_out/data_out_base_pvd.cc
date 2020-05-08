// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// write the pvd master record for parallel visualization through the
// vtu file format

#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

#include "patches.h"


template <int dim, int spacedim>
void
check(std::ostream &out)
{
  std::vector<std::pair<double, std::string>> names(5);
  names[0] = std::make_pair(0, "x1");
  names[1] = std::make_pair(1, "x2");
  names[2] = std::make_pair(1e1, "x3");
  names[3] = std::make_pair(3.141, "d");
  names[4] = std::make_pair(42e19, "i");

  DataOutBase::write_pvd_record(out, names);
}



int
main()
{
  std::ofstream logfile("output");
  check<2, 2>(logfile);
}
