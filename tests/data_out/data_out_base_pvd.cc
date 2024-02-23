// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// write the pvd primary record for parallel visualization through the
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
