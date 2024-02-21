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


// write the pvtu primary record for parallel visualization through the
// vtu file format

#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

#include "patches.h"



std::vector<DataOutBase::Patch<2, 2>> patches;
std::vector<std::string>              names;

class DataOutX : public DataOutInterface<2, 2>
{
  virtual const std::vector<::DataOutBase::Patch<2, 2>> &
  get_patches() const
  {
    return patches;
  }

  virtual std::vector<std::string>
  get_dataset_names() const
  {
    return names;
  }
};


template <int dim, int spacedim>
void
check(std::ostream &out)
{
  names.resize(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";

  std::vector<std::string> filenames = names;

  DataOutX x;
  x.write_pvtu_record(out, filenames);
}



int
main()
{
  std::ofstream logfile("output");
  check<2, 2>(logfile);
}
