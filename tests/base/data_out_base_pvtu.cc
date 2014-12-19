// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


// write the pvtu master record for parallel visualization through the
// vtu file format

#include "../tests.h"
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/logstream.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>

#include "patches.h"



std::vector<DataOutBase::Patch<2,2> > patches;
std::vector<std::string> names;

class DataOutX : public DataOutInterface<2,2>
{
  virtual
  const std::vector< ::DataOutBase::Patch<2,2> > &
  get_patches () const
  {
    return patches;
  }

  virtual
  std::vector<std::string>
  get_dataset_names () const
  {
    return names;
  }
};


template <int dim, int spacedim>
void check(std::ostream &out)
{
  names.resize (5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";

  std::vector<std::string> filenames = names;

  DataOutX x;
  x.write_pvtu_record (out, filenames);
}



int main()
{
  std::ofstream logfile("output");
  check<2,2>(logfile);
}
