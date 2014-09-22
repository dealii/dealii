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


// write the pvd master record for parallel visualization through the
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
    return std::vector<std::string>();
  }
};


template <int dim, int spacedim>
void check(std::ostream &out)
{
  std::vector<std::pair<double,std::string> > names(5);
  names[0] = std::make_pair(0,"x1");
  names[1] = std::make_pair(1,"x2");
  names[2] = std::make_pair(1e1,"x3");
  names[3] = std::make_pair(3.141,"d");
  names[4] = std::make_pair(42e19,"i");

  DataOutX x;
  x.write_pvd_record (out, names);
}



int main()
{
  std::ofstream logfile("output");
  check<2,2>(logfile);
}
