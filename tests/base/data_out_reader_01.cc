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


#include "../tests.h"
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/logstream.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "patches.h"

// test DataOutReader::merge

template <int dim, int spacedim>
void check()
{
  const unsigned int np = 1;

  std::vector<DataOutBase::Patch<dim, spacedim> > patches(np);

  create_patches(patches);

  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";
  std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > vectors;

  std::ostringstream old_data;
  DataOutBase::write_deal_II_intermediate(patches, names, vectors,
                                          DataOutBase::Deal_II_IntermediateFlags(),
                                          old_data);

  DataOutReader<dim,spacedim> data;
  {
    std::istringstream input(old_data.str());
    data.read (input);
  }
  DataOutReader<dim,spacedim> additional_data;
  {
    std::istringstream input(old_data.str());
    additional_data.read (input);
  }

  data.merge (additional_data);

  {
    std::ofstream out2( "outfile");
    data.write_deal_II_intermediate (out2);
  }

  cat_file("outfile");
  std::remove ("outfile");
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  check<1,1>();
  check<1,2>();
  check<2,2>();
  check<2,3>();
  check<3,3>();
}
