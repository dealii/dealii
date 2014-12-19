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
#include <sstream>
#include <string>
#include <stdio.h>

#include "patches.h"

// Output data on repetitions of the unit hypercube

template <int dim, int spacedim>
void check(DataOutBase::TecplotFlags flags,
           std::ostream &out)
{
  const unsigned int np = 4;

  std::vector<DataOutBase::Patch<dim, spacedim> > patches(np);

  create_patches(patches);

  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";
  std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > vectors;
  DataOutBase::write_tecplot_binary(patches, names, vectors, flags, out);
}


template<int dim, int spacedim>
void check_all()
{
  char name[100];
  DataOutBase::TecplotFlags flags;
  if (true)
    {
      sprintf(name, "%d%d.tecplot", dim, spacedim);
      flags.tecplot_binary_file_name=name;

      check<dim,spacedim>(flags, deallog.get_file_stream());
    }
}

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  check_all<1,1>();
  check_all<1,2>();
  check_all<2,2>();
  check_all<2,3>();
  check_all<3,3>();
}
