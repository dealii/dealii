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

// Identical to data_out_base_gnuplot, but checks output of zero-dimensional
// features in a higher dimensional space

#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

#include "patches.h"

// Output data on repetitions of the unit hypercube

// define this as 1 to get output into a separate file for each testcase
#define SEPARATE_FILES 0


template <int dim, int spacedim>
void
check(DataOutBase::GnuplotFlags flags, std::ostream &out)
{
  const unsigned int np = 4;

  std::vector<DataOutBase::Patch<dim, spacedim>> patches(np);

  create_patches(patches);

  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vectors;
  DataOutBase::write_gnuplot(patches, names, vectors, flags, out);
}


template <int dim, int spacedim>
void
check_all(std::ostream &log)
{
#if SEPARATE_FILES == 0
  std::ostream &out = log;
#endif

  char                      name[100];
  const char *              format = "%d%d.gnuplot";
  DataOutBase::GnuplotFlags flags;
  for (unsigned int i = 0; i < 5; ++i)
    {
      sprintf(name, format, dim, spacedim, "");
#if SEPARATE_FILES == 1
      std::ofstream out(name);
#else
      out << "==============================\n"
          << name << "\n==============================\n";
#endif
      check<dim, spacedim>(flags, out);
    }
}

int
main()
{
  std::ofstream logfile("output");
  check_all<0, 1>(logfile);
  check_all<0, 2>(logfile);
  check_all<0, 3>(logfile);
}
