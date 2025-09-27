// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that one can output physical units along with output
// fields. This test is for scalar quantities and checks the
// pvtu writer.


#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

#include "patches.h"

// Output data on repetitions of the unit hypercube

template <int dim, int spacedim>
void
check(DataOutBase::VtkFlags flags, std::ostream &out)
{
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
  DataOutBase::write_pvtu_record(out, {"abc.vtu"}, names, vectors, flags);
}


template <int dim, int spacedim>
void
check_all(std::ostream &log)
{
  DataOutBase::VtkFlags flags;

  // Set the units of some but not all output quantities
  flags.physical_units["x1"] = "kg/s";
  flags.physical_units["x3"] = "N/m";
  flags.physical_units["i"]  = "J/K";

  log << "==============================" << std::endl
      << dim << spacedim << ".vtu" << std::endl
      << "==============================" << std::endl;
  check<dim, spacedim>(flags, log);
}

int
main()
{
  initlog();

  check_all<1, 1>(deallog.get_file_stream());
  check_all<1, 2>(deallog.get_file_stream());
  check_all<2, 2>(deallog.get_file_stream());
  check_all<2, 3>(deallog.get_file_stream());
  check_all<3, 3>(deallog.get_file_stream());
}
