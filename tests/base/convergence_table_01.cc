// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/convergence_table.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>

#include <string>
#include <vector>

#include "../tests.h"

// test the method evaluate_convergence_rates with argument
// reference_column_key, assuming data coming from 2d

int
main()
{
  initlog();

  ConvergenceTable table;
  ConvergenceTable table_2;

  for (unsigned int i = 1; i <= 4; ++i)
    {
      table.add_value("error", 3. * i);
      table.add_value("key", i);
      table_2.add_value("error", std::pow(1. / sqrt(i), 4.));
      table_2.add_value("key", i);
    }

  table.evaluate_convergence_rates("error",
                                   "key",
                                   ConvergenceTable::reduction_rate);
  table_2.evaluate_convergence_rates("error",
                                     "key",
                                     ConvergenceTable::reduction_rate_log2);

  // output
  table.write_text(deallog.get_file_stream());
  table_2.write_text(deallog.get_file_stream());
}
