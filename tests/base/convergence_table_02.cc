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

// test the method evaluate_all_convergence_rates with argument
// reference_column_key,

int
main()
{
  initlog();

  ConvergenceTable table;

  for (unsigned int i = 1; i <= 4; ++i)
    {
      table.add_value("error", 1. / sqrt(i));
      table.add_value("error2", std::pow(1. / sqrt(i), 2.));
      table.add_value("error3", std::pow(1. / sqrt(i), 3.));
      table.add_value("error4", std::pow(1. / sqrt(i), 4.));
      table.add_value("key", i);
    }
  table.omit_column_from_convergence_rate_evaluation("key");
  table.evaluate_all_convergence_rates("key",
                                       ConvergenceTable::reduction_rate_log2);

  // output
  table.write_text(deallog.get_file_stream());
}
