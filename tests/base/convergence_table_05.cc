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


#include <deal.II/base/convergence_table.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>

#include <string>
#include <vector>

#include "../tests.h"

// test the method evaluate_convergence_rates with key column given by
// 64 bit dof indices

int
main()
{
  initlog();

  ConvergenceTable t;

  const std::uint64_t t1 = 100;
  const std::uint64_t t2 = 400;
  t.add_value("cells", t1);
  t.add_value("error", 0.1);
  t.add_value("cells", t2);
  t.add_value("error", 0.025);
  t.evaluate_convergence_rates("error",
                               "cells",
                               ConvergenceTable::reduction_rate_log2,
                               2);
  t.write_text(deallog.get_file_stream());
}
