// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2017 by the deal.II authors
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
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>

#include <string>
#include <vector>

// test the method evaluate_convergence_rates with key column given by
// unsigned long long int

int
main()
{
  initlog();

  ConvergenceTable t;

  const unsigned long long int t1 = 100;
  const unsigned long long int t2 = 400;
  t.add_value("cells", t1);
  t.add_value("error", 0.1);
  t.add_value("cells", t2);
  t.add_value("error", 0.025);
  t.evaluate_convergence_rates(
    "error", "cells", ConvergenceTable::reduction_rate_log2, 2);
  t.write_text(deallog.get_file_stream());
}
