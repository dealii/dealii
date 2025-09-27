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


// test output of tables with columns in a variety of data types


#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>

#include <string>
#include <vector>

#include "../tests.h"


int
main()
{
  initlog();

  TableHandler table;

  // have keys that we rotate in each row to
  // make sure that individual entries in the
  // same column really have different data
  // types
  std::string keys[4] = {"key1", "key2", "key3", "key4"};


  for (unsigned int i = 0; i < 10; ++i)
    {
      table.add_value(keys[(0 + i) % 4], i);
      table.add_value(keys[(1 + i) % 4], sqrt(i));
      table.add_value(keys[(2 + i) % 4], 'a' + i);
      table.add_value(keys[(3 + i) % 4], std::string("abc-") + "0123456789"[i]);
    }

  // output
  table.write_text(deallog.get_file_stream());
}
