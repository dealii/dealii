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


// test output of supercolumns


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
  table.set_auto_fill_mode(true);

  std::string keys[] = {"key1", "key2", "key3", "key4", "key5"};

  table.add_value(keys[0], 0.0);
  table.add_value(keys[0], 1.0);
  table.add_value(keys[1], 123.456);
  table.add_value(keys[2], std::string("abc"));
  table.add_value(keys[3], std::string("A"));
  table.add_value(keys[4], 123456789.0);

  table.add_column_to_supercolumn("key1", "short");
  table.add_column_to_supercolumn("key2", "short");
  table.add_column_to_supercolumn("key3", "very_long_supercolumn");
  table.add_column_to_supercolumn("key4", "very_long_supercolumn");

  table.write_text(deallog.get_file_stream(),
                   TableHandler::table_with_separate_column_description);

  deallog << std::endl;

  table.write_text(deallog.get_file_stream(), TableHandler::table_with_headers);
}
