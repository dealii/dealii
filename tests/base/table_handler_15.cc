// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test TableHandler::declare_column()


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

  // declare columns either directly, or by adding a value
  table.declare_column("key1");
  table.add_value("key2", 0);
  table.declare_column("key3");

  // finish that row (even though there is something only in column 2)
  // and output the whole shebang
  table.start_new_row();

  table.write_text(deallog.get_file_stream(), TableHandler::org_mode_table);
}
