// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// test TableHandler::declare_column()


#include "../tests.h"
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/logstream.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>


int main ()
{
  initlog ();

  TableHandler table;

  // declare columns either directly, or by adding a value
  table.declare_column ("key1");
  table.add_value("key2", 0);
  table.declare_column ("key3");

  // finish that row (even though there is something only in column 2)
  // and output the whole shebang
  table.start_new_row ();

  table.write_text(deallog.get_file_stream(),
                   TableHandler::org_mode_table);
}
