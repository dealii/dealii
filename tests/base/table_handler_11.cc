// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


// test output of supercolumns


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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  TableHandler table;
  table.set_auto_fill_mode (true);

  std::string keys[] = { "key1", "key2", "key3", "key4", "key5" };

  table.add_value(keys[0], 0.0);
  table.add_value(keys[0], 1.0);
  table.add_value(keys[1], 123.456);
  table.add_value(keys[2], std::string("abc"));
  table.add_value(keys[3], std::string("A"));
  table.add_value(keys[4], 123456789.0);

  table.add_column_to_supercolumn("key1","short");
  table.add_column_to_supercolumn("key2","short");
  table.add_column_to_supercolumn("key3","very_long_supercolumn");
  table.add_column_to_supercolumn("key4","very_long_supercolumn");

  table.write_text(deallog.get_file_stream(),
                   TableHandler::table_with_separate_column_description);

  deallog << std::endl;

  table.write_text(deallog.get_file_stream(),
                   TableHandler::table_with_headers);
}
