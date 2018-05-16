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


// pad columns if we don't fill them in each iteration. also fill them below
// the last element if we create output


#include "../tests.h"
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>

#include <vector>
#include <string>


int
main ()
{
  initlog();

  TableHandler table;
  table.set_auto_fill_mode (true);

  std::string keys[3] = { "key1", "key2", "key3" };

  // fill rows 1 and 2 partially
  table.add_value(keys[0], 0);
  table.add_value(keys[0], 1);
  // now fill row 3 completely
  table.add_value(keys[0], 2);
  table.add_value(keys[1], 13);
  table.add_value(keys[2], std::string("a"));

  // now again fill row 4 partially
  table.add_value(keys[0], 1);

  // produce output. hope that row 4 is
  // completely padded
  table.write_text(deallog.get_file_stream());
}
