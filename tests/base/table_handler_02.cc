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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


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
