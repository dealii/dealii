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


// verify that the flags we set for precision when printing stuff from
// a table do not affect the precision flags set for the stream to
// which we print


#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>

#include <string>
#include <vector>

#include "../tests.h"


int
main()
{
  initlog();

  // set precision on the output
  // stream to 4 digits
  deallog.get_file_stream() << std::setprecision(4);

  // but then set precision on the
  // table output to 2
  TableHandler table;
  table.add_value("key", 0.123456789);
  table.set_precision("key", 2);

  // now output the table...
  table.write_text(deallog.get_file_stream());
  // ...and then output some other
  // number, hopefully with 4 digits
  // of precision
  deallog.get_file_stream() << 0.123456789 << std::endl;
}
