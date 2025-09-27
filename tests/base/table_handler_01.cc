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


#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>

#include <string>
#include <vector>

#include "../tests.h"

// test the method set_tex_table_caption
// it creates a caption for the whole table

int
main()
{
  initlog();

  TableHandler table;

  for (unsigned int i = 1; i <= 10; ++i)
    {
      table.add_value("numbers", i);
      table.add_value("squares", i * i);
      table.add_value("square roots", sqrt(i));
    }

  table.set_tex_table_caption("This is a caption text with \\LaTeX{} symbols");

  // output
  table.write_tex(deallog.get_file_stream());
}
