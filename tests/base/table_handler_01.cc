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


#include "../tests.h"
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/logstream.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

// test the method set_tex_table_caption
// it creates a caption for the whole table

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  TableHandler table;

  for (unsigned int i=1; i<=10; ++i)
    {
      table.add_value("numbers", i);
      table.add_value("squares", i*i);
      table.add_value("square roots", sqrt(i));
    }

  table.set_tex_table_caption("This is a caption text with \\LaTeX{} symbols");

  // output
  table.write_tex(deallog.get_file_stream());
}
