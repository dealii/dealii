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


// we test the performance of writing a large table

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

using namespace dealii;


#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>


int main ()
{
  TimerOutput timer (std::cout, TimerOutput::summary, TimerOutput::cpu_times);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  TableHandler table;

  std::string keys[] = { "key1", "key2", "key3", "key4", "key5", "key6", "key7", "key8", "key9", "key10", "key11", "key12", "key13", "key14", "key15"};

  unsigned int n_keys = 15;
  unsigned int n_rows = 40000;

  for (unsigned int j=0; j<n_rows; ++j)
    {
      table.add_value("begin", std::string("this is some text"));
      for (unsigned int i=0; i<n_keys; ++i)
        {
          table.add_value(keys[i], j*1.0+i/100.0);
        }
    }

  timer.enter_section("write");
  {
    std::ofstream data("datatable.txt");
    table.write_text(data);
//       TableHandler::table_with_separate_column_description);
  }
  timer.exit_section("write");
}
