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
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/logstream.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

// test the method evaluate_convergence_rates with argument
// reference_column_key and data from various dimension where the reference
// data must be scaled by dimension

// data generated from a run of step-7 with Q2 elements for 1D, 2D and Q1
// elements for 3D

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  ConvergenceTable table_1;
  ConvergenceTable table_2;
  ConvergenceTable table_3;

  deallog << "Testing 1D data" << std::endl;
  table_1.add_value("n_cells", 16);
  table_1.add_value("error", 4.841e-04);
  table_1.add_value("n_cells", 24);
  table_1.add_value("error", 1.446e-04);
  table_1.add_value("n_cells", 32);
  table_1.add_value("error", 6.120e-05);
  table_1.add_value("n_cells", 48);
  table_1.add_value("error", 1.817e-05);
  table_1.add_value("n_cells", 64);
  table_1.add_value("error", 7.672e-06);
  table_1.add_value("n_cells", 96);
  table_1.add_value("error", 2.274e-06);
  table_1.add_value("n_cells", 128);
  table_1.add_value("error", 9.597e-07);
  table_1.set_scientific("error", true);
  table_1.evaluate_convergence_rates("error", "n_cells", ConvergenceTable::reduction_rate);
  table_1.evaluate_convergence_rates("error", "n_cells", ConvergenceTable::reduction_rate_log2, 1);
  table_1.write_text(deallog.get_file_stream());

  deallog << std::endl << "Testing 2d data" << std::endl;
  table_2.add_value("n_cells", 64);
  table_2.add_value("error", 7.755e-03);
  table_2.add_value("n_cells", 144);
  table_2.add_value("error", 2.326e-03);
  table_2.add_value("n_cells", 256);
  table_2.add_value("error", 9.969e-04);
  table_2.add_value("n_cells", 576);
  table_2.add_value("error", 2.987e-04);
  table_2.add_value("n_cells", 1024);
  table_2.add_value("error", 1.265e-04);
  table_2.add_value("n_cells", 2304);
  table_2.add_value("error", 3.759e-05);
  table_2.add_value("n_cells", 4096);
  table_2.add_value("error", 1.587e-05);
  table_2.set_scientific("error", true);
  table_2.evaluate_convergence_rates("error", "n_cells", ConvergenceTable::reduction_rate);
  table_2.evaluate_convergence_rates("error", "n_cells", ConvergenceTable::reduction_rate_log2, 2);
  table_2.write_text(deallog.get_file_stream());

  deallog << std::endl << "Testing 3d data" << std::endl;
  table_3.add_value("n_cells", 512);
  table_3.add_value("error", 5.791e-02);
  table_3.add_value("n_cells", 1728);
  table_3.add_value("error", 2.624e-02);
  table_3.add_value("n_cells", 4096);
  table_3.add_value("error", 1.480e-02);
  table_3.add_value("n_cells", 13824);
  table_3.add_value("error", 6.587e-03);
  table_3.add_value("n_cells", 32768);
  table_3.add_value("error", 3.708e-03);
  table_3.add_value("n_cells", 110592);
  table_3.add_value("error", 1.649e-03);
  table_3.add_value("n_cells", 262144);
  table_3.add_value("error", 9.275e-04);
  table_3.set_scientific("error", true);
  table_3.evaluate_convergence_rates("error", "n_cells", ConvergenceTable::reduction_rate);
  table_3.evaluate_convergence_rates("error", "n_cells", ConvergenceTable::reduction_rate_log2, 3);
  table_3.write_text(deallog.get_file_stream());


}
