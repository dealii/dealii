//-----------------------------------------------------------------------------
//    $Id: convergence_table_01.cc 20952 2010-04-06 15:02:46Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include "../tests.h"
#include <base/data_out_base.h>
#include <base/table_handler.h>
#include <base/convergence_table.h>
#include <base/logstream.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

// test the method evaluate_convergence_rates with argument reference_column_key, 

int main ()
{
  std::ofstream logfile("convergence_table_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  ConvergenceTable table;
  ConvergenceTable table_2;

  for (unsigned int i=1; i<=4; ++i)
  {
    table.add_value("error", 3.*i);
    table.add_value("key", i);
    table_2.add_value("error", std::pow(1./sqrt(i),4.));
    table_2.add_value("key", i);
  }

  table.evaluate_convergence_rates("error", "key", ConvergenceTable::reduction_rate);
  table_2.evaluate_convergence_rates("error", "key", ConvergenceTable::reduction_rate_log2);

                                 // output
  table.write_text(deallog.get_file_stream());
  table_2.write_text(deallog.get_file_stream());
}
