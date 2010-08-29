//-----------------------------------------------------------------------------
//    $Id$
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

// test the method evaluate_all_convergence_rates with argument reference_column_key, 

int main ()
{
  std::ofstream logfile("convergence_table_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  ConvergenceTable table;

  for (unsigned int i=1; i<=4; ++i)
  {
    table.add_value("error", 1./sqrt(i));
    table.add_value("error2", std::pow(1./sqrt(i),2.));
    table.add_value("error3", std::pow(1./sqrt(i),3.));
    table.add_value("error4", std::pow(1./sqrt(i),4.));
    table.add_value("key", i);
  }
  table.omit_column_from_convergence_rate_evaluation ( "key" );
  table.evaluate_all_convergence_rates("key", ConvergenceTable::reduction_rate_log2);

                                 // output
  table.write_text(deallog.get_file_stream());
}
