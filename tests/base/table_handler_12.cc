//-----------------------------------------------------------------------------
//    $Id: table_handler_07.cc -1   $
//    Version: $Name$
//
//    Copyright (C) 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// verify that the flags we set for precision when printing stuff from
// a table do not affect the precision flags set for the stream to
// which we print


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
  std::ofstream logfile("table_handler_12/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

				   // set precision on the output
				   // stream to 4 digits
  logfile << std::setprecision(4);

				   // but then set precision on the
				   // table output to 2
  TableHandler table;
  table.add_value("key", 0.123456789);
  table.set_precision ("key", 2);

				   // now output the table...
  table.write_text(logfile);
				   // ...and then output some other
				   // number, hopefully with 4 digits
				   // of precision
  logfile << 0.123456789 << std::endl;
}
