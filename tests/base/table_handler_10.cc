//-----------------------------------------------------------------------------
//    $Id: table_handler_06.cc 24924 2012-01-25 12:35:17Z kormann $
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// test changing header style after adding entries


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
  std::ofstream logfile("table_handler_10/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  TableHandler table;
  table.set_auto_fill_mode (true);

  std::string keys[3] = { "key1", "key2", "key3" };

				   // fill rows 1 and 2 partially
  table.add_value(keys[0], 0.0);
  table.add_value(keys[0], 1.0);
				   // now fill row 3 completely
  table.add_value(keys[0], 2.0);
  table.add_value(keys[1], 123.456);
  table.add_value(keys[2], std::string("abc"));

				   // now again fill row 4 partially
  table.add_value(keys[0], 1.0);

				   // produce output. hope that row 4 is
				   // completely padded
  table.write_text(deallog.get_file_stream(),
		   TableHandler::table_with_separate_column_description);

  table.set_scientific(keys[0], true);
  table.set_precision(keys[1], 1);

  table.write_text(deallog.get_file_stream(),
		   TableHandler::table_with_separate_column_description);
  
}
