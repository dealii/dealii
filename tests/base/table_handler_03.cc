//-----------------------------------------------------------------------------
//    $Id$
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

// make sure we print something at all when an entry in the table corresponds
// to the empty string


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
  std::ofstream logfile("table_handler_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  TableHandler table;

  std::string keys[3] = { "key1", "key2", "key3" };

  table.add_value(keys[0], 0);
  table.add_value(keys[1], std::string(""));
  table.add_value(keys[2], std::string("a"));

  table.write_text(deallog.get_file_stream());
}
