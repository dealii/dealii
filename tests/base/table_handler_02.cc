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

// test output of tables with columns in a variety of data types


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
  std::ofstream logfile("table_handler_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  TableHandler table;

				   // have keys that we rotate in each row to
				   // make sure that individual entries in the
				   // same column really have different data
				   // types
  std::string keys[4] = { "key1", "key2", "key3", "key4" };


  for (unsigned int i=0; i<10; ++i)
  {
    table.add_value(keys[(0+i)%4], i);
    table.add_value(keys[(1+i)%4], sqrt(i));
    table.add_value(keys[(2+i)%4], 'a'+i);
    table.add_value(keys[(3+i)%4], std::string("abc-")+"0123456789"[i]);
  }

                                 // output
  table.write_text(deallog.get_file_stream());
}
