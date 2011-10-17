//----------------------------  table_1.cc  ---------------------------
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
//----------------------------  table_1.cc  ---------------------------

// check serialization for TableHandler

#include "serialization.h"
#include <deal.II/base/table_handler.h>
#include <boost/serialization/vector.hpp>

namespace dealii
{
  bool compare (const TableHandler &t1,
		const TableHandler &t2)
  {
    std::ostringstream o1, o2;
    t1.write_tex (o1);
    t1.write_text (o1);
    t2.write_tex (o2);
    t2.write_text (o2);

    return (o1.str() == o2.str());
  }
}


void test ()
{
  TableHandler t1, t2;
  std::string keys[4] = { "key1", "key2", "key3", "key4" };


  for (unsigned int i=0; i<10; ++i)
  {
    t1.add_value(keys[(0+i)%4], i);
    t1.add_value(keys[(1+i)%4], sqrt(i));
    t1.add_value(keys[(2+i)%4], 'a'+i);
    t1.add_value(keys[(3+i)%4], std::string("abc-")+"0123456789"[i]);
  }
  t1.set_tex_table_caption("This is a caption text with \\LaTeX{} symbols");

  verify (t1, t2);
}


int main ()
{
  std::ofstream logfile("table_handler/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
