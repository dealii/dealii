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
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
