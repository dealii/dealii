// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check serialization for TableHandler

#include <deal.II/base/table_handler.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"

namespace dealii
{
  bool
  compare(const TableHandler &t1, const TableHandler &t2)
  {
    std::ostringstream o1, o2;
    t1.write_tex(o1);
    t1.write_text(o1);
    t2.write_tex(o2);
    t2.write_text(o2);

    return (o1.str() == o2.str());
  }
} // namespace dealii


void
test()
{
  TableHandler t1, t2;
  std::string  keys[4] = {"key1", "key2", "key3", "key4"};


  for (unsigned int i = 0; i < 10; ++i)
    {
      t1.add_value(keys[(0 + i) % 4], i);
      t1.add_value(keys[(1 + i) % 4], sqrt(i));
      t1.add_value(keys[(2 + i) % 4], 'a' + i);
      t1.add_value(keys[(3 + i) % 4], std::string("abc-") + "0123456789"[i]);
    }
  t1.set_tex_table_caption("This is a caption text with \\LaTeX{} symbols");

  verify(t1, t2);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
