// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"


// bug in AffineConstraints<double>

// index=18466 line_index=652 lines_cache[line_index]=919940456 lines.size()=21

void
test()
{
  IndexSet      rel;
  std::ifstream f(SOURCE_DIR "/constraints_01/is.23");
  rel.read(f);

  AffineConstraints<double> cm;
  cm.clear();
  cm.reinit(rel);

  unsigned int inhoms[] = {8385,  8386,  8391,  17886, 17892, 17895, 18066,
                           18069, 18072, 18075, 18086, 18089, 18092, 18095,
                           18138, 18141, 18144, 18147, 18158, 18161, 18164};

  for (unsigned int i = 0; i < sizeof(inhoms) / sizeof(inhoms[0]); ++i)
    {
      deallog << inhoms[i] << std::endl;
      cm.add_line(inhoms[i]);
      cm.set_inhomogeneity(inhoms[i], 1.0);
    }

  cm.print(deallog.get_file_stream());

  bool is = cm.is_inhomogeneously_constrained(18466);
  deallog << "constraint 18466 inhom? " << is << std::endl;
  Assert(!is, ExcInternalError());
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test();

  deallog << "OK" << std::endl;
}
