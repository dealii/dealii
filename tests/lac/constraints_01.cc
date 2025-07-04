// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
  cm.reinit(rel, rel);

  unsigned int inhoms[] = {8385,  8386,  8391,  17886, 17892, 17895, 18066,
                           18069, 18072, 18075, 18086, 18089, 18092, 18095,
                           18138, 18141, 18144, 18147, 18158, 18161, 18164};

  for (unsigned int i = 0; i < sizeof(inhoms) / sizeof(inhoms[0]); ++i)
    {
      deallog << inhoms[i] << std::endl;
      cm.constrain_dof_to_zero(inhoms[i]);
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
