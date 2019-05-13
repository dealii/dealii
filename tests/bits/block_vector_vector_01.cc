// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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



// check existence of
// BlockVector<double>::BlockVector(BlockVector<float>). this conversion
// constructor was disabled previously altogether because of a compiler defect
// that did not honor the 'explicit' keyword on template constructors. this is
// now autoconf'ed.

#include <deal.II/lac/block_vector.h>

#include "../tests.h"


void
test(BlockVector<double> &v)
{
  for (unsigned int i = 0; i < v.size(); ++i)
    v(i) = i + 1.;
  BlockVector<float> w(v);

  AssertThrow(w == v, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      std::vector<types::global_dof_index> block_sizes(2, 50);
      BlockVector<double>                  v(block_sizes);
      test(v);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
