// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// make sure that block vector iterator allows reading and writing correctly

#include <deal.II/lac/block_vector.h>

#include "../tests.h"


void
test()
{
  BlockVector<double> v(2, 1);
  v(0) = 1;
  v(1) = 2;

  // first check reading through a const
  // iterator
  {
    BlockVector<double>::const_iterator i = v.begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // same, but create iterator in a different
  // way
  {
    BlockVector<double>::const_iterator i =
      const_cast<const BlockVector<double> &>(v).begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // read through a read-write iterator
  {
    BlockVector<double>::iterator i = v.begin();
    AssertThrow(*i == 1, ExcInternalError());
    ++i;
    AssertThrow(*i == 2, ExcInternalError());
  }

  // write through a read-write iterator
  {
    BlockVector<double>::iterator i = v.begin();

    *i = 2;
    ++i;
    *i = 3;
  }

  // and read again
  {
    BlockVector<double>::iterator i = v.begin();
    AssertThrow(*i == 2, ExcInternalError());
    ++i;
    AssertThrow(*i == 3, ExcInternalError());
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      test();
    }
  catch (const std::exception &exc)
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
