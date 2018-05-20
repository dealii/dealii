// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
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

// make sure that block vector iterator allows reading and writing correctly

#include "../tests.h"
#include <deal.II/lac/block_vector.h>

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
    BlockVector<double>::const_iterator i
      = const_cast<const BlockVector<double>&>(v).begin();
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
  catch(std::exception& exc)
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
  catch(...)
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
