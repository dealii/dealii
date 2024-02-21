// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// document crash in deallog related to missing newline


#include <limits>

#include "../tests.h"


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  try
    {
      initlog();
      deallog << "OK" << std::endl;
      deallog << "no newline here!";
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}
