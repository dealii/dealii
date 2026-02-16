// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test that AssertNothrow is printing what we expect.

#include "../tests.h"

int
main()
{
  initlog();
  deal_II_exceptions::disable_abort_on_exception();
  AssertNothrow(1 == 2, ExcInternalError());
}
