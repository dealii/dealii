// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Make sure that Utilities::pack/unpack can be used on objects that
// are empty, and that the resulting packed string has length zero.


#include <deal.II/base/point.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <tuple>

#include "../tests.h"



void
test()
{
  using Empty = std::tuple<>;

  Empty e;
  deallog << "Size = " << sizeof(e) << std::endl;

  // Pack the object and make sure the packed length is zero
  const std::vector<char> packed = Utilities::pack(e, false);
  deallog << "Packed size = " << packed.size() << std::endl;

  // Then unpack again. The two should be the same -- though one may
  // question what equality of objects of size zero might mean -- and
  // we should check so:
  Empty e2 = Utilities::unpack<Empty>(packed, false);
  Assert(e2 == e, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
