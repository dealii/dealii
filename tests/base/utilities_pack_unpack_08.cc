// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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
