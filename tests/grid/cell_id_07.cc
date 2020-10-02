// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// testing CellId from / to_string

#include <deal.II/grid/cell_id.h>

#include "../tests.h"

int
main()
{
  initlog();

  std::vector<CellId> ids;

  ids.emplace_back(CellId(0, {0}));
  ids.emplace_back(CellId(0, {1}));
  ids.emplace_back(CellId(0, {0, 0}));
  ids.emplace_back(CellId(0, {1, 2, 0}));
  ids.emplace_back(CellId(23, {}));

  for (const auto &id : ids)
    {
      std::string s = id.to_string();
      deallog << s << '\n';
      AssertThrow(CellId(s) == id, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}
