// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
