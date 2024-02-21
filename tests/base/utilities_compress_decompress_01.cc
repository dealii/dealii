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


// test Utilities::compress/decompress on a string

#include <deal.II/base/utilities.h>

#include "../tests.h"

int
main()
{
  initlog();
  std::string input        = "deal.II Rocks!";
  auto        compressed   = Utilities::compress(input);
  auto        decompressed = Utilities::decompress(compressed);
  deallog << (decompressed == input ? "OK" : "NOT OK") << std::endl;
}
