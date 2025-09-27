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


// test Utilities::encode_base64/decode_base64 on a few strings

#include <deal.II/base/utilities.h>

#include "../tests.h"

void
test(const std::string &input)
{
  deallog << "Encoding: " << input << std::endl;
  auto encoded = Utilities::encode_base64({input.begin(), input.end()});
  deallog << "Encoded: " << std::string(encoded.begin(), encoded.end())
          << std::endl;
  deallog << "Decoded: " << Utilities::decode_base64(encoded) << std::endl;
}

int
main()
{
  initlog();
  test("Hello");
  test("deal.II rocks!");
  test("1234");
}
