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
