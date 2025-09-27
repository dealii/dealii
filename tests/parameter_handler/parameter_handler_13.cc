// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check the Patterns::Map pattern

#include <deal.II/base/parameter_handler.h>

#include "../tests.h"

void
check(const char *p)
{
  ParameterHandler prm;
  prm.declare_entry("test_13",
                    "-1:a, 0:b, 1:c",
                    Patterns::Map(Patterns::Integer(-1, 1),
                                  Patterns::Selection("a|b|c"),
                                  2,
                                  3));

  std::ifstream in(p);
  prm.parse_input(in);

  deallog << "test_13=" << prm.get("test_13") << std::endl;

  const std::vector<std::string> split_entries =
    Utilities::split_string_list(prm.get("test_13"), ',');
  for (const std::string &entry : split_entries)
    {
      const std::vector<std::string> parts =
        Utilities::split_string_list(entry, ':');
      const std::string key   = parts[0];
      const std::string value = parts[1];
      deallog << " found key = '" << key << "' value = '" << value << "'"
              << std::endl;
    }
}


int
main()
{
  initlog();

  check(SOURCE_DIR "/prm/parameter_handler_13.prm");

  return 0;
}
