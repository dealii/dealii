// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test functions in namespace Utilities

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  deallog << Utilities::int_to_string(42, 4) << std::endl;
  deallog << Utilities::int_to_string(42) << std::endl;
  deallog << Utilities::needed_digits(424) << std::endl;
  deallog << Utilities::needed_digits(0) << std::endl;
  deallog << Utilities::needed_digits(1) << std::endl;
  deallog << Utilities::needed_digits(2) << std::endl;
  deallog << Utilities::needed_digits(10) << std::endl;
  deallog << Utilities::needed_digits(1000000000) << std::endl;
  deallog << Utilities::string_to_int(" 413 ") << std::endl;

  std::vector<std::string> v;
  v.push_back("1");
  v.push_back(" -12");
  v.push_back("+125 ");
  AssertThrow(Utilities::string_to_int(v).size() == 3, ExcInternalError());
  deallog << Utilities::string_to_int(v)[0] << std::endl;
  deallog << Utilities::string_to_int(v)[1] << std::endl;
  deallog << Utilities::string_to_int(v)[2] << std::endl;

  {
    const char *p = "alpha, beta, gamma ";
    AssertThrow(Utilities::split_string_list(p).size() == 3,
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p)[0] == "alpha",
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p)[1] == "beta",
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p)[2] == "gamma",
                ExcInternalError());
  }

  {
    const char *p = "alpha; beta; gamma ";
    AssertThrow(Utilities::split_string_list(p, ';').size() == 3,
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p, ';')[0] == "alpha",
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p, ';')[1] == "beta",
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p, ';')[2] == "gamma",
                ExcInternalError());
  }

  {
    const char *p = "alpha;; beta;; gamma ";
    AssertThrow(Utilities::split_string_list(p, ";;").size() == 3,
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p, ";;")[0] == "alpha",
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p, ";;")[1] == "beta",
                ExcInternalError());
    AssertThrow(Utilities::split_string_list(p, ";;")[2] == "gamma",
                ExcInternalError());
  }

  deallog << Utilities::generate_normal_random_number(13, 44) << ' ';
  deallog << Utilities::generate_normal_random_number(13, 44) << ' ';
  deallog << Utilities::generate_normal_random_number(13, 44) << ' ';
  deallog << std::endl;
}



int
main()
{
  initlog();

  test();
}
