// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  // create a pattern and let it
  // output its description
  std::vector<std::unique_ptr<Patterns::PatternBase>> ps;
  ps.push_back(std::make_unique<Patterns::Integer>());
  ps.push_back(std::make_unique<Patterns::Double>());
  ps.push_back(std::make_unique<Patterns::Anything>());

  Patterns::Tuple   pattern(ps, ";");
  const std::string desc = pattern.description();

  // now let the same class re-create
  // a pattern object from the
  // description and verify that the
  // result is the same as what we
  // started out with
  std::unique_ptr<Patterns::Tuple> pattern2 = Patterns::Tuple::create(desc);

  AssertThrow(pattern2 != nullptr, ExcInternalError());
  AssertThrow(desc == pattern2->description(), ExcInternalError());

  deallog << desc << std::endl;
}
