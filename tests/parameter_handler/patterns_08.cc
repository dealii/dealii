// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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
