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

// Check that the description of a patterns works for all OutputStyles
// and number of PatternStyles.

#include <deal.II/base/parameter_handler.h>

#include <memory>

#include "../tests.h"

int
main()
{
  initlog();

  // one pattern
  {
    const auto &pattern = Patterns::Tuple(Patterns::Double());

    deallog << pattern.description(Patterns::PatternBase::OutputStyle::Machine)
            << '\n'
            << pattern.description(Patterns::PatternBase::OutputStyle::Text)
            << '\n'
            << pattern.description(Patterns::PatternBase::OutputStyle::LaTeX)
            << '\n'
            << std::endl;
  }

  // two patterns
  {
    const auto &pattern =
      Patterns::Tuple(Patterns::Double(), Patterns::Anything());

    deallog << pattern.description(Patterns::PatternBase::OutputStyle::Machine)
            << '\n'
            << pattern.description(Patterns::PatternBase::OutputStyle::Text)
            << '\n'
            << pattern.description(Patterns::PatternBase::OutputStyle::LaTeX)
            << '\n'
            << std::endl;
  }

  // three patterns
  {
    const auto &pattern = Patterns::Tuple(Patterns::Double(),
                                          Patterns::Anything(),
                                          Patterns::Bool());

    deallog << pattern.description(Patterns::PatternBase::OutputStyle::Machine)
            << '\n'
            << pattern.description(Patterns::PatternBase::OutputStyle::Text)
            << '\n'
            << pattern.description(Patterns::PatternBase::OutputStyle::LaTeX)
            << '\n'
            << std::endl;
  }
}
