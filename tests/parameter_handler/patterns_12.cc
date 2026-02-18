// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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
