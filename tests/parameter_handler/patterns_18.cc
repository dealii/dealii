// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2023 by the deal.II authors
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

// test that lists of lists are parsed correctly
#include <deal.II/base/patterns.h>

#include "../tests.h"

int
main()
{
  initlog();

  // Example parameter pattern for a list of lists
  const std::string description =
    "[List of <[List of <[Double -MAX_DOUBLE...MAX_DOUBLE (inclusive)]> of length 3...3 (inclusive) separated by < >]> of length 0...4294967295 (inclusive)]";
  // The above pattern describes a list of lists of 3 doubles, e.g.:
  // "1 2 3,1.2 3.2 1.2" is valid
  // "1.2,1.3,1.4,1.6 1.2,1.3,1.0,1.5,6.5,5.1 7.6,13.5,1.5" is not valid

  // Parse the description in deal.ii
  std::unique_ptr<dealii::Patterns::PatternBase> pattern =
    dealii::Patterns::pattern_factory(description);

  // Write the original and parsed description
  deallog << "Original description: " << std::endl
          << description << std::endl
          << std::endl
          << "Parsed description: " << std::endl
          << pattern->description(dealii::Patterns::PatternBase::Machine)
          << std::endl;

  // First check, should be a match
  deallog << "First check: \"\" is "
          << (pattern->match("") ? "a match" : "not a match") << std::endl;

  // Second check, should be a match
  deallog << "Second check: \"1 2 3,1.2 3.2 1.2\" is "
          << (pattern->match("1 2 3,1.2 3.2 1.2") ? "a match" : "not a match")
          << std::endl;

  // Third check, should not be a match
  deallog
    << "Third check: \"1.2,1.3,1.4,1.6 1.2,1.3,1.0,1.5,6.5,5.1 7.6,13.5,1.5\" is "
    << (pattern->match("1.2,1.3,1.4,1.6 1.2,1.3,1.0,1.5,6.5,5.1 7.6,13.5,1.5") ?
          "a match" :
          "not a match")
    << std::endl;
}
