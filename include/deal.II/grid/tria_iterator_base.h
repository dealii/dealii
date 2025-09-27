// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tria_iterator_base_h
#define dealii_tria_iterator_base_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace in which an enumeration is declared that denotes the states in
 * which an iterator can be in.
 *
 * @ingroup Iterators
 */
namespace IteratorState
{
  /**
   * The three states an iterator can be in: valid, past-the-end and invalid.
   */
  enum IteratorStates
  {
    /// Iterator points to a valid object
    valid,
    /// Iterator reached end of container
    past_the_end,
    /// Iterator is invalid, probably due to an error
    invalid
  };
} // namespace IteratorState



DEAL_II_NAMESPACE_CLOSE

#endif
