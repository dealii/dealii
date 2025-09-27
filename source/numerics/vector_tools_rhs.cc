// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/numerics/vector_tools_rhs.templates.h>

#include <set>


DEAL_II_NAMESPACE_OPEN


namespace VectorTools
{
  // separate implementation for 1d because otherwise we get linker errors since
  // (hp::)FEFaceValues<1> is not compiled
  template <>
  void
  create_boundary_right_hand_side(const Mapping<1, 1> &,
                                  const DoFHandler<1, 1> &,
                                  const Quadrature<0> &,
                                  const Function<1> &,
                                  Vector<double> &,
                                  const std::set<types::boundary_id> &)
  {
    Assert(false, ExcImpossibleInDim(1));
  }



  template <>
  void
  create_boundary_right_hand_side(const Mapping<1, 2> &,
                                  const DoFHandler<1, 2> &,
                                  const Quadrature<0> &,
                                  const Function<2> &,
                                  Vector<double> &,
                                  const std::set<types::boundary_id> &)
  {
    Assert(false, ExcImpossibleInDim(1));
  }



} // namespace VectorTools

// ---------------------------- explicit instantiations --------------------
#include "numerics/vector_tools_rhs.inst"

DEAL_II_NAMESPACE_CLOSE
