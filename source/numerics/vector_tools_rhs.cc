// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/numerics/vector_tools.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace VectorTools
{

// separate implementation for 1D because otherwise we get linker errors since
// (hp::)FEFaceValues<1> is not compiled
  template <>
  void
  create_boundary_right_hand_side (const Mapping<1,1> &,
                                   const DoFHandler<1,1> &,
                                   const Quadrature<0> &,
                                   const Function<1> &,
                                   Vector<double> &,
                                   const std::set<types::boundary_id> &)
  {
    Assert (false, ExcImpossibleInDim(1));
  }



  template <>
  void
  create_boundary_right_hand_side (const Mapping<1,2> &,
                                   const DoFHandler<1,2> &,
                                   const Quadrature<0> &,
                                   const Function<2> &,
                                   Vector<double> &,
                                   const std::set<types::boundary_id> &)
  {
    Assert (false, ExcImpossibleInDim(1));
  }



  template <>
  void
  create_boundary_right_hand_side (const hp::MappingCollection<1,1> &,
                                   const hp::DoFHandler<1,1> &,
                                   const hp::QCollection<0> &,
                                   const Function<1> &,
                                   Vector<double> &,
                                   const std::set<types::boundary_id> &)
  {
    Assert (false, ExcImpossibleInDim(1));
  }



  template <>
  void
  create_boundary_right_hand_side (const hp::MappingCollection<1,2> &,
                                   const hp::DoFHandler<1,2> &,
                                   const hp::QCollection<0> &,
                                   const Function<2> &,
                                   Vector<double> &,
                                   const std::set<types::boundary_id> &)
  {
    Assert (false, ExcImpossibleInDim(1));
  }
}

// ---------------------------- explicit instantiations --------------------
#include "vector_tools_rhs.inst"

DEAL_II_NAMESPACE_CLOSE
