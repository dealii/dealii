// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values_extractors.h>

DEAL_II_NAMESPACE_OPEN

namespace FEValuesExtractors
{
  std::string
  Scalar::get_name() const
  {
    return "Scalar(" + Utilities::int_to_string(component) + ")";
  }


  std::string
  Vector::get_name() const
  {
    return "Vector(" + Utilities::int_to_string(first_vector_component) + ")";
  }


  template <int rank>
  std::string
  Tensor<rank>::get_name() const
  {
    return "Tensor<" + Utilities::int_to_string(rank) + ">(" +
           Utilities::int_to_string(first_tensor_component) + ")";
  }


  template <int rank>
  std::string
  SymmetricTensor<rank>::get_name() const
  {
    return "SymmetricTensor<" + Utilities::int_to_string(rank) + ">(" +
           Utilities::int_to_string(first_tensor_component) + ")";
  }

  // Explicit instantiations
  template struct Tensor<0>;
  template struct Tensor<1>;
  template struct Tensor<2>;
  template struct Tensor<3>;
  template struct Tensor<4>;
  template struct SymmetricTensor<2>;
  template struct SymmetricTensor<4>;

} // namespace FEValuesExtractors



DEAL_II_NAMESPACE_CLOSE
