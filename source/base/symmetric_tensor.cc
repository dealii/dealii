// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

// Required for instantiation of template functions
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/symmetric_tensor.templates.h>

#include <deal.II/differentiation/ad/adolc_product_types.h>
#include <deal.II/differentiation/ad/sacado_product_types.h>

#ifdef DEAL_II_WITH_ADOLC
#  include <adolc/adouble.h>
#  include <adolc/adtl.h>
#endif


DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_ADOLC
#  ifdef DEAL_II_ADOLC_WITH_ADVANCED_BRANCHING

// Specializations of the above functions for taped ADOL-C numbers
// when the advanced branching feature is activated.
// We could copy-paste all of these functions and add the appropriate
// conditional assignments (see the ADOL-C manual, section 1.8).
// However, some of the conditions are quite complicated (with possibly
// no 1-1 correspondence for the operations along each branch) so this
// needs some careful attention. For the sake of simplicity and until
// this can be rigourously checked, at this point in time we simply
// disable these functions.
namespace internal
{
  namespace SymmetricTensorImplementation
  {
    template <>
    struct Inverse<4, 3, adouble>
    {
      static dealii::SymmetricTensor<4, 3, adouble>
      value(const dealii::SymmetricTensor<4, 3, adouble> & /*t*/)
      {
        AssertThrow(false, ExcADOLCAdvancedBranching());
        return dealii::SymmetricTensor<4, 3, adouble>();
      }
    };
  } // namespace SymmetricTensorImplementation
} // namespace internal

template <>
std::array<adouble, 1>
eigenvalues(const SymmetricTensor<2, 1, adouble> & /*T*/)
{
  AssertThrow(false, ExcADOLCAdvancedBranching());
  return std::array<adouble, 1>();
}



template <>
std::array<adouble, 2>
eigenvalues(const SymmetricTensor<2, 2, adouble> & /*T*/)
{
  AssertThrow(false, ExcADOLCAdvancedBranching());
  return std::array<adouble, 2>();
}



template <>
std::array<adouble, 3>
eigenvalues(const SymmetricTensor<2, 3, adouble> & /*T*/)
{
  AssertThrow(false, ExcADOLCAdvancedBranching());
  return std::array<adouble, 3>();
}



template <>
std::array<std::pair<adouble, Tensor<1, 1, adouble>>, 1>
eigenvectors<1, adouble>(const SymmetricTensor<2, 1, adouble> & /*T*/,
                         const SymmetricTensorEigenvectorMethod /*method*/)
{
  AssertThrow(false, ExcADOLCAdvancedBranching());
  return std::array<std::pair<adouble, Tensor<1, 1, adouble>>, 1>();
}



template <>
std::array<std::pair<adouble, Tensor<1, 2, adouble>>, 2>
eigenvectors<2, adouble>(const SymmetricTensor<2, 2, adouble> & /*T*/,
                         const SymmetricTensorEigenvectorMethod /*method*/)
{
  AssertThrow(false, ExcADOLCAdvancedBranching());
  return std::array<std::pair<adouble, Tensor<1, 2, adouble>>, 2>();
}



template <>
std::array<std::pair<adouble, Tensor<1, 3, adouble>>, 3>
eigenvectors<3, adouble>(const SymmetricTensor<2, 3, adouble> & /*T*/,
                         const SymmetricTensorEigenvectorMethod /*method*/)
{
  AssertThrow(false, ExcADOLCAdvancedBranching());
  return std::array<std::pair<adouble, Tensor<1, 3, adouble>>, 3>();
}

#  else

template std::array<adouble, 1>
eigenvalues(const SymmetricTensor<2, 1, adouble> &);

template std::array<adouble, 2>
eigenvalues(const SymmetricTensor<2, 2, adouble> &);

template std::array<adouble, 3>
eigenvalues(const SymmetricTensor<2, 3, adouble> &);

template std::array<std::pair<adouble, Tensor<1, 1, adouble>>, 1>
eigenvectors<1, adouble>(const SymmetricTensor<2, 1, adouble> &,
                         const SymmetricTensorEigenvectorMethod);

template std::array<std::pair<adouble, Tensor<1, 2, adouble>>, 2>
eigenvectors<2, adouble>(const SymmetricTensor<2, 2, adouble> &,
                         const SymmetricTensorEigenvectorMethod);

template std::array<std::pair<adouble, Tensor<1, 3, adouble>>, 3>
eigenvectors<3, adouble>(const SymmetricTensor<2, 3, adouble> &,
                         const SymmetricTensorEigenvectorMethod);
#  endif

template std::array<adtl::adouble, 1>
eigenvalues(const SymmetricTensor<2, 1, adtl::adouble> &);

template std::array<adtl::adouble, 2>
eigenvalues(const SymmetricTensor<2, 2, adtl::adouble> &);

template std::array<adtl::adouble, 3>
eigenvalues(const SymmetricTensor<2, 3, adtl::adouble> &);

template std::array<std::pair<adtl::adouble, Tensor<1, 1, adtl::adouble>>, 1>
eigenvectors<1, adtl::adouble>(const SymmetricTensor<2, 1, adtl::adouble> &,
                               const SymmetricTensorEigenvectorMethod);

template std::array<std::pair<adtl::adouble, Tensor<1, 2, adtl::adouble>>, 2>
eigenvectors<2, adtl::adouble>(const SymmetricTensor<2, 2, adtl::adouble> &,
                               const SymmetricTensorEigenvectorMethod);

template std::array<std::pair<adtl::adouble, Tensor<1, 3, adtl::adouble>>, 3>
eigenvectors<3, adtl::adouble>(const SymmetricTensor<2, 3, adtl::adouble> &,
                               const SymmetricTensorEigenvectorMethod);
#endif

// explicit instantiations
#include "base/symmetric_tensor.inst"


DEAL_II_NAMESPACE_CLOSE
