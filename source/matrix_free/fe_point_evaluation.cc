// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_piecewise.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_q_dg0.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/shape_info.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace FEPointEvaluation
  {
    /**
     * Similar to internal::MatrixFreeFunctions::ShapeInfo but not supporting
     * non-tensor product elements and FE_DGQHermite
     */
    template <int dim, int spacedim>
    bool
    is_fast_path_supported(const FiniteElement<dim, spacedim> &fe,
                           const unsigned int base_element_number)
    {
      // check if supported
      const bool flag = [&]() {
        const FiniteElement<dim, spacedim> *fe_ptr =
          &(fe.base_element(base_element_number));
        if (fe_ptr->n_components() != 1)
          return false;

        // then check if the base element is supported or not
        if (dynamic_cast<const FE_Poly<dim, spacedim> *>(fe_ptr) != nullptr)
          {
            const FE_Poly<dim, spacedim> *fe_poly_ptr =
              dynamic_cast<const FE_Poly<dim, spacedim> *>(fe_ptr);

            if (dynamic_cast<const TensorProductPolynomials<dim> *>(
                  &fe_poly_ptr->get_poly_space()) == nullptr)
              return false;
          }
        else
          return false;

        return true;
      }();

      // make sure that if supported also ShapeInfo is supporting it
      if (flag)
        Assert(internal::MatrixFreeFunctions::ShapeInfo<double>::is_supported(
                 fe),
               ExcInternalError());

      return flag;
    }



    template <int dim, int spacedim>
    bool
    is_fast_path_supported(const Mapping<dim, spacedim> &mapping)
    {
      if (dynamic_cast<const MappingQ<dim, spacedim> *>(&mapping))
        {
          return true;
        }
      else if (dynamic_cast<const MappingCartesian<dim, spacedim> *>(&mapping))
        {
          return true;
        }
      return false;
    }



    template <int dim, int spacedim>
    std::vector<Polynomials::Polynomial<double>>
    get_polynomial_space(const FiniteElement<dim, spacedim> &fe)
    {
      Assert(fe.n_components() == 1, ExcNotImplemented());
      const FE_Poly<dim, spacedim> *fe_poly_ptr =
        dynamic_cast<const FE_Poly<dim, spacedim> *>(&fe);

      // we should catch the case that we cannot dynamic cast in
      // is_fast_path_supported
      Assert(fe_poly_ptr != nullptr, ExcNotImplemented());
      if (const auto polyspace =
            dynamic_cast<const TensorProductPolynomials<dim> *>(
              &fe_poly_ptr->get_poly_space()))
        return polyspace->get_underlying_polynomials();
      else
        DEAL_II_NOT_IMPLEMENTED();
      return {};
    }
  } // namespace FEPointEvaluation
} // namespace internal


// explicit instantiations
#include "matrix_free/fe_point_evaluation.inst"


DEAL_II_NAMESPACE_CLOSE
