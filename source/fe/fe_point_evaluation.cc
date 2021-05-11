// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_piecewise.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_point_evaluation.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_q_dg0.h>

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
    is_fast_path_supported(const FiniteElement<dim, spacedim> &fe)
    {
      // check if supported
      const bool flag = [&]() {
        if (dim != spacedim)
          return false;

        for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
          {
            const FiniteElement<dim, spacedim> *fe_ptr =
              &(fe.base_element(base));
            if (fe_ptr->n_components() != 1)
              return false;

            // then check if the base element is supported or not
            if (dynamic_cast<const FE_Poly<dim, spacedim> *>(fe_ptr) != nullptr)
              {
                const FE_Poly<dim, spacedim> *fe_poly_ptr =
                  dynamic_cast<const FE_Poly<dim, spacedim> *>(fe_ptr);

                if (dynamic_cast<const TensorProductPolynomials<dim> *>(
                      &fe_poly_ptr->get_poly_space()) == nullptr &&
                    dynamic_cast<const TensorProductPolynomials<
                        dim,
                        Polynomials::PiecewisePolynomial<double>> *>(
                      &fe_poly_ptr->get_poly_space()) == nullptr &&
                    dynamic_cast<const FE_DGP<dim, spacedim> *>(fe_ptr) ==
                      nullptr &&
                    dynamic_cast<const FE_Q_DG0<dim, spacedim> *>(fe_ptr) ==
                      nullptr &&
                    dynamic_cast<const FE_DGQHermite<dim, spacedim> *>(
                      fe_ptr) == nullptr)
                  return false;
              }
            else
              return false;
          }

        // if we arrived here, all base elements were supported so we can
        // support the present element
        return true;
      }();

      // make sure that if supported also ShapeInfo is supporting it
      if (flag)
        Assert(internal::MatrixFreeFunctions::ShapeInfo<double>::is_supported(
                 fe),
               ExcInternalError());

      return flag;
    }
  } // namespace FEPointEvaluation
} // namespace internal


// explicit instantiations
#include "fe_point_evaluation.inst"


DEAL_II_NAMESPACE_CLOSE
