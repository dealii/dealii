// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_mapping_related_data_h
#define dealii_fe_mapping_related_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_update_flags.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

/** @addtogroup feaccess */
/** @{ */

namespace internal
{
  namespace FEValuesImplementation
  {
    /**
     * A class that stores all of the mapping related data used in
     * dealii::FEValues, dealii::FEFaceValues, and dealii::FESubfaceValues
     * objects. Objects of this kind will be given as <i>output</i> argument
     * when dealii::FEValues::reinit() calls Mapping::fill_fe_values() for a
     * given cell, face, or subface.
     *
     * The data herein will then be provided as <i>input</i> argument in the
     * following call to FiniteElement::fill_fe_values().
     *
     * @ingroup feaccess
     */
    template <int dim, int spacedim = dim>
    class MappingRelatedData
    {
    public:
      /**
       * Initialize all vectors to correct size.
       */
      void
      initialize(const unsigned int n_quadrature_points,
                 const UpdateFlags  flags);

      /**
       * Compute and return an estimate for the memory consumption (in bytes)
       * of this object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Store an array of weights times the Jacobi determinant at the
       * quadrature points. This function is reset each time reinit() is
       * called. The Jacobi determinant is actually the reciprocal value of
       * the Jacobi matrices stored in this class, see the general
       * documentation of this class for more information.
       *
       * However, if this object refers to an FEFaceValues or FESubfaceValues
       * object, then the JxW_values correspond to the Jacobian of the
       * transformation of the face, not the cell, i.e. the dimensionality is
       * that of a surface measure, not of a volume measure. In this case, it
       * is computed from the boundary forms, rather than the Jacobian matrix.
       */
      std::vector<double> JxW_values;

      /**
       * Array of the Jacobian matrices at the quadrature points.
       */
      std::vector<DerivativeForm<1, dim, spacedim>> jacobians;

      /**
       * Array of the derivatives of the Jacobian matrices at the quadrature
       * points.
       */
      std::vector<DerivativeForm<2, dim, spacedim>> jacobian_grads;

      /**
       * Array of the inverse Jacobian matrices at the quadrature points.
       */
      std::vector<DerivativeForm<1, spacedim, dim>> inverse_jacobians;

      /**
       * Array of the derivatives of the Jacobian matrices at the quadrature
       * points, pushed forward to the real cell coordinates.
       */
      std::vector<Tensor<3, spacedim>> jacobian_pushed_forward_grads;

      /**
       * Array of the second derivatives of the Jacobian matrices at the
       * quadrature points.
       */
      std::vector<DerivativeForm<3, dim, spacedim>> jacobian_2nd_derivatives;

      /**
       * Array of the  second derivatives of the Jacobian matrices at the
       * quadrature points, pushed forward to the real cell coordinates.
       */
      std::vector<Tensor<4, spacedim>> jacobian_pushed_forward_2nd_derivatives;

      /**
       * Array of the  third derivatives of the Jacobian matrices at the
       * quadrature points.
       */
      std::vector<DerivativeForm<4, dim, spacedim>> jacobian_3rd_derivatives;

      /**
       * Array of the  third derivatives of the Jacobian matrices at the
       * quadrature points, pushed forward to the real cell coordinates.
       */
      std::vector<Tensor<5, spacedim>> jacobian_pushed_forward_3rd_derivatives;

      /**
       * Array of quadrature points. This array is set up upon calling
       * reinit() and contains the quadrature points on the real element,
       * rather than on the reference element.
       */
      std::vector<Point<spacedim>> quadrature_points;

      /**
       * List of outward normal vectors at the quadrature points.
       */
      std::vector<Tensor<1, spacedim>> normal_vectors;

      /**
       * List of boundary forms at the quadrature points.
       */
      std::vector<Tensor<1, spacedim>> boundary_forms;
    };
  } // namespace FEValuesImplementation
} // namespace internal


/** @} */



DEAL_II_NAMESPACE_CLOSE

#endif
