// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_distributed_field_transfer_h
#define dealii_distributed_field_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/cell_data_transfer.templates.h>

#include <deal.II/dofs/dof_handler.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    /**
     * Classes and functions in the experimental namespace can be modified or
     * removed without being deprecated first.
     */
    namespace experimental
    {
      /**
       * This class is similar to SolutionTransfer but it supports the case
       * where elements have been activated during refinement, i.e., FE_Nothing
       * elements have been associated with a finite elements during refinement.
       */
      template <int dim, typename VectorType, int spacedim = dim>
      class FieldTransfer
      {
      private:
        using Number = typename VectorType::value_type;

      public:
        /**
         * Constructor.
         *
         * @param[in] dof_handler The DoFHandler on which all the operations
         * will happen. This constructor must be called before the underlying
         * Triangulation is coarsened/refined.
         */
        FieldTransfer(const DoFHandler<dim, spacedim> &dof_handler);

        /**
         * Prepare the current object for coarsening and refinement.
         * @param[in] in The vector that will be interpolated
         * @param[in] fe_nothing_index The finite element index associated with
         * FE_Nothing
         */
        void
        prepare_for_coarsening_and_refinement(
          const VectorType  &in,
          const unsigned int fe_nothing_index);

        /**
         * Interpolate the data previously stored in this object before the mesh
         * was refined or coarsened onto the current set of cells. @p new_value
         * is the value associated to the new degrees of freedom that where
         * created during the element activation. @p affine_constraints is the
         * AffineConstraints after refinement.
         */
        void
        interpolate(const Number                    &new_value,
                    const AffineConstraints<Number> &affine_constraints,
                    VectorType                      &out);

      private:
        /**
         * DoFHandler associated with the object.
         */
        const DoFHandler<dim, spacedim> &dof_handler;

        /**
         * Data transferred by cell_data_transfer.
         */
        std::vector<Vector<Number>> data_to_transfer;

        /**
         * CellDataTransfer used to perform the field transfer.
         */
        std::unique_ptr<
          CellDataTransfer<dim, spacedim, std::vector<Vector<Number>>>>
          cell_data_transfer;
      };
    } // namespace experimental
  }   // namespace distributed
} // namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif
