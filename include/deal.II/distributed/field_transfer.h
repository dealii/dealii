// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2023 by the deal.II authors
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
