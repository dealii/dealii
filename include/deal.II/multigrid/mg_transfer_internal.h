// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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


#ifndef dealii_mg_transfer_internal_h
#define dealii_mg_transfer_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MGTransfer
  {
    /**
     * Internal function for filling the copy indices from global to level
     * indices
     *
     * If @p skip_interface_dofs is false, the mapping will also contain
     * DoFs at the interface between levels. This is desirable when
     * transferring solution vectors instead of residuals.
     */
    template <int dim, int spacedim>
    void
    fill_copy_indices(
      const DoFHandler<dim, spacedim> &dof_handler,
      const MGConstrainedDoFs *        mg_constrained_dofs,
      std::vector<std::vector<
        std::pair<types::global_dof_index, types::global_dof_index>>>
        &copy_indices,
      std::vector<std::vector<
        std::pair<types::global_dof_index, types::global_dof_index>>>
        &copy_indices_global_mine,
      std::vector<std::vector<
        std::pair<types::global_dof_index, types::global_dof_index>>>
        &        copy_indices_level_mine,
      const bool skip_interface_dofs = true);



    /**
     * Given the collection of child cells in lexicographic ordering as seen
     * from the parent, this function computes the first index of the given
     * child
     */
    template <int dim>
    unsigned int
    compute_shift_within_children(const unsigned int child,
                                  const unsigned int fe_shift_1d,
                                  const unsigned int fe_degree);

    /**
     * A structure that stores data related to the finite element contained in
     * the DoFHandler. Used only for the initialization using
     * <tt>setup_transfer</tt>.
     */
    template <typename Number>
    struct ElementInfo
    {
      /**
       * A variable storing the degree of the finite element. The selection of
       * the computational kernel is based on this number.
       */
      unsigned int fe_degree;

      /**
       * A variable storing whether the element is continuous and there is a
       * joint degree of freedom in the center of the 1D line.
       */
      bool element_is_continuous;

      /**
       * A variable storing the number of components in the finite element.
       */
      unsigned int n_components;

      /**
       * A variable storing the number of degrees of freedom on all child cells.
       * It is <tt>2<sup>dim</sup>*fe.dofs_per_cell</tt> for DG elements and
       * somewhat less for continuous elements.
       */
      unsigned int n_child_cell_dofs;

      /**
       * An array that holds the numbering between the numbering of degrees of
       * freedom in the finite element and the lexicographic numbering needed
       * for the tensor product application.
       */
      std::vector<unsigned int> lexicographic_numbering;

      /**
       * This variable holds the one-dimensional embedding (prolongation) matrix
       * from mother element to all the children.
       */
      std::vector<Number> prolongation_matrix_1d;
    };

    /**
     * Set up most of the internal data structures of MGTransferMatrixFree
     */
    template <int dim, typename Number>
    void
    setup_transfer(
      const DoFHandler<dim> &  dof_handler,
      const MGConstrainedDoFs *mg_constrained_dofs,
      const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        &                                     external_partitioners,
      ElementInfo<Number> &                   elem_info,
      std::vector<std::vector<unsigned int>> &level_dof_indices,
      std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
        &                        parent_child_connect,
      std::vector<unsigned int> &n_owned_level_cells,
      std::vector<std::vector<std::vector<unsigned short>>> &dirichlet_indices,
      std::vector<std::vector<Number>> &                     weights_on_refined,
      std::vector<Table<2, unsigned int>> &copy_indices_global_mine,
      MGLevelObject<std::shared_ptr<const Utilities::MPI::Partitioner>>
        &vector_partitioners);

  } // namespace MGTransfer
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
