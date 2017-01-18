// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


#ifndef dealii__mg_transfer_internal_h
#define dealii__mg_transfer_internal_h

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
     */
    template <int dim, int spacedim>
    void fill_copy_indices(const dealii::DoFHandler<dim,spacedim>                                                  &mg_dof,
                           const MGConstrainedDoFs                                                                 *mg_constrained_dofs,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_global_mine,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_level_mine);



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
     * Stores data related to the finite element contained in the
     * DoFHandler. Used only for the initialization using
     * <tt>setup_transfer</tt>.
     */
    template <typename Number>
    struct ElementInfo
    {

      /**
       * Stores the degree of the finite element. The selection of the
       * computational kernel is based on this number.
       */
      unsigned int fe_degree;

      /**
       * Stores whether the element is continuous and there is a joint degree of
       * freedom in the center of the 1D line.
       */
      bool element_is_continuous;

      /**
       * Stores the number of components in the finite element.
       */
      unsigned int n_components;

      /**
       * Stores the number of degrees of freedom on all child cells. It is
       * <tt>2<sup>dim</sup>*fe.dofs_per_cell</tt> for DG elements and somewhat
       * less for continuous elements.
       */
      unsigned int n_child_cell_dofs;

      /**
       * Holds the numbering between the numbering of degrees of freedom in
       * the finite element and the lexicographic numbering needed for the
       * tensor product application.
       */
      std::vector<unsigned int> lexicographic_numbering;

      /**
       * Holds the one-dimensional embedding (prolongation) matrix from mother
       * element to all the children.
       */
      AlignedVector<VectorizedArray<Number> > prolongation_matrix_1d;

    };

    /**
     * Sets up most of the internal data structures of MGTransferMatrixFree
     */
    template <int dim, typename Number>
    void setup_transfer(const dealii::DoFHandler<dim>                                       &mg_dof,
                        const MGConstrainedDoFs                                             *mg_constrained_dofs,
                        ElementInfo<Number>                                                 &elem_info,
                        std::vector<std::vector<unsigned int> >                             &level_dof_indices,
                        std::vector<std::vector<std::pair<unsigned int,unsigned int> > >    &parent_child_connect,
                        std::vector<unsigned int>                                           &n_owned_level_cells,
                        std::vector<std::vector<std::vector<unsigned short> > >             &dirichlet_indices,
                        std::vector<std::vector<Number> >                                   &weights_on_refined,
                        std::vector<std::vector<std::pair<unsigned int, unsigned int> > >   &copy_indices_global_mine,
                        MGLevelObject<LinearAlgebra::distributed::Vector<Number> >          &ghosted_level_vector);

  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
