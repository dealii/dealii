// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2022 by the deal.II authors
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

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_levels.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_iterator.templates.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <limits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

template <typename Number>
DeclException2(ExcNonMatchingElementsSetDofValuesByInterpolation,
               Number,
               Number,
               << "Called set_dof_values_by_interpolation(), but"
               << " the element to be set, value " << std::setprecision(16)
               << arg1 << ", does not match with the non-zero value "
               << std::setprecision(16) << arg2 << " already set before.");

namespace internal
{
#ifdef DEBUG
  /**
   * In the set_dof_values(), we need to invoke abs() also on unsigned data
   * types, which is ill-formed on newer C++ standards. To avoid this, we use
   * std::abs on default types, but simply return the number on unsigned types.
   */
  template <typename Number>
  std::enable_if_t<!std::is_unsigned_v<Number>,
                   typename numbers::NumberTraits<Number>::real_type>
  get_abs(const Number a)
  {
    return std::abs(a);
  }

  template <typename Number>
  std::enable_if_t<std::is_unsigned_v<Number>, Number>
  get_abs(const Number a)
  {
    return a;
  }

  /**
   * Check if a vector is a deal.II vector.
   */
  template <typename VectorType>
  constexpr bool is_dealii_vector =
    std::is_same_v<VectorType,
                   dealii::Vector<typename VectorType::value_type>> ||
    std::is_same_v<VectorType,
                   dealii::BlockVector<typename VectorType::value_type>> ||
    std::is_same_v<VectorType,
                   dealii::LinearAlgebra::distributed::Vector<
                     typename VectorType::value_type>> ||
    std::is_same_v<VectorType,
                   dealii::LinearAlgebra::distributed::BlockVector<
                     typename VectorType::value_type>>;

  /**
   * Helper functions that call set_ghost_state() if the vector supports this
   * operation.
   */
  template <typename T>
  using set_ghost_state_t =
    decltype(std::declval<const T>().set_ghost_state(std::declval<bool>()));

  template <typename T>
  constexpr bool has_set_ghost_state =
    is_supported_operation<set_ghost_state_t, T>;

  template <
    typename VectorType,
    std::enable_if_t<has_set_ghost_state<VectorType>, VectorType> * = nullptr>
  void
  set_ghost_state(VectorType &vector, const bool ghosted)
  {
    vector.set_ghost_state(ghosted);
  }

  template <
    typename VectorType,
    std::enable_if_t<!has_set_ghost_state<VectorType>, VectorType> * = nullptr>
  void
  set_ghost_state(VectorType &, const bool)
  {
    // serial vector: nothing to do
  }
#endif

  /**
   * Helper function that sets the values on a cell, but also checks if the
   * new values are similar to the old values.
   */
  template <int  dim,
            int  spacedim,
            bool lda,
            class OutputVector,
            typename number>
  void
  set_dof_values(const DoFCellAccessor<dim, spacedim, lda> &cell,
                 const Vector<number>                      &local_values,
                 OutputVector                              &values,
                 const bool                                 perform_check)
  {
    (void)perform_check;

#ifdef DEBUG
    if (perform_check && is_dealii_vector<OutputVector>)
      {
        const bool old_ghost_state = values.has_ghost_elements();
        set_ghost_state(values, true);

        Vector<number> local_values_old(cell.get_fe().n_dofs_per_cell());
        cell.get_dof_values(values, local_values_old);

        for (unsigned int i = 0; i < cell.get_fe().n_dofs_per_cell(); ++i)
          {
            // a check consistent with the one in
            // Utilities::MPI::Partitioner::import_from_ghosted_array_finish()
            Assert(local_values_old[i] == number() ||
                     get_abs(local_values_old[i] - local_values[i]) <=
                       get_abs(local_values_old[i] + local_values[i]) *
                         100000. *
                         std::numeric_limits<typename numbers::NumberTraits<
                           number>::real_type>::epsilon(),
                   ExcNonMatchingElementsSetDofValuesByInterpolation<number>(
                     local_values[i], local_values_old[i]));
          }

        set_ghost_state(values, old_ghost_state);
      }
#endif

    cell.set_dof_values(local_values, values);
  }


  template <int  dim,
            int  spacedim,
            bool lda,
            class OutputVector,
            typename number>
  void
  process_by_interpolation(
    const DoFCellAccessor<dim, spacedim, lda>       &cell,
    const Vector<number>                            &local_values,
    OutputVector                                    &values,
    const types::fe_index                            fe_index_,
    const std::function<void(const DoFCellAccessor<dim, spacedim, lda> &cell,
                             const Vector<number> &local_values,
                             OutputVector         &values)> &processor)
  {
    const types::fe_index fe_index =
      (cell.get_dof_handler().has_hp_capabilities() == false &&
       fe_index_ == numbers::invalid_fe_index) ?
        DoFHandler<dim, spacedim>::default_fe_index :
        fe_index_;

    if (cell.is_active() && !cell.is_artificial())
      {
        if ((cell.get_dof_handler().has_hp_capabilities() == false) ||
            // for hp-DoFHandlers, we need to require that on
            // active cells, you either don't specify an fe_index,
            // or that you specify the correct one
            (fe_index == cell.active_fe_index()) ||
            (fe_index == numbers::invalid_fe_index))
          // simply set the values on this cell
          processor(cell, local_values, values);
        else
          {
            Assert(local_values.size() ==
                     cell.get_dof_handler().get_fe(fe_index).n_dofs_per_cell(),
                   ExcMessage("Incorrect size of local_values vector."));

            FullMatrix<double> interpolation(
              cell.get_fe().n_dofs_per_cell(),
              cell.get_dof_handler().get_fe(fe_index).n_dofs_per_cell());

            cell.get_fe().get_interpolation_matrix(
              cell.get_dof_handler().get_fe(fe_index), interpolation);

            // do the interpolation to the target space. for historical
            // reasons, matrices are set to size 0x0 internally even if
            // we reinit as 4x0, so we have to treat this case specially
            Vector<number> tmp(cell.get_fe().n_dofs_per_cell());
            if ((tmp.size() > 0) && (local_values.size() > 0))
              interpolation.vmult(tmp, local_values);

            // now set the dof values in the global vector
            processor(cell, tmp, values);
          }
      }
    else
      // otherwise distribute them to the children
      {
        Assert((cell.get_dof_handler().has_hp_capabilities() == false) ||
                 (fe_index != numbers::invalid_fe_index),
               ExcMessage(
                 "You cannot call this function on non-active cells "
                 "of DoFHandler objects unless you provide an explicit "
                 "finite element index because they do not have naturally "
                 "associated finite element spaces associated: degrees "
                 "of freedom are only distributed on active cells for which "
                 "the active FE index has been set."));

        const FiniteElement<dim, spacedim> &fe =
          cell.get_dof_handler().get_fe(fe_index);
        const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

        Assert(local_values.size() == dofs_per_cell,
               (typename DoFCellAccessor<dim, spacedim, lda>::BaseClass::
                  ExcVectorDoesNotMatch()));
        Assert(values.size() == cell.get_dof_handler().n_dofs(),
               (typename DoFCellAccessor<dim, spacedim, lda>::BaseClass::
                  ExcVectorDoesNotMatch()));

        Vector<number> tmp(dofs_per_cell);

        for (unsigned int child = 0; child < cell.n_children(); ++child)
          {
            if (tmp.size() > 0)
              fe.get_prolongation_matrix(child, cell.refinement_case())
                .vmult(tmp, local_values);
            process_by_interpolation(
              *cell.child(child), tmp, values, fe_index, processor);
          }
      }
  }

} // namespace internal



template <int dim, int spacedim, bool lda>
template <class OutputVector, typename number>
void
DoFCellAccessor<dim, spacedim, lda>::set_dof_values_by_interpolation(
  const Vector<number> &local_values,
  OutputVector         &values,
  const types::fe_index fe_index_,
  const bool            perform_check) const
{
  internal::process_by_interpolation<dim, spacedim, lda, OutputVector, number>(
    *this,
    local_values,
    values,
    fe_index_,
    [perform_check](const DoFCellAccessor<dim, spacedim, lda> &cell,
                    const Vector<number>                      &local_values,
                    OutputVector                              &values) {
      internal::set_dof_values(cell, local_values, values, perform_check);
    });
}


template <int dim, int spacedim, bool lda>
template <class OutputVector, typename number>
void
DoFCellAccessor<dim, spacedim, lda>::
  distribute_local_to_global_by_interpolation(
    const Vector<number> &local_values,
    OutputVector         &values,
    const types::fe_index fe_index_) const
{
  internal::process_by_interpolation<dim, spacedim, lda, OutputVector, number>(
    *this,
    local_values,
    values,
    fe_index_,
    [](const DoFCellAccessor<dim, spacedim, lda> &cell,
       const Vector<number>                      &local_values,
       OutputVector                              &values) {
      std::vector<types::global_dof_index> dof_indices(
        cell.get_fe().n_dofs_per_cell());
      cell.get_dof_indices(dof_indices);
      AffineConstraints<number>().distribute_local_to_global(local_values,
                                                             dof_indices,
                                                             values);
    });
}


// --------------------------------------------------------------------------
// explicit instantiations
#include "dof_accessor_set.inst"

DEAL_II_NAMESPACE_CLOSE
