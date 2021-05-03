// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_distributed_cell_data_transfer_templates_h
#define dealii_distributed_cell_data_transfer_templates_h


#include <deal.II/base/config.h>

#include <deal.II/distributed/cell_data_transfer.h>

#ifdef DEAL_II_WITH_P4EST

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/la_vector.h>
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/trilinos_epetra_vector.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace parallel
  {
    namespace distributed
    {
      namespace CellDataTransferImplementation
      {
        template <typename VectorType>
        void
        post_unpack_action(std::vector<VectorType *> &all_out)
        {
          for (auto &out : all_out)
            out->compress(::dealii::VectorOperation::insert);
        }

        template <typename value_type>
        void
        post_unpack_action(std::vector<std::vector<value_type> *> &)
        {
          // Do nothing for std::vector as VectorType.
        }
      } // namespace CellDataTransferImplementation
    }   // namespace distributed
  }     // namespace parallel
} // namespace internal



namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim, typename VectorType>
    CellDataTransfer<dim, spacedim, VectorType>::CellDataTransfer(
      const parallel::distributed::Triangulation<dim, spacedim> &triangulation,
      const bool                        transfer_variable_size_data,
      const std::function<std::vector<value_type>(
        const typename dealii::Triangulation<dim, spacedim>::cell_iterator
          &              parent,
        const value_type parent_value)> refinement_strategy,
      const std::function<value_type(
        const typename dealii::Triangulation<dim, spacedim>::cell_iterator
          &                            parent,
        const std::vector<value_type> &children_values)> coarsening_strategy)
      : triangulation(&triangulation, typeid(*this).name())
      , transfer_variable_size_data(transfer_variable_size_data)
      , refinement_strategy(refinement_strategy)
      , coarsening_strategy(coarsening_strategy)
      , handle(numbers::invalid_unsigned_int)
    {}

    template <int dim, int spacedim, typename VectorType>
    DEAL_II_DEPRECATED
    CellDataTransfer<dim, spacedim, VectorType>::CellDataTransfer(
      const parallel::distributed::Triangulation<dim, spacedim> &triangulation,
      const bool transfer_variable_size_data,
      const std::function<
        value_type(const typename parallel::distributed::
                     Triangulation<dim, spacedim>::cell_iterator &parent,
                   const VectorType &input_vector)> coarsening_strategy)
      : triangulation(&triangulation, typeid(*this).name())
      , transfer_variable_size_data(transfer_variable_size_data)
      , refinement_strategy(&dealii::AdaptationStrategies::Refinement::
                              preserve<dim, spacedim, value_type>)
      , handle(numbers::invalid_unsigned_int)
    {
      value_type (*const *old_strategy)(
        const typename parallel::distributed::Triangulation<dim, spacedim>::
          cell_iterator &,
        const VectorType &) =
        coarsening_strategy.template target<
          value_type (*)(const typename parallel::distributed::
                           Triangulation<dim, spacedim>::cell_iterator &,
                         const VectorType &)>();

      if (*old_strategy == CoarseningStrategies::check_equality)
        const_cast<std::function<value_type(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
          const std::vector<value_type> &)> &>(this->coarsening_strategy) =
          &dealii::AdaptationStrategies::Coarsening::
            check_equality<dim, spacedim, value_type>;
      else if (*old_strategy == CoarseningStrategies::sum)
        const_cast<std::function<value_type(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
          const std::vector<value_type> &)> &>(this->coarsening_strategy) =
          &dealii::AdaptationStrategies::Coarsening::
            sum<dim, spacedim, value_type>;
      else if (*old_strategy == CoarseningStrategies::mean)
        const_cast<std::function<value_type(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
          const std::vector<value_type> &)> &>(this->coarsening_strategy) =
          &dealii::AdaptationStrategies::Coarsening::
            mean<dim, spacedim, value_type>;
      else
        Assert(
          false,
          ExcMessage(
            "The constructor using the former function type of the "
            "coarsening_strategy parameter is no longer supported. Please use "
            "the latest function type instead"));
    }


    // Interface for packing
    // ---------------------

    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::register_data_attach()
    {
      // TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim, spacedim> *tria =
        const_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
          &(*triangulation));
      Assert(tria != nullptr, ExcInternalError());

      handle = tria->register_data_attach(
        [this](const typename parallel::distributed::
                 Triangulation<dim, spacedim>::cell_iterator &cell,
               const typename parallel::distributed::
                 Triangulation<dim, spacedim>::CellStatus status) {
          return this->pack_callback(cell, status);
        },
        /*returns_variable_size_data=*/transfer_variable_size_data);
    }



    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::
      prepare_for_coarsening_and_refinement(
        const std::vector<const VectorType *> &all_in)
    {
      input_vectors = all_in;
      register_data_attach();
    }



    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::
      prepare_for_coarsening_and_refinement(const VectorType &in)
    {
      std::vector<const VectorType *> all_in(1, &in);
      prepare_for_coarsening_and_refinement(all_in);
    }



    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::prepare_for_serialization(
      const std::vector<const VectorType *> &all_in)
    {
      prepare_for_coarsening_and_refinement(all_in);
    }



    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::prepare_for_serialization(
      const VectorType &in)
    {
      std::vector<const VectorType *> all_in(1, &in);
      prepare_for_serialization(all_in);
    }



    // Interface for unpacking
    // -----------------------

    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::unpack(
      std::vector<VectorType *> &all_out)
    {
      Assert(input_vectors.size() == all_out.size(),
             ExcDimensionMismatch(input_vectors.size(), all_out.size()));

      // TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim, spacedim> *tria =
        const_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
          &(*triangulation));
      Assert(tria != nullptr, ExcInternalError());

      tria->notify_ready_to_unpack(
        handle,
        [this, &all_out](
          const typename parallel::distributed::Triangulation<dim, spacedim>::
            cell_iterator &cell,
          const typename parallel::distributed::Triangulation<dim, spacedim>::
            CellStatus status,
          const boost::iterator_range<std::vector<char>::const_iterator>
            &data_range) {
          this->unpack_callback(cell, status, data_range, all_out);
        });

      dealii::internal::parallel::distributed::CellDataTransferImplementation::
        post_unpack_action(all_out);

      input_vectors.clear();
    }



    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::unpack(VectorType &out)
    {
      std::vector<VectorType *> all_out(1, &out);
      unpack(all_out);
    }



    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::deserialize(
      std::vector<VectorType *> &all_out)
    {
      // For deserialization, we need to register this object
      // to the triangulation first to get a valid handle for
      // data access.
      register_data_attach();

      // This makes unpack() happy.
      input_vectors.resize(all_out.size());

      unpack(all_out);
    }



    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::deserialize(VectorType &out)
    {
      std::vector<VectorType *> all_out(1, &out);
      deserialize(all_out);
    }



    // Callback functions
    // ------------------

    template <int dim, int spacedim, typename VectorType>
    std::vector<char>
    CellDataTransfer<dim, spacedim, VectorType>::pack_callback(
      const typename parallel::distributed::Triangulation<dim, spacedim>::
        cell_iterator &cell,
      const typename parallel::distributed::Triangulation<dim,
                                                          spacedim>::CellStatus
        status)
    {
      std::vector<value_type> cell_data(input_vectors.size());

      // Extract data from input_vectors for this particular cell.
      auto it_input  = input_vectors.cbegin();
      auto it_output = cell_data.begin();
      for (; it_input != input_vectors.cend(); ++it_input, ++it_output)
        {
          switch (status)
            {
              case parallel::distributed::Triangulation<dim,
                                                        spacedim>::CELL_PERSIST:
              case parallel::distributed::Triangulation<dim,
                                                        spacedim>::CELL_REFINE:
                // Cell either persists, or will be refined, and its children do
                // not exist yet in the latter case.
                *it_output = (**it_input)[cell->active_cell_index()];
                break;

              case parallel::distributed::Triangulation<dim,
                                                        spacedim>::CELL_COARSEN:
                {
                  // Cell is parent whose children will get coarsened to.
                  // Decide data to store on parent by provided strategy.
                  std::vector<value_type> children_values(cell->n_children());
                  for (unsigned int child_index = 0;
                       child_index < cell->n_children();
                       ++child_index)
                    {
                      const auto &child = cell->child(child_index);
                      Assert(child->is_active() && child->coarsen_flag_set(),
                             typename dealii::Triangulation<
                               dim>::ExcInconsistentCoarseningFlags());

                      children_values[child_index] =
                        (**it_input)[child->active_cell_index()];
                    }

                  *it_output = coarsening_strategy(cell, children_values);
                  break;
                }

              default:
                Assert(false, ExcInternalError());
                break;
            }
        }

      // We don't have to pack the whole container if there is just one entry.
      if (input_vectors.size() == 1)
        return Utilities::pack(
          cell_data[0], /*allow_compression=*/transfer_variable_size_data);
      else
        return Utilities::pack(
          cell_data, /*allow_compression=*/transfer_variable_size_data);
    }



    template <int dim, int spacedim, typename VectorType>
    void
    CellDataTransfer<dim, spacedim, VectorType>::unpack_callback(
      const typename parallel::distributed::Triangulation<dim, spacedim>::
        cell_iterator &cell,
      const typename parallel::distributed::Triangulation<dim,
                                                          spacedim>::CellStatus
        status,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &                        data_range,
      std::vector<VectorType *> &all_out)
    {
      std::vector<value_type> cell_data;

      // We have to unpack the corresponding datatype that has been packed
      // beforehand.
      if (all_out.size() == 1)
        cell_data.push_back(Utilities::unpack<value_type>(
          data_range.begin(),
          data_range.end(),
          /*allow_compression=*/transfer_variable_size_data));
      else
        cell_data = Utilities::unpack<std::vector<value_type>>(
          data_range.begin(),
          data_range.end(),
          /*allow_compression=*/transfer_variable_size_data);

      // Check if sizes match.
      Assert(cell_data.size() == all_out.size(), ExcInternalError());

      auto it_input  = cell_data.cbegin();
      auto it_output = all_out.begin();
      for (; it_input != cell_data.cend(); ++it_input, ++it_output)
        switch (status)
          {
            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_PERSIST:
            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_COARSEN:
              // Cell either persists, or has been coarsened.
              // Thus, cell has no (longer) children.
              (**it_output)[cell->active_cell_index()] = *it_input;
              break;

            case parallel::distributed::Triangulation<dim,
                                                      spacedim>::CELL_REFINE:
              {
                // Cell has been refined, and is now parent of its children.
                // Thus, distribute parent's data on its children.
                const std::vector<value_type> children_values =
                  refinement_strategy(cell, *it_input);
                for (unsigned int child_index = 0;
                     child_index < cell->n_children();
                     ++child_index)
                  (**it_output)[cell->child(child_index)->active_cell_index()] =
                    children_values[child_index];
                break;
              }

            default:
              Assert(false, ExcInternalError());
              break;
          }
    }
  } // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif /* DEAL_II_WITH_P4EST */

#endif /* dealii_distributed_cell_data_transfer_templates_h */
