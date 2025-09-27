// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solution_transfer_templates_h
#define dealii_solution_transfer_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_element_access.h>

#include <deal.II/numerics/solution_transfer.h>

#include <boost/range/iterator_range_core.hpp>

#include <functional>
#include <numeric>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace SolutionTransferImplementation
  {
    /**
     * Optimized pack function for values assigned on degrees of freedom.
     *
     * Given that the elements of @p dof_values are stored in consecutive
     * locations, we can just memcpy them. Since floating point values don't
     * compress well, we also waive the compression that the default
     * Utilities::pack() and Utilities::unpack() functions offer.
     */
    template <typename value_type>
    std::vector<char>
    pack_dof_values(std::vector<Vector<value_type>> &dof_values,
                    const unsigned int               dofs_per_cell)
    {
      for (const auto &values : dof_values)
        {
          AssertDimension(values.size(), dofs_per_cell);
          (void)values;
        }

      const std::size_t bytes_per_entry = sizeof(value_type) * dofs_per_cell;

      std::vector<char> buffer(dof_values.size() * bytes_per_entry);
      for (unsigned int i = 0; i < dof_values.size(); ++i)
        std::memcpy(&buffer[i * bytes_per_entry],
                    &dof_values[i](0),
                    bytes_per_entry);

      return buffer;
    }



    /**
     * Optimized unpack function for values assigned on degrees of freedom.
     */
    template <typename value_type>
    std::vector<Vector<value_type>>
    unpack_dof_values(
      const boost::iterator_range<std::vector<char>::const_iterator>
                        &data_range,
      const unsigned int dofs_per_cell)
    {
      const std::size_t  bytes_per_entry = sizeof(value_type) * dofs_per_cell;
      const unsigned int n_elements      = data_range.size() / bytes_per_entry;

      Assert((data_range.size() % bytes_per_entry == 0), ExcInternalError());

      std::vector<Vector<value_type>> unpacked_data;
      unpacked_data.reserve(n_elements);
      for (unsigned int i = 0; i < n_elements; ++i)
        {
          Vector<value_type> dof_values(dofs_per_cell);
          std::memcpy(&dof_values(0),
                      &(*std::next(data_range.begin(), i * bytes_per_entry)),
                      bytes_per_entry);
          unpacked_data.emplace_back(std::move(dof_values));
        }

      return unpacked_data;
    }
  } // namespace SolutionTransferImplementation
} // namespace internal



template <int dim, typename VectorType, int spacedim>
SolutionTransfer<dim, VectorType, spacedim>::SolutionTransfer(
  const DoFHandler<dim, spacedim> &dof,
  const bool                       average_values)
  : dof_handler(&dof, typeid(*this).name())
  , average_values(average_values)
  , handle(numbers::invalid_unsigned_int)
{}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::
  prepare_for_coarsening_and_refinement(
    const std::vector<const VectorType *> &all_in)
{
  const dealii::internal::parallel::shared::
    TemporarilyRestoreSubdomainIds<dim, spacedim>
      subdomain_modifier(dof_handler->get_triangulation());

  for (unsigned int i = 0; i < all_in.size(); ++i)
    Assert(all_in[i]->size() == dof_handler->n_dofs(),
           ExcDimensionMismatch(all_in[i]->size(), dof_handler->n_dofs()));

  input_vectors = all_in;
  register_data_attach();
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::
  prepare_for_coarsening_and_refinement(const std::vector<VectorType> &all_in)
{
  std::vector<const VectorType *> temp(all_in.size());

  for (std::size_t i = 0; i < temp.size(); ++i)
    temp[i] = &(all_in[i]);

  this->prepare_for_coarsening_and_refinement(temp);
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::register_data_attach()
{
  // TODO: casting away constness is bad
  auto tria = const_cast<dealii::Triangulation<dim, spacedim> *>(
    &dof_handler->get_triangulation());
  Assert(tria != nullptr, ExcInternalError());

  Assert(handle == numbers::invalid_unsigned_int,
         ExcMessage("You can only add one solution per "
                    "SolutionTransfer object."));

  handle = tria->register_data_attach(
    [this](const typename Triangulation<dim, spacedim>::cell_iterator &cell_,
           const CellStatus                                            status) {
      return this->pack_callback(cell_, status);
    },
    /*returns_variable_size_data=*/dof_handler->has_hp_capabilities());
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::
  prepare_for_coarsening_and_refinement(const VectorType &in)
{
  std::vector<const VectorType *> all_in(1, &in);
  prepare_for_coarsening_and_refinement(all_in);
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::prepare_for_serialization(
  const VectorType &in)
{
  std::vector<const VectorType *> all_in(1, &in);
  prepare_for_serialization(all_in);
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::prepare_for_serialization(
  const std::vector<const VectorType *> &all_in)
{
  prepare_for_coarsening_and_refinement(all_in);
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::deserialize(VectorType &in)
{
  std::vector<VectorType *> all_in(1, &in);
  deserialize(all_in);
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::deserialize(
  std::vector<VectorType *> &all_in)
{
  register_data_attach();

  // this makes interpolate() happy
  input_vectors.resize(all_in.size());

  interpolate(all_in);
}


template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::interpolate(
  std::vector<VectorType *> &all_out)
{
  const dealii::internal::parallel::shared::
    TemporarilyRestoreSubdomainIds<dim, spacedim>
      subdomain_modifier(dof_handler->get_triangulation());

  Assert(input_vectors.size() == all_out.size(),
         ExcDimensionMismatch(input_vectors.size(), all_out.size()));
  for (unsigned int i = 0; i < all_out.size(); ++i)
    Assert(all_out[i]->size() == dof_handler->n_dofs(),
           ExcDimensionMismatch(all_out[i]->size(), dof_handler->n_dofs()));

  // TODO: casting away constness is bad
  auto tria = const_cast<dealii::Triangulation<dim, spacedim> *>(
    &dof_handler->get_triangulation());
  Assert(tria != nullptr, ExcInternalError());
  Assert(
    handle != numbers::invalid_unsigned_int,
    ExcMessage(
      "You can only call interpolate() once per SolutionTransfer object."));

  using Number = typename VectorType::value_type;

  if (average_values)
    for (auto *const vec : all_out)
      *vec = Number();

  VectorType valence;

  // initialize valence vector only if we need to average
  if (average_values)
    valence.reinit(*all_out[0]);

  tria->notify_ready_to_unpack(
    handle,
    [this, &all_out, &valence](
      const typename Triangulation<dim, spacedim>::cell_iterator &cell_,
      const CellStatus                                            status,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &data_range) {
      this->unpack_callback(cell_, status, data_range, all_out, valence);
    });

  if (average_values)
    {
      // finalize valence: compress and invert
      valence.compress(VectorOperation::add);
      for (const auto i : valence.locally_owned_elements())
        valence[i] = (static_cast<Number>(valence[i]) == Number() ?
                        Number() :
                        (Number(1.0) / static_cast<Number>(valence[i])));
      valence.compress(VectorOperation::insert);

      for (auto *const vec : all_out)
        {
          // compress and weight with valence
          vec->compress(VectorOperation::add);
          vec->scale(valence);
        }
    }
  else
    {
      for (auto *const vec : all_out)
        vec->compress(VectorOperation::insert);
    }

  input_vectors.clear();
  handle = numbers::invalid_unsigned_int;
}


template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::interpolate(
  std::vector<VectorType> &all_out)
{
  std::vector<VectorType *> temp(all_out.size());

  for (std::size_t i = 0; i < temp.size(); ++i)
    temp[i] = &(all_out[i]);

  this->interpolate(temp);
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::interpolate(VectorType &out)
{
  std::vector<VectorType *> all_out(1, &out);
  interpolate(all_out);
}



template <int dim, typename VectorType, int spacedim>
std::vector<char>
SolutionTransfer<dim, VectorType, spacedim>::pack_callback(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell_,
  const CellStatus                                            status)
{
  typename DoFHandler<dim, spacedim>::cell_iterator cell(*cell_, dof_handler);

  // create buffer for each individual object
  std::vector<::dealii::Vector<typename VectorType::value_type>> dof_values(
    input_vectors.size());

  unsigned int fe_index = 0;
  if (dof_handler->has_hp_capabilities())
    {
      switch (status)
        {
          case CellStatus::cell_will_persist:
          case CellStatus::cell_will_be_refined:
            {
              fe_index = cell->future_fe_index();
              break;
            }

          case CellStatus::children_will_be_coarsened:
            {
              // In case of coarsening, we need to find a suitable FE index
              // for the parent cell. We choose the 'least dominant fe'
              // on all children from the associated FECollection.
              if constexpr (running_in_debug_mode())
                {
                  for (const auto &child : cell->child_iterators())
                    Assert(child->is_active() && child->coarsen_flag_set(),
                           typename dealii::Triangulation<
                             dim>::ExcInconsistentCoarseningFlags());
                }

              fe_index = dealii::internal::hp::DoFHandlerImplementation::
                dominated_future_fe_on_children<dim, spacedim>(cell);
              break;
            }

          default:
            Assert(false, ExcInternalError());
            break;
        }
    }

  const unsigned int dofs_per_cell =
    dof_handler->get_fe(fe_index).n_dofs_per_cell();

  if (dofs_per_cell == 0)
    return std::vector<char>(); // nothing to do for FE_Nothing

  auto it_input  = input_vectors.cbegin();
  auto it_output = dof_values.begin();
  for (; it_input != input_vectors.cend(); ++it_input, ++it_output)
    {
      it_output->reinit(dofs_per_cell);
      cell->get_interpolated_dof_values(*(*it_input), *it_output, fe_index);
    }

  return internal::SolutionTransferImplementation::pack_dof_values<
    typename VectorType::value_type>(dof_values, dofs_per_cell);
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::unpack_callback(
  const typename Triangulation<dim, spacedim>::cell_iterator     &cell_,
  const CellStatus                                                status,
  const boost::iterator_range<std::vector<char>::const_iterator> &data_range,
  std::vector<VectorType *>                                      &all_out,
  VectorType                                                     &valence)
{
  typename DoFHandler<dim, spacedim>::cell_iterator cell(*cell_, dof_handler);

  unsigned int fe_index = 0;
  if (dof_handler->has_hp_capabilities())
    {
      switch (status)
        {
          case CellStatus::cell_will_persist:
          case CellStatus::children_will_be_coarsened:
            {
              fe_index = cell->active_fe_index();
              break;
            }

          case CellStatus::cell_will_be_refined:
            {
              // After refinement, this particular cell is no longer active,
              // and its children have inherited its FE index. However, to
              // unpack the data on the old cell, we need to recover its FE
              // index from one of the children. Just to be sure, we also
              // check if all children have the same FE index.
              fe_index = cell->child(0)->active_fe_index();
              for (unsigned int child_index = 1;
                   child_index < cell->n_children();
                   ++child_index)
                Assert(cell->child(child_index)->active_fe_index() == fe_index,
                       ExcInternalError());
              break;
            }

          default:
            Assert(false, ExcInternalError());
            break;
        }
    }

  const unsigned int dofs_per_cell =
    dof_handler->get_fe(fe_index).n_dofs_per_cell();

  if (dofs_per_cell == 0)
    return; // nothing to do for FE_Nothing

  const std::vector<::dealii::Vector<typename VectorType::value_type>>
    dof_values = internal::SolutionTransferImplementation::unpack_dof_values<
      typename VectorType::value_type>(data_range, dofs_per_cell);

  // check if sizes match
  AssertDimension(dof_values.size(), all_out.size());

  // check if we have enough dofs provided by the FE object
  // to interpolate the transferred data correctly
  for (auto it_dof_values = dof_values.begin();
       it_dof_values != dof_values.end();
       ++it_dof_values)
    Assert(
      dofs_per_cell == it_dof_values->size(),
      ExcMessage(
        "The transferred data was packed with a different number of dofs than the "
        "currently registered FE object assigned to the DoFHandler has."));

  // distribute data for each registered vector on mesh
  auto it_input  = dof_values.cbegin();
  auto it_output = all_out.begin();
  for (; it_input != dof_values.cend(); ++it_input, ++it_output)
    if (average_values)
      cell->distribute_local_to_global_by_interpolation(*it_input,
                                                        *(*it_output),
                                                        fe_index);
    else
      cell->set_dof_values_by_interpolation(*it_input,
                                            *(*it_output),
                                            fe_index,
                                            true);

  if (average_values)
    {
      // compute valence vector if averaging should be performed
      Vector<typename VectorType::value_type> ones(dofs_per_cell);
      ones = 1.0;
      cell->distribute_local_to_global_by_interpolation(ones,
                                                        valence,
                                                        fe_index);
    }
}



template <int dim, typename VectorType, int spacedim>
void
SolutionTransfer<dim, VectorType, spacedim>::clear()
{
  // nothing to do
}



namespace Legacy
{

  template <int dim, typename VectorType, int spacedim>
  SolutionTransfer<dim, VectorType, spacedim>::SolutionTransfer(
    const DoFHandler<dim, spacedim> &dof)
    : dof_handler(&dof, typeid(*this).name())
    , n_dofs_old(0)
    , prepared_for(none)
  {
    Assert((dynamic_cast<
              const parallel::distributed::Triangulation<dim, spacedim> *>(
              &dof_handler->get_triangulation()) == nullptr),
           ExcMessage(
             "You are calling the dealii::SolutionTransfer class "
             "with a DoFHandler that is built on a "
             "parallel::distributed::Triangulation. This will not "
             "work for parallel computations. You probably want to "
             "use the parallel::distributed::SolutionTransfer class."));
  }



  template <int dim, typename VectorType, int spacedim>
  SolutionTransfer<dim, VectorType, spacedim>::~SolutionTransfer()
  {
    clear();
  }



  template <int dim, typename VectorType, int spacedim>
  void
  SolutionTransfer<dim, VectorType, spacedim>::clear()
  {
    indices_on_cell.clear();
    dof_values_on_cell.clear();
    cell_map.clear();

    prepared_for = none;
  }



  template <int dim, typename VectorType, int spacedim>
  void
  SolutionTransfer<dim, VectorType, spacedim>::prepare_for_pure_refinement()
  {
    Assert(prepared_for != pure_refinement, ExcAlreadyPrepForRef());
    Assert(prepared_for != coarsening_and_refinement,
           ExcAlreadyPrepForCoarseAndRef());

    clear();

    // We need to access dof indices on the entire domain. For
    // parallel::shared::Triangulations, ownership of cells might change. If
    // they allow artificial cells, we need to restore the "true" cell owners
    // temporarily.
    // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
    // the current set of subdomain ids, set subdomain ids to the "true" owner
    // of each cell upon construction of the TemporarilyRestoreSubdomainIds
    // object, and later restore these flags when it is destroyed.
    const dealii::internal::parallel::shared::
      TemporarilyRestoreSubdomainIds<dim, spacedim>
        subdomain_modifier(dof_handler->get_triangulation());

    const unsigned int n_active_cells =
      dof_handler->get_triangulation().n_active_cells();
    n_dofs_old = dof_handler->n_dofs();

    // efficient reallocation of indices_on_cell
    std::vector<std::vector<types::global_dof_index>>(n_active_cells)
      .swap(indices_on_cell);

    for (const auto &cell : dof_handler->active_cell_iterators())
      {
        const unsigned int i = cell->active_cell_index();
        indices_on_cell[i].resize(cell->get_fe().n_dofs_per_cell());
        // on each cell store the indices of the
        // dofs. after refining we get the values
        // on the children by taking these
        // indices, getting the respective values
        // out of the data vectors and prolonging
        // them to the children
        cell->get_dof_indices(indices_on_cell[i]);
        cell_map[std::make_pair(cell->level(), cell->index())] =
          Pointerstruct(&indices_on_cell[i], cell->active_fe_index());
      }
    prepared_for = pure_refinement;
  }



  template <int dim, typename VectorType, int spacedim>
  void
  SolutionTransfer<dim, VectorType, spacedim>::refine_interpolate(
    const VectorType &in,
    VectorType       &out) const
  {
    Assert(prepared_for == pure_refinement, ExcNotPrepared());
    Assert(in.size() == n_dofs_old,
           ExcDimensionMismatch(in.size(), n_dofs_old));
    Assert(out.size() == dof_handler->n_dofs(),
           ExcDimensionMismatch(out.size(), dof_handler->n_dofs()));
    Assert(&in != &out,
           ExcMessage("Vectors cannot be used as input and output"
                      " at the same time!"));

    // We need to access dof indices on the entire domain. For
    // parallel::shared::Triangulations, ownership of cells might change. If
    // they allow artificial cells, we need to restore the "true" cell owners
    // temporarily.
    // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
    // the current set of subdomain ids, set subdomain ids to the "true" owner
    // of each cell upon construction of the TemporarilyRestoreSubdomainIds
    // object, and later restore these flags when it is destroyed.
    const dealii::internal::parallel::shared::
      TemporarilyRestoreSubdomainIds<dim, spacedim>
        subdomain_modifier(dof_handler->get_triangulation());

    Vector<typename VectorType::value_type> local_values(0);

    typename std::map<std::pair<unsigned int, unsigned int>,
                      Pointerstruct>::const_iterator pointerstruct,
      cell_map_end = cell_map.end();

    for (const auto &cell : dof_handler->cell_iterators())
      {
        pointerstruct =
          cell_map.find(std::make_pair(cell->level(), cell->index()));

        if (pointerstruct != cell_map_end)
          // this cell was refined or not
          // touched at all, so we can get
          // the new values by just setting
          // or interpolating to the children,
          // which is both done by one
          // function
          {
            const unsigned int this_fe_index =
              pointerstruct->second.active_fe_index;
            const unsigned int dofs_per_cell =
              cell->get_dof_handler().get_fe(this_fe_index).n_dofs_per_cell();
            local_values.reinit(dofs_per_cell, true);

            // make sure that the size of the stored indices is the same as
            // dofs_per_cell. since we store the desired fe_index, we know
            // what this size should be
            Assert(dofs_per_cell == (*pointerstruct->second.indices_ptr).size(),
                   ExcInternalError());
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              local_values(i) =
                dealii::internal::ElementAccess<VectorType>::get(
                  in, (*pointerstruct->second.indices_ptr)[i]);
            cell->set_dof_values_by_interpolation(local_values,
                                                  out,
                                                  this_fe_index,
                                                  true);
          }
      }
  }



  namespace internal
  {
    /**
     * Generate a table that contains
     * interpolation matrices between
     * each combination of finite
     * elements used in a DoFHandler of
     * some kind. Since not all
     * elements can be interpolated
     * onto each other, the table may
     * contain empty matrices for those
     * combinations of elements for
     * which no such interpolation is
     * implemented.
     */
    template <int dim, int spacedim>
    void
    extract_interpolation_matrices(
      const DoFHandler<dim, spacedim>      &dof,
      dealii::Table<2, FullMatrix<double>> &matrices)
    {
      if (dof.has_hp_capabilities() == false)
        return;

      const dealii::hp::FECollection<dim, spacedim> &fe =
        dof.get_fe_collection();
      matrices.reinit(fe.size(), fe.size());
      for (unsigned int i = 0; i < fe.size(); ++i)
        for (unsigned int j = 0; j < fe.size(); ++j)
          if (i != j)
            {
              matrices(i, j).reinit(fe[i].n_dofs_per_cell(),
                                    fe[j].n_dofs_per_cell());

              // see if we can get the interpolation matrices for this
              // combination of elements. if not, reset the matrix sizes to zero
              // to indicate that this particular combination isn't
              // supported. this isn't an outright error right away since we may
              // never need to actually interpolate between these two elements
              // on actual cells; we simply have to trigger an error if someone
              // actually tries
              try
                {
                  fe[i].get_interpolation_matrix(fe[j], matrices(i, j));
                }
              catch (const typename FiniteElement<dim, spacedim>::
                       ExcInterpolationNotImplemented &)
                {
                  matrices(i, j).reinit(0, 0);
                }
            }
    }


    template <int dim, int spacedim>
    void
    restriction_additive(const FiniteElement<dim, spacedim> &,
                         std::vector<std::vector<bool>> &)
    {}

    template <int dim, int spacedim>
    void
    restriction_additive(
      const dealii::hp::FECollection<dim, spacedim> &fe,
      std::vector<std::vector<bool>>                &restriction_is_additive)
    {
      restriction_is_additive.resize(fe.size());
      for (unsigned int f = 0; f < fe.size(); ++f)
        {
          restriction_is_additive[f].resize(fe[f].n_dofs_per_cell());
          for (unsigned int i = 0; i < fe[f].n_dofs_per_cell(); ++i)
            restriction_is_additive[f][i] = fe[f].restriction_is_additive(i);
        }
    }
  } // namespace internal



  template <int dim, typename VectorType, int spacedim>
  void
  SolutionTransfer<dim, VectorType, spacedim>::
    prepare_for_coarsening_and_refinement(const std::vector<VectorType> &all_in)
  {
    Assert(prepared_for != pure_refinement, ExcAlreadyPrepForRef());
    Assert(prepared_for != coarsening_and_refinement,
           ExcAlreadyPrepForCoarseAndRef());

    clear();
    n_dofs_old                 = dof_handler->n_dofs();
    const unsigned int in_size = all_in.size();

    if constexpr (running_in_debug_mode())
      {
        Assert(in_size != 0,
               ExcMessage("The array of input vectors you pass to this "
                          "function has no elements. This is not useful."));
        for (unsigned int i = 0; i < in_size; ++i)
          {
            Assert(all_in[i].size() == n_dofs_old,
                   ExcDimensionMismatch(all_in[i].size(), n_dofs_old));
          }
      }

    // We need to access dof indices on the entire domain. For
    // parallel::shared::Triangulations, ownership of cells might change. If
    // they allow artificial cells, we need to restore the "true" cell owners
    // temporarily.
    // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
    // the current set of subdomain ids, set subdomain ids to the "true" owner
    // of each cell upon construction of the TemporarilyRestoreSubdomainIds
    // object, and later restore these flags when it is destroyed.
    const dealii::internal::parallel::shared::
      TemporarilyRestoreSubdomainIds<dim, spacedim>
        subdomain_modifier(dof_handler->get_triangulation());

    // first count the number
    // of cells that will be coarsened
    // and that'll stay or be refined
    unsigned int n_cells_to_coarsen        = 0;
    unsigned int n_cells_to_stay_or_refine = 0;
    for (const auto &act_cell : dof_handler->active_cell_iterators())
      {
        if (act_cell->coarsen_flag_set())
          ++n_cells_to_coarsen;
        else
          ++n_cells_to_stay_or_refine;
      }
    Assert((n_cells_to_coarsen + n_cells_to_stay_or_refine) ==
             dof_handler->get_triangulation().n_active_cells(),
           ExcInternalError());

    unsigned int n_coarsen_fathers = 0;
    for (const auto &cell : dof_handler->cell_iterators())
      if (!cell->is_active() && cell->child(0)->coarsen_flag_set())
        ++n_coarsen_fathers;
    Assert(n_cells_to_coarsen >= 2 * n_coarsen_fathers, ExcInternalError());
    (void)n_cells_to_coarsen;

    // allocate the needed memory. initialize
    // the following arrays in an efficient
    // way, without copying much
    std::vector<std::vector<types::global_dof_index>>(n_cells_to_stay_or_refine)
      .swap(indices_on_cell);

    std::vector<std::vector<Vector<typename VectorType::value_type>>>(
      n_coarsen_fathers,
      std::vector<Vector<typename VectorType::value_type>>(in_size))
      .swap(dof_values_on_cell);

    Table<2, FullMatrix<double>>   interpolation_hp;
    std::vector<std::vector<bool>> restriction_is_additive;

    internal::extract_interpolation_matrices(*dof_handler, interpolation_hp);
    internal::restriction_additive(dof_handler->get_fe_collection(),
                                   restriction_is_additive);

    // we need counters for
    // the 'to_stay_or_refine' cells 'n_sr' and
    // the 'coarsen_fathers' cells 'n_cf',
    unsigned int n_sr = 0, n_cf = 0;
    for (const auto &cell : dof_handler->cell_iterators())
      {
        // CASE 1: active cell that remains as it is
        if (cell->is_active() && !cell->coarsen_flag_set())
          {
            const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
            indices_on_cell[n_sr].resize(dofs_per_cell);
            // cell will not be coarsened,
            // so we get away by storing the
            // dof indices and later
            // interpolating to the children
            cell->get_dof_indices(indices_on_cell[n_sr]);
            cell_map[std::make_pair(cell->level(), cell->index())] =
              Pointerstruct(&indices_on_cell[n_sr], cell->active_fe_index());
            ++n_sr;
          }

        // CASE 2: cell is inactive but will become active
        else if (cell->has_children() && cell->child(0)->coarsen_flag_set())
          {
            // we will need to interpolate from the children of this cell
            // to the current one. in the hp-context, this also means
            // we need to figure out which finite element space to interpolate
            // to since that is not implied by the global FE as in the non-hp-
            // case. we choose the 'least dominant fe' on all children from
            // the associated FECollection.
            std::set<unsigned int> fe_indices_children;
            for (const auto &child : cell->child_iterators())
              {
                Assert(child->is_active() && child->coarsen_flag_set(),
                       typename dealii::Triangulation<
                         dim>::ExcInconsistentCoarseningFlags());

                fe_indices_children.insert(child->active_fe_index());
              }
            Assert(!fe_indices_children.empty(), ExcInternalError());

            const unsigned int target_fe_index =
              dof_handler->get_fe_collection().find_dominated_fe_extended(
                fe_indices_children, /*codim=*/0);

            Assert(target_fe_index != numbers::invalid_unsigned_int,
                   dealii::internal::hp::DoFHandlerImplementation::
                     ExcNoDominatedFiniteElementOnChildren());

            const unsigned int dofs_per_cell =
              dof_handler->get_fe(target_fe_index).n_dofs_per_cell();

            std::vector<Vector<typename VectorType::value_type>>(
              in_size, Vector<typename VectorType::value_type>(dofs_per_cell))
              .swap(dof_values_on_cell[n_cf]);


            // store the data of each of the input vectors. get this data
            // as interpolated onto a finite element space that encompasses
            // that of all the children. note that
            // cell->get_interpolated_dof_values already does all of the
            // interpolations between spaces
            for (unsigned int j = 0; j < in_size; ++j)
              cell->get_interpolated_dof_values(all_in[j],
                                                dof_values_on_cell[n_cf][j],
                                                target_fe_index);
            cell_map[std::make_pair(cell->level(), cell->index())] =
              Pointerstruct(&dof_values_on_cell[n_cf], target_fe_index);
            ++n_cf;
          }
      }
    Assert(n_sr == n_cells_to_stay_or_refine, ExcInternalError());
    Assert(n_cf == n_coarsen_fathers, ExcInternalError());

    prepared_for = coarsening_and_refinement;
  }



  template <int dim, typename VectorType, int spacedim>
  void
  SolutionTransfer<dim, VectorType, spacedim>::
    prepare_for_coarsening_and_refinement(const VectorType &in)
  {
    std::vector<VectorType> all_in(1, in);
    prepare_for_coarsening_and_refinement(all_in);
  }



  template <int dim, typename VectorType, int spacedim>
  void
  SolutionTransfer<dim, VectorType, spacedim>::interpolate(
    const std::vector<VectorType> &all_in,
    std::vector<VectorType>       &all_out) const
  {
    const unsigned int size = all_in.size();
    if constexpr (running_in_debug_mode())
      {
        Assert(prepared_for == coarsening_and_refinement, ExcNotPrepared());
        Assert(all_out.size() == size,
               ExcDimensionMismatch(all_out.size(), size));
        for (unsigned int i = 0; i < size; ++i)
          Assert(all_in[i].size() == n_dofs_old,
                 ExcDimensionMismatch(all_in[i].size(), n_dofs_old));
        for (unsigned int i = 0; i < all_out.size(); ++i)
          Assert(all_out[i].size() == dof_handler->n_dofs(),
                 ExcDimensionMismatch(all_out[i].size(),
                                      dof_handler->n_dofs()));
        for (unsigned int i = 0; i < size; ++i)
          for (unsigned int j = 0; j < size; ++j)
            Assert(&all_in[i] != &all_out[j],
                   ExcMessage("Vectors cannot be used as input and output"
                              " at the same time!"));
      }

    // We need to access dof indices on the entire domain. For
    // parallel::shared::Triangulations, ownership of cells might change. If
    // they allow artificial cells, we need to restore the "true" cell owners
    // temporarily.
    // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
    // the current set of subdomain ids, set subdomain ids to the "true" owner
    // of each cell upon construction of the TemporarilyRestoreSubdomainIds
    // object, and later restore these flags when it is destroyed.
    const dealii::internal::parallel::shared::
      TemporarilyRestoreSubdomainIds<dim, spacedim>
        subdomain_modifier(dof_handler->get_triangulation());

    Vector<typename VectorType::value_type> local_values;
    std::vector<types::global_dof_index>    dofs;

    typename std::map<std::pair<unsigned int, unsigned int>,
                      Pointerstruct>::const_iterator pointerstruct,
      cell_map_end = cell_map.end();

    Table<2, FullMatrix<double>> interpolation_hp;
    internal::extract_interpolation_matrices(*dof_handler, interpolation_hp);
    Vector<typename VectorType::value_type> tmp, tmp2;

    for (const auto &cell : dof_handler->cell_iterators())
      {
        pointerstruct =
          cell_map.find(std::make_pair(cell->level(), cell->index()));

        if (pointerstruct != cell_map_end)
          {
            const std::vector<types::global_dof_index> *const indexptr =
              pointerstruct->second.indices_ptr;

            const std::vector<Vector<typename VectorType::value_type>>
              *const valuesptr = pointerstruct->second.dof_values_ptr;

            // cell stayed as it was or was refined
            if (indexptr != nullptr)
              {
                Assert(valuesptr == nullptr, ExcInternalError());

                const unsigned int old_fe_index =
                  pointerstruct->second.active_fe_index;

                // get the values of each of the input data vectors on this cell
                // and prolong it to its children
                unsigned int in_size = indexptr->size();
                for (unsigned int j = 0; j < size; ++j)
                  {
                    tmp.reinit(in_size, true);
                    for (unsigned int i = 0; i < in_size; ++i)
                      tmp(i) = dealii::internal::ElementAccess<VectorType>::get(
                        all_in[j], (*indexptr)[i]);

                    cell->set_dof_values_by_interpolation(tmp,
                                                          all_out[j],
                                                          old_fe_index,
                                                          true);
                  }
              }
            else if (valuesptr)
              // the children of this cell were deleted
              {
                Assert(!cell->has_children(), ExcInternalError());
                Assert(indexptr == nullptr, ExcInternalError());

                const unsigned int dofs_per_cell =
                  cell->get_fe().n_dofs_per_cell();
                dofs.resize(dofs_per_cell);
                // get the local
                // indices
                cell->get_dof_indices(dofs);

                // distribute the stored data to the new vectors
                for (unsigned int j = 0; j < size; ++j)
                  {
                    // make sure that the size of the stored indices is the same
                    // as dofs_per_cell. this is kind of a test if we use the
                    // same FE in the hp-case. to really do that test we would
                    // have to store the fe_index of all cells
                    const Vector<typename VectorType::value_type> *data =
                      nullptr;
                    const unsigned int active_fe_index =
                      cell->active_fe_index();
                    if (active_fe_index !=
                        pointerstruct->second.active_fe_index)
                      {
                        const unsigned int old_index =
                          pointerstruct->second.active_fe_index;
                        const FullMatrix<double> &interpolation_matrix =
                          interpolation_hp(active_fe_index, old_index);
                        // The interpolation matrix might be empty when using
                        // FE_Nothing.
                        if (interpolation_matrix.empty())
                          tmp.reinit(dofs_per_cell, false);
                        else
                          {
                            tmp.reinit(dofs_per_cell, true);
                            AssertDimension((*valuesptr)[j].size(),
                                            interpolation_matrix.n());
                            AssertDimension(tmp.size(),
                                            interpolation_matrix.m());
                            interpolation_matrix.vmult(tmp, (*valuesptr)[j]);
                          }
                        data = &tmp;
                      }
                    else
                      data = &(*valuesptr)[j];


                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      dealii::internal::ElementAccess<VectorType>::set(
                        (*data)(i), dofs[i], all_out[j]);
                  }
              }
            // undefined status
            else
              Assert(false, ExcInternalError());
          }
      }

    // We have written into the output vectors. If this was a PETSc vector, for
    // example, then we need to compress these to make future operations safe:
    for (auto &vec : all_out)
      vec.compress(VectorOperation::insert);
  }



  template <int dim, typename VectorType, int spacedim>
  void
  SolutionTransfer<dim, VectorType, spacedim>::interpolate(
    const VectorType &in,
    VectorType       &out) const
  {
    Assert(in.size() == n_dofs_old,
           ExcDimensionMismatch(in.size(), n_dofs_old));
    Assert(out.size() == dof_handler->n_dofs(),
           ExcDimensionMismatch(out.size(), dof_handler->n_dofs()));

    std::vector<VectorType> all_in  = {in};
    std::vector<VectorType> all_out = {out};

    interpolate(all_in, all_out);

    out = all_out[0];
  }



  template <int dim, typename VectorType, int spacedim>
  std::size_t
  SolutionTransfer<dim, VectorType, spacedim>::memory_consumption() const
  {
    // at the moment we do not include the memory
    // consumption of the cell_map as we have no
    // real idea about memory consumption of a
    // std::map
    return (MemoryConsumption::memory_consumption(dof_handler) +
            MemoryConsumption::memory_consumption(n_dofs_old) +
            sizeof(prepared_for) +
            MemoryConsumption::memory_consumption(indices_on_cell) +
            MemoryConsumption::memory_consumption(dof_values_on_cell));
  }



  template <int dim, typename VectorType, int spacedim>
  std::size_t
  SolutionTransfer<dim, VectorType, spacedim>::Pointerstruct::
    memory_consumption() const
  {
    return sizeof(*this);
  }

} // namespace Legacy

DEAL_II_NAMESPACE_CLOSE


#endif
