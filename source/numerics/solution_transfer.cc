// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#include <deal.II/base/memory_consumption.h>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_element_access.h>

#include <deal.II/numerics/solution_transfer.h>

DEAL_II_NAMESPACE_OPEN

template <int dim, typename VectorType, typename DoFHandlerType>
SolutionTransfer<dim, VectorType, DoFHandlerType>::SolutionTransfer(
  const DoFHandlerType &dof)
  : dof_handler(&dof, typeid(*this).name())
  , n_dofs_old(0)
  , prepared_for(none)
{
  Assert((dynamic_cast<const parallel::distributed::Triangulation<
            DoFHandlerType::dimension,
            DoFHandlerType::space_dimension> *>(
            &dof_handler->get_triangulation()) == nullptr),
         ExcMessage("You are calling the dealii::SolutionTransfer class "
                    "with a DoF handler that is built on a "
                    "parallel::distributed::Triangulation. This will not "
                    "work for parallel computations. You probably want to "
                    "use the parallel::distributed::SolutionTransfer class."));
}



template <int dim, typename VectorType, typename DoFHandlerType>
SolutionTransfer<dim, VectorType, DoFHandlerType>::~SolutionTransfer()
{
  clear();
}



template <int dim, typename VectorType, typename DoFHandlerType>
void
SolutionTransfer<dim, VectorType, DoFHandlerType>::clear()
{
  indices_on_cell.clear();
  dof_values_on_cell.clear();
  cell_map.clear();

  prepared_for = none;
}



template <int dim, typename VectorType, typename DoFHandlerType>
void
SolutionTransfer<dim, VectorType, DoFHandlerType>::prepare_for_pure_refinement()
{
  Assert(prepared_for != pure_refinement, ExcAlreadyPrepForRef());
  Assert(prepared_for != coarsening_and_refinement,
         ExcAlreadyPrepForCoarseAndRef());

  clear();

  // We need to access dof indices on the entire domain. For
  // parallel::shared::Triangulations, ownership of cells might change. If they
  // allow artificial cells, we need to restore the "true" cell owners
  // temporarily.
  // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
  // the current set of subdomain ids, set subdomain ids to the "true" owner of
  // each cell upon construction of the TemporarilyRestoreSubdomainIds object,
  // and later restore these flags when it is destroyed.
  const internal::parallel::shared::
    TemporarilyRestoreSubdomainIds<dim, DoFHandlerType::space_dimension>
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



template <int dim, typename VectorType, typename DoFHandlerType>
void
SolutionTransfer<dim, VectorType, DoFHandlerType>::refine_interpolate(
  const VectorType &in,
  VectorType &      out) const
{
  Assert(prepared_for == pure_refinement, ExcNotPrepared());
  Assert(in.size() == n_dofs_old, ExcDimensionMismatch(in.size(), n_dofs_old));
  Assert(out.size() == dof_handler->n_dofs(),
         ExcDimensionMismatch(out.size(), dof_handler->n_dofs()));
  Assert(&in != &out,
         ExcMessage("Vectors cannot be used as input and output"
                    " at the same time!"));

  // We need to access dof indices on the entire domain. For
  // parallel::shared::Triangulations, ownership of cells might change. If they
  // allow artificial cells, we need to restore the "true" cell owners
  // temporarily.
  // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
  // the current set of subdomain ids, set subdomain ids to the "true" owner of
  // each cell upon construction of the TemporarilyRestoreSubdomainIds object,
  // and later restore these flags when it is destroyed.
  const internal::parallel::shared::
    TemporarilyRestoreSubdomainIds<dim, DoFHandlerType::space_dimension>
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
            local_values(i) = internal::ElementAccess<VectorType>::get(
              in, (*pointerstruct->second.indices_ptr)[i]);
          cell->set_dof_values_by_interpolation(local_values,
                                                out,
                                                this_fe_index);
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
  extract_interpolation_matrices(const dealii::DoFHandler<dim, spacedim> &dof,
                                 dealii::Table<2, FullMatrix<double>> &matrices)
  {
    if (dof.has_hp_capabilities() == false)
      return;

    const dealii::hp::FECollection<dim, spacedim> &fe = dof.get_fe_collection();
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
  restriction_additive(const dealii::hp::FECollection<dim, spacedim> &fe,
                       std::vector<std::vector<bool>> &restriction_is_additive)
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



template <int dim, typename VectorType, typename DoFHandlerType>
void
SolutionTransfer<dim, VectorType, DoFHandlerType>::
  prepare_for_coarsening_and_refinement(const std::vector<VectorType> &all_in)
{
  Assert(prepared_for != pure_refinement, ExcAlreadyPrepForRef());
  Assert(prepared_for != coarsening_and_refinement,
         ExcAlreadyPrepForCoarseAndRef());

  clear();
  n_dofs_old                 = dof_handler->n_dofs();
  const unsigned int in_size = all_in.size();

#ifdef DEBUG
  Assert(in_size != 0,
         ExcMessage("The array of input vectors you pass to this "
                    "function has no elements. This is not useful."));
  for (unsigned int i = 0; i < in_size; ++i)
    {
      Assert(all_in[i].size() == n_dofs_old,
             ExcDimensionMismatch(all_in[i].size(), n_dofs_old));
    }
#endif

  // We need to access dof indices on the entire domain. For
  // parallel::shared::Triangulations, ownership of cells might change. If they
  // allow artificial cells, we need to restore the "true" cell owners
  // temporarily.
  // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
  // the current set of subdomain ids, set subdomain ids to the "true" owner of
  // each cell upon construction of the TemporarilyRestoreSubdomainIds object,
  // and later restore these flags when it is destroyed.
  const internal::parallel::shared::
    TemporarilyRestoreSubdomainIds<dim, DoFHandlerType::space_dimension>
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
                 internal::hp::DoFHandlerImplementation::
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



template <int dim, typename VectorType, typename DoFHandlerType>
void
SolutionTransfer<dim, VectorType, DoFHandlerType>::
  prepare_for_coarsening_and_refinement(const VectorType &in)
{
  std::vector<VectorType> all_in(1, in);
  prepare_for_coarsening_and_refinement(all_in);
}



template <int dim, typename VectorType, typename DoFHandlerType>
void
SolutionTransfer<dim, VectorType, DoFHandlerType>::interpolate(
  const std::vector<VectorType> &all_in,
  std::vector<VectorType> &      all_out) const
{
  const unsigned int size = all_in.size();
#ifdef DEBUG
  Assert(prepared_for == coarsening_and_refinement, ExcNotPrepared());
  Assert(all_out.size() == size, ExcDimensionMismatch(all_out.size(), size));
  for (unsigned int i = 0; i < size; ++i)
    Assert(all_in[i].size() == n_dofs_old,
           ExcDimensionMismatch(all_in[i].size(), n_dofs_old));
  for (unsigned int i = 0; i < all_out.size(); ++i)
    Assert(all_out[i].size() == dof_handler->n_dofs(),
           ExcDimensionMismatch(all_out[i].size(), dof_handler->n_dofs()));
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      Assert(&all_in[i] != &all_out[j],
             ExcMessage("Vectors cannot be used as input and output"
                        " at the same time!"));
#endif

  // We need to access dof indices on the entire domain. For
  // parallel::shared::Triangulations, ownership of cells might change. If they
  // allow artificial cells, we need to restore the "true" cell owners
  // temporarily.
  // We use the TemporarilyRestoreSubdomainIds class for this purpose: we save
  // the current set of subdomain ids, set subdomain ids to the "true" owner of
  // each cell upon construction of the TemporarilyRestoreSubdomainIds object,
  // and later restore these flags when it is destroyed.
  const internal::parallel::shared::
    TemporarilyRestoreSubdomainIds<dim, DoFHandlerType::space_dimension>
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
          if (indexptr)
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
                    tmp(i) =
                      internal::ElementAccess<VectorType>::get(all_in[j],
                                                               (*indexptr)[i]);

                  cell->set_dof_values_by_interpolation(tmp,
                                                        all_out[j],
                                                        old_fe_index);
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
                  // as dofs_per_cell. this is kind of a test if we use the same
                  // FE in the hp-case. to really do that test we would have to
                  // store the fe_index of all cells
                  const Vector<typename VectorType::value_type> *data = nullptr;
                  const unsigned int active_fe_index = cell->active_fe_index();
                  if (active_fe_index != pointerstruct->second.active_fe_index)
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
                          AssertDimension(tmp.size(), interpolation_matrix.m());
                          interpolation_matrix.vmult(tmp, (*valuesptr)[j]);
                        }
                      data = &tmp;
                    }
                  else
                    data = &(*valuesptr)[j];


                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    internal::ElementAccess<VectorType>::set((*data)(i),
                                                             dofs[i],
                                                             all_out[j]);
                }
            }
          // undefined status
          else
            Assert(false, ExcInternalError());
        }
    }
}



template <int dim, typename VectorType, typename DoFHandlerType>
void
SolutionTransfer<dim, VectorType, DoFHandlerType>::interpolate(
  const VectorType &in,
  VectorType &      out) const
{
  Assert(in.size() == n_dofs_old, ExcDimensionMismatch(in.size(), n_dofs_old));
  Assert(out.size() == dof_handler->n_dofs(),
         ExcDimensionMismatch(out.size(), dof_handler->n_dofs()));

  std::vector<VectorType> all_in(1);
  all_in[0] = in;
  std::vector<VectorType> all_out(1);
  all_out[0] = out;
  interpolate(all_in, all_out);
  out = all_out[0];
}



template <int dim, typename VectorType, typename DoFHandlerType>
std::size_t
SolutionTransfer<dim, VectorType, DoFHandlerType>::memory_consumption() const
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



template <int dim, typename VectorType, typename DoFHandlerType>
std::size_t
SolutionTransfer<dim, VectorType, DoFHandlerType>::Pointerstruct::
  memory_consumption() const
{
  return sizeof(*this);
}


/*-------------- Explicit Instantiations -------------------------------*/
#define SPLIT_INSTANTIATIONS_COUNT 4
#ifndef SPLIT_INSTANTIATIONS_INDEX
#  define SPLIT_INSTANTIATIONS_INDEX 0
#endif
#include "solution_transfer.inst"

DEAL_II_NAMESPACE_CLOSE
