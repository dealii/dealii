// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/numerics/solution_transfer.h>

DEAL_II_NAMESPACE_OPEN



template<int dim, typename VECTOR, class DH>
SolutionTransfer<dim, VECTOR, DH>::SolutionTransfer(const DH &dof)
  :
  dof_handler(&dof, typeid(*this).name()),
  n_dofs_old(0),
  prepared_for(none)
{}



template<int dim, typename VECTOR, class DH>
SolutionTransfer<dim, VECTOR, DH>::~SolutionTransfer()
{
  clear ();
}



template<int dim, typename VECTOR, class DH>
void SolutionTransfer<dim, VECTOR, DH>::clear ()
{
  indices_on_cell.clear();
  dof_values_on_cell.clear();
  cell_map.clear();

  prepared_for=none;
}



template<int dim, typename VECTOR, class DH>
void SolutionTransfer<dim, VECTOR, DH>::prepare_for_pure_refinement()
{
  Assert(prepared_for!=pure_refinement, ExcAlreadyPrepForRef());
  Assert(prepared_for!=coarsening_and_refinement,
         ExcAlreadyPrepForCoarseAndRef());

  clear();

  const unsigned int n_active_cells = dof_handler->get_tria().n_active_cells();
  n_dofs_old=dof_handler->n_dofs();

  // efficient reallocation of indices_on_cell
  std::vector<std::vector<types::global_dof_index> > (n_active_cells)
  .swap(indices_on_cell);

  typename DH::active_cell_iterator cell = dof_handler->begin_active(),
                                    endc = dof_handler->end();

  for (unsigned int i=0; cell!=endc; ++cell, ++i)
    {
      indices_on_cell[i].resize(cell->get_fe().dofs_per_cell);
      // on each cell store the indices of the
      // dofs. after refining we get the values
      // on the children by taking these
      // indices, getting the respective values
      // out of the data vectors and prolonging
      // them to the children
      cell->get_dof_indices(indices_on_cell[i]);
      cell_map[std::make_pair(cell->level(),cell->index())].indices_ptr=&indices_on_cell[i];
    }
  prepared_for=pure_refinement;
}



template<int dim, typename VECTOR, class DH>
void
SolutionTransfer<dim, VECTOR, DH>::refine_interpolate(const VECTOR &in,
                                                      VECTOR       &out) const
{
  Assert(prepared_for==pure_refinement, ExcNotPrepared());
  Assert(in.size()==n_dofs_old, ExcDimensionMismatch(in.size(),n_dofs_old));
  Assert(out.size()==dof_handler->n_dofs(),
         ExcDimensionMismatch(out.size(),dof_handler->n_dofs()));
  Assert(&in != &out,
         ExcMessage ("Vectors cannot be used as input and output"
                     " at the same time!"));

  Vector<typename VECTOR::value_type> local_values(0);

  typename DH::cell_iterator cell = dof_handler->begin(),
                             endc = dof_handler->end();

  typename std::map<std::pair<unsigned int, unsigned int>, Pointerstruct>::const_iterator
  pointerstruct,
  cell_map_end=cell_map.end();

  for (; cell!=endc; ++cell)
    {
      pointerstruct=cell_map.find(std::make_pair(cell->level(),cell->index()));

      if (pointerstruct!=cell_map_end)
        // this cell was refined or not
        // touched at all, so we can get
        // the new values by just setting
        // or interpolating to the children,
        // which is both done by one
        // function
        {
          const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;
          local_values.reinit(dofs_per_cell, true); // fast reinit, i.e. the
          // entries are not set to 0.
          // make sure that the size of the
          // stored indices is the same as
          // dofs_per_cell. this is kind of a
          // test if we use the same fe in the
          // hp case. to really do that test we
          // would have to store the fe_index
          // of all cells
          Assert(dofs_per_cell==(*pointerstruct->second.indices_ptr).size(),
                 ExcNumberOfDoFsPerCellHasChanged());
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            local_values(i)=in((*pointerstruct->second.indices_ptr)[i]);
          cell->set_dof_values_by_interpolation(local_values, out);
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
  template <typename DH>
  void extract_interpolation_matrices (const DH &,
                                       Table<2,FullMatrix<double> > &)
  {}

  template <int dim, int spacedim>
  void extract_interpolation_matrices (const dealii::hp::DoFHandler<dim,spacedim> &dof,
                                       Table<2,FullMatrix<double> > &matrices)
  {
    const dealii::hp::FECollection<dim,spacedim> &fe = dof.get_fe();
    matrices.reinit (fe.size(), fe.size());
    for (unsigned int i=0; i<fe.size(); ++i)
      for (unsigned int j=0; j<fe.size(); ++j)
        if (i != j)
          {
            matrices(i,j).reinit (fe[i].dofs_per_cell, fe[j].dofs_per_cell);

            // see if we can get the interpolation matrices for this
            // combination of elements. if not, reset the matrix sizes to zero
            // to indicate that this particular combination isn't
            // supported. this isn't an outright error right away since we may
            // never need to actually interpolate between these two elements
            // on actual cells; we simply have to trigger an error if someone
            // actually tries
            try
              {
                fe[i].get_interpolation_matrix (fe[j], matrices(i,j));
              }
            catch (const typename FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented &)
              {
                matrices(i,j).reinit (0,0);
              }
          }
  }


  template <int dim, int spacedim>
  void restriction_additive (const FiniteElement<dim,spacedim> &,
                             std::vector<std::vector<bool> > &)
  {}

  template <int dim, int spacedim>
  void restriction_additive (const dealii::hp::FECollection<dim,spacedim> &fe,
                             std::vector<std::vector<bool> > &restriction_is_additive)
  {
    restriction_is_additive.resize (fe.size());
    for (unsigned int f=0; f<fe.size(); ++f)
      {
        restriction_is_additive[f].resize (fe[f].dofs_per_cell);
        for (unsigned int i=0; i<fe[f].dofs_per_cell; ++i)
          restriction_is_additive[f][i] = fe[f].restriction_is_additive(i);
      }
  }
}



template<int dim, typename VECTOR, class DH>
void
SolutionTransfer<dim, VECTOR, DH>::
prepare_for_coarsening_and_refinement(const std::vector<VECTOR> &all_in)
{
  Assert(prepared_for!=pure_refinement, ExcAlreadyPrepForRef());
  Assert(prepared_for!=coarsening_and_refinement,
         ExcAlreadyPrepForCoarseAndRef());

  const unsigned int in_size=all_in.size();
  Assert(in_size!=0, ExcNoInVectorsGiven());

  clear();

  const unsigned int n_active_cells = dof_handler->get_tria().n_active_cells();
  n_dofs_old=dof_handler->n_dofs();

  for (unsigned int i=0; i<in_size; ++i)
    {
      Assert(all_in[i].size()==n_dofs_old,
             ExcDimensionMismatch(all_in[i].size(),n_dofs_old));
    }

  // first count the number
  // of cells that'll be coarsened
  // and that'll stay or be refined
  unsigned int n_cells_to_coarsen=0;
  unsigned int n_cells_to_stay_or_refine=0;
  typename DH::active_cell_iterator
  act_cell = dof_handler->begin_active(),
  endc = dof_handler->end();
  for (; act_cell!=endc; ++act_cell)
    {
      if (act_cell->coarsen_flag_set())
        ++n_cells_to_coarsen;
      else
        ++n_cells_to_stay_or_refine;
    }
  Assert((n_cells_to_coarsen+n_cells_to_stay_or_refine)==n_active_cells,
         ExcInternalError());

  unsigned int n_coarsen_fathers=0;
  typename DH::cell_iterator
  cell=dof_handler->begin();
  for (; cell!=endc; ++cell)
    if (!cell->active() && cell->child(0)->coarsen_flag_set())
      ++n_coarsen_fathers;

  if (n_cells_to_coarsen)
    Assert(n_cells_to_coarsen>=2*n_coarsen_fathers, ExcInternalError());

  // allocate the needed memory. initialize
  // the following arrays in an efficient
  // way, without copying much
  std::vector<std::vector<types::global_dof_index> >
  (n_cells_to_stay_or_refine)
  .swap(indices_on_cell);

  std::vector<std::vector<Vector<typename VECTOR::value_type> > >
  (n_coarsen_fathers,
   std::vector<Vector<typename VECTOR::value_type> > (in_size))
  .swap(dof_values_on_cell);

  typename VECTOR::value_type zero_val = typename VECTOR::value_type();
  Table<2,FullMatrix<double> > interpolation_hp;
  std::vector<std::vector<bool> > restriction_is_additive;
  internal::extract_interpolation_matrices (*dof_handler, interpolation_hp);
  internal::restriction_additive (dof_handler->get_fe(), restriction_is_additive);
  Vector<typename VECTOR::value_type> tmp, tmp2;

  // we need counters for
  // the 'to_stay_or_refine' cells 'n_sr' and
  // the 'coarsen_fathers' cells 'n_cf',
  unsigned int n_sr=0, n_cf=0;
  cell = dof_handler->begin();
  for (; cell!=endc; ++cell)
    {
      if (cell->active() && !cell->coarsen_flag_set())
        {
          const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;
          indices_on_cell[n_sr].resize(dofs_per_cell);
          // cell will not be coarsened,
          // so we get away by storing the
          // dof indices and later
          // interpolating to the children
          cell->get_dof_indices(indices_on_cell[n_sr]);
          cell_map[std::make_pair(cell->level(), cell->index())]
            = Pointerstruct(&indices_on_cell[n_sr],cell->active_fe_index());
          ++n_sr;
        }
      else if (cell->has_children() && cell->child(0)->coarsen_flag_set())
        {
          // note that if one child has the
          // coarsen flag, then all should
          // have if Tria::prepare_* has
          // worked correctly
          for (unsigned int i=1; i<cell->n_children(); ++i)
            Assert(cell->child(i)->coarsen_flag_set(),
                   ExcTriaPrepCoarseningNotCalledBefore());
          const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;

          std::vector<Vector<typename VECTOR::value_type> >(in_size,
                                                            Vector<typename VECTOR::value_type>(dofs_per_cell))
          .swap(dof_values_on_cell[n_cf]);

          unsigned int fe_index = cell->active_fe_index();
          unsigned int most_general_child = 0;
          bool different_elements = false;
          for (unsigned int child=0; child<cell->n_children(); ++child)
            {
              if (cell->child(child)->active_fe_index() != fe_index)
                different_elements = true;
              // take FE index from the child with most
              // degrees of freedom locally
              if (cell->child(child)->get_fe().dofs_per_cell >
                  cell->child(most_general_child)->get_fe().dofs_per_cell)
                most_general_child = child;
            }
          if (different_elements == true)
            fe_index = cell->child(most_general_child)->active_fe_index();

          for (unsigned int j=0; j<in_size; ++j)
            {
              // store the data of each of
              // the input vectors
              if (different_elements == false)
                cell->get_interpolated_dof_values(all_in[j],
                                                  dof_values_on_cell[n_cf][j]);
              else
                {
                  // if we have different elements, first
                  // interpolate the children's contribution to
                  // the most general FE on the child
                  // level. Then we manually write the
                  // interpolation operation to the coarser
                  // level
                  dof_values_on_cell[n_cf][j].reinit (cell->child(most_general_child)->get_fe().dofs_per_cell);
                  const unsigned int fe_ind_general =
                    cell->child(most_general_child)->active_fe_index();
                  for (unsigned int child=0; child<cell->n_children(); ++child)
                    {
                      tmp.reinit (cell->child(child)->get_fe().dofs_per_cell,
                                  true);
                      cell->child(child)->get_dof_values (all_in[j],
                                                          tmp);
                      const unsigned int child_ind =
                        cell->child(child)->active_fe_index();
                      if (child_ind != fe_ind_general)
                        {
                          tmp2.reinit (cell->child(most_general_child)->get_fe().dofs_per_cell,
                                       true);

                          // if the
                          // matrix
                          // has size
                          // 0 and
                          // the
                          // corresponding
                          // elements
                          // have
                          // more
                          // than
                          // zero
                          // DoFs,
                          // then
                          // this
                          // means
                          // that the
                          // internal::extract_interpolation_matrices
                          // function
                          // above
                          // couldn't
                          // get the
                          // interpolation
                          // matrix
                          // between
                          // this
                          // pair of
                          // elements. since
                          // we need
                          // it here,
                          // this is
                          // an error
                          if ((interpolation_hp(fe_ind_general, child_ind).m() == 0)
                              &&
                              (interpolation_hp(fe_ind_general, child_ind).n() == 0)
                              &&
                              ((tmp2.size() > 0) || (tmp.size() > 0)))
                            Assert (false,
                                    ExcMessage (std::string("No interpolation from element ")
                                                + cell->child(child)->get_fe().get_name()
                                                + " to element "
                                                + cell->child(most_general_child)->get_fe().get_name()
                                                + " is available, but it is needed here."));
                          interpolation_hp(fe_ind_general, child_ind).vmult (tmp2, tmp);
                        }
                      else
                        tmp2.swap (tmp);

                      // now do the interpolation
                      // operation. FullMatrix only wants us
                      // to call vmult if the matrix size is
                      // actually non-zero, so check that
                      // case
                      const unsigned int dofs_per_cell = tmp2.size();
                      tmp.reinit (dofs_per_cell, true);
                      if (dofs_per_cell > 0)
                        cell->child(most_general_child)->get_fe().
                        get_restriction_matrix(child, cell->refinement_case()).vmult (tmp, tmp2);
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        if (restriction_is_additive[fe_ind_general][i])
                          dof_values_on_cell[n_cf][j](i) += tmp(i);
                        else if (tmp(i) != zero_val)
                          dof_values_on_cell[n_cf][j](i) = tmp(i);
                    }
                }
            }
          cell_map[std::make_pair(cell->level(), cell->index())]
            = Pointerstruct(&dof_values_on_cell[n_cf], fe_index);
          ++n_cf;
        }
    }
  Assert(n_sr==n_cells_to_stay_or_refine, ExcInternalError());
  Assert(n_cf==n_coarsen_fathers, ExcInternalError());

  prepared_for=coarsening_and_refinement;
}



template<int dim, typename VECTOR, class DH>
void
SolutionTransfer<dim, VECTOR, DH>::prepare_for_coarsening_and_refinement(const VECTOR &in)
{
  std::vector<VECTOR> all_in=std::vector<VECTOR>(1, in);
  prepare_for_coarsening_and_refinement(all_in);
}



template<int dim, typename VECTOR, class DH>
void SolutionTransfer<dim, VECTOR, DH>::
interpolate (const std::vector<VECTOR> &all_in,
             std::vector<VECTOR>       &all_out) const
{
  Assert(prepared_for==coarsening_and_refinement, ExcNotPrepared());
  const unsigned int size=all_in.size();
  Assert(all_out.size()==size, ExcDimensionMismatch(all_out.size(), size));
  for (unsigned int i=0; i<size; ++i)
    Assert (all_in[i].size() == n_dofs_old,
            ExcDimensionMismatch(all_in[i].size(), n_dofs_old));
  for (unsigned int i=0; i<all_out.size(); ++i)
    Assert (all_out[i].size() == dof_handler->n_dofs(),
            ExcDimensionMismatch(all_out[i].size(), dof_handler->n_dofs()));
  for (unsigned int i=0; i<size; ++i)
    for (unsigned int j=0; j<size; ++j)
      Assert(&all_in[i] != &all_out[j],
             ExcMessage ("Vectors cannot be used as input and output"
                         " at the same time!"));

  Vector<typename VECTOR::value_type> local_values;
  std::vector<types::global_dof_index> dofs;

  typename std::map<std::pair<unsigned int, unsigned int>, Pointerstruct>::const_iterator
  pointerstruct,
  cell_map_end=cell_map.end();

  Table<2,FullMatrix<double> > interpolation_hp;
  internal::extract_interpolation_matrices (*dof_handler, interpolation_hp);
  Vector<typename VECTOR::value_type> tmp, tmp2;

  typename DH::cell_iterator cell = dof_handler->begin(),
                             endc = dof_handler->end();
  for (; cell!=endc; ++cell)
    {
      pointerstruct=cell_map.find(std::make_pair(cell->level(),cell->index()));

      if (pointerstruct!=cell_map_end)
        {
          const std::vector<types::global_dof_index> *const indexptr
            =pointerstruct->second.indices_ptr;

          const std::vector<Vector<typename VECTOR::value_type> > *const valuesptr
            =pointerstruct->second.dof_values_ptr;

          const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;

          // cell stayed or is
          // refined
          if (indexptr)
            {
              Assert (valuesptr == 0,
                      ExcInternalError());

              // get the values of
              // each of the input
              // data vectors on this
              // cell and prolong it
              // to its children
              unsigned int in_size = indexptr->size();
              local_values.reinit(dofs_per_cell, true);
              for (unsigned int j=0; j<size; ++j)
                {
                  // check the FE index of the new element. if
                  // we have children and not all of the
                  // children have the same FE index, need to
                  // manually implement
                  // set_dof_values_by_interpolation
                  unsigned int new_fe_index = cell->active_fe_index();
                  bool different_elements = false;
                  if (cell->has_children())
                    {
                      new_fe_index = cell->child(0)->active_fe_index();
                      for (unsigned int child=1; child<cell->n_children(); ++child)
                        if (cell->child(child)->active_fe_index() != new_fe_index)
                          {
                            different_elements = true;
                            break;
                          }
                    }
                  if (different_elements == true ||
                      new_fe_index != pointerstruct->second.active_fe_index)
                    {
                      const unsigned int old_index =
                        pointerstruct->second.active_fe_index;
                      tmp.reinit (in_size, true);
                      for (unsigned int i=0; i<in_size; ++i)
                        tmp(i)=all_in[j]((*indexptr)[i]);
                      if (different_elements == false)
                        {
                          AssertDimension (tmp.size(),
                                           interpolation_hp(new_fe_index,old_index).n());
                          local_values.reinit (cell->has_children() ?
                                               cell->child(0)->get_fe().dofs_per_cell
                                               : cell->get_fe().dofs_per_cell, true);
                          // do the interpolation. we get into trouble if the
                          // interpolation_hp(new,old) matrix hasn't been computed.
                          // this can happen if the respective elements don't support
                          // the corresponding interpolation; if that's the case, then
                          // the computation of the matrix simply sets the matrix
                          // back to size zero. so if we get here and that is
                          // the wrong size, then this may be because the elements
                          // haven't implemented the correct function yet
                          //
                          // there is one wrinkle. we would like to only error out if
                          // the size of the matrix is 0 times 0 but at least one
                          // of the elements has more than one dof per cell. the
                          // problem is that if you reinit a matrix to 4x0, it automatically
                          // sets its size to 0x0. so we can only execute the following
                          // test if *both* elements have dofs_per_cell>0, not if *at
                          // least one* have.
                          Assert (! ((interpolation_hp(new_fe_index,old_index).m() == 0)
                                     &&
                                     (interpolation_hp(new_fe_index,old_index).n() == 0)
                                     &&
                                     ((dof_handler->get_fe()[new_fe_index].dofs_per_cell > 0)
                                      &&
                                      (dof_handler->get_fe()[old_index].dofs_per_cell > 0))),
                                  ExcMessage ("The interpolation between two different "
                                              "elements you are trying to use here has "
                                              "not been implemented for this pair of "
                                              "elements!"));

                          // simple case where all children have the
                          // same FE index: just interpolate to their FE
                          // first and then use the standard routines
                          if (tmp.size() > 0)
                            interpolation_hp(new_fe_index,old_index).vmult (local_values, tmp);
                          else
                            local_values = 0;
                        }

                      if (cell->has_children() == false)
                        cell->set_dof_values (local_values, all_out[j]);
                      else
                        for (unsigned int child=0; child<cell->n_children(); ++child)
                          {
                            if (different_elements == true)
                              {
                                const unsigned int c_index =
                                  cell->child(child)->active_fe_index();
                                if (c_index != old_index)
                                  {
                                    AssertDimension (tmp.size(),
                                                     interpolation_hp(c_index,old_index).n());
                                    local_values.reinit(cell->child(child)->get_fe().dofs_per_cell, true);

                                    // do the interpolation. same problem as above
                                    Assert (! ((interpolation_hp(c_index,old_index).m() == 0)
                                               &&
                                               (interpolation_hp(c_index,old_index).n() == 0)
                                               &&
                                               ((dof_handler->get_fe()[c_index].dofs_per_cell > 0)
                                                &&
                                                (dof_handler->get_fe()[old_index].dofs_per_cell > 0))),
                                            ExcMessage ("The interpolation between two different "
                                                        "elements you are trying to use here has "
                                                        "not been implemented for this pair of "
                                                        "elements!"));

                                    if (tmp.size() > 0)
                                      interpolation_hp(c_index,old_index).vmult (local_values, tmp);
                                    else
                                      local_values = 0;
                                  }
                                else
                                  local_values = tmp;
                              }
                            tmp2.reinit (cell->child(child)->get_fe().dofs_per_cell, true);
                            cell->child(child)->get_fe().get_prolongation_matrix(child, cell->refinement_case())
                            .vmult (tmp2, local_values);
                            cell->child(child)->set_dof_values(tmp2,all_out[j]);
                          }
                    }
                  else
                    {
                      AssertDimension (dofs_per_cell, indexptr->size());
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        local_values(i)=all_in[j]((*indexptr)[i]);
                      cell->set_dof_values_by_interpolation(local_values,
                                                            all_out[j]);
                    }
                }
            }
          // children of cell were
          // deleted
          else if (valuesptr)
            {
              Assert (!cell->has_children(), ExcInternalError());
              Assert (indexptr == 0,
                      ExcInternalError());
              dofs.resize(dofs_per_cell);
              // get the local
              // indices
              cell->get_dof_indices(dofs);

              // distribute the
              // stored data to the
              // new vectors
              for (unsigned int j=0; j<size; ++j)
                {
                  // make sure that the size of
                  // the stored indices is the
                  // same as
                  // dofs_per_cell. this is
                  // kind of a test if we use
                  // the same fe in the hp
                  // case. to really do that
                  // test we would have to
                  // store the fe_index of all
                  // cells
                  const Vector<typename VECTOR::value_type> *data = 0;
                  const unsigned int active_fe_index = cell->active_fe_index();
                  if (active_fe_index != pointerstruct->second.active_fe_index)
                    {
                      const unsigned int old_index = pointerstruct->second.active_fe_index;
                      tmp.reinit (dofs_per_cell, true);
                      AssertDimension ((*valuesptr)[j].size(),
                                       interpolation_hp(active_fe_index,old_index).n());
                      AssertDimension (tmp.size(),
                                       interpolation_hp(active_fe_index,old_index).m());
                      interpolation_hp(active_fe_index,old_index).vmult (tmp, (*valuesptr)[j]);
                      data = &tmp;
                    }
                  else
                    data = &(*valuesptr)[j];


                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    all_out[j](dofs[i])=(*data)(i);
                }
            }
          // undefined status
          else
            Assert(false, ExcInternalError());
        }
    }
}



template<int dim, typename VECTOR, class DH>
void SolutionTransfer<dim, VECTOR, DH>::interpolate(const VECTOR &in,
                                                    VECTOR       &out) const
{
  Assert (in.size()==n_dofs_old,
          ExcDimensionMismatch(in.size(), n_dofs_old));
  Assert (out.size()==dof_handler->n_dofs(),
          ExcDimensionMismatch(out.size(), dof_handler->n_dofs()));

  std::vector<VECTOR> all_in(1);
  all_in[0] = in;
  std::vector<VECTOR> all_out(1);
  all_out[0] = out;
  interpolate(all_in,
              all_out);
  out=all_out[0];
}



template<int dim, typename VECTOR, class DH>
std::size_t
SolutionTransfer<dim, VECTOR, DH>::memory_consumption () const
{
  // at the moment we do not include the memory
  // consumption of the cell_map as we have no
  // real idea about memory consumption of a
  // std::map
  return (MemoryConsumption::memory_consumption (dof_handler) +
          MemoryConsumption::memory_consumption (n_dofs_old) +
          sizeof (prepared_for) +
          MemoryConsumption::memory_consumption (indices_on_cell) +
          MemoryConsumption::memory_consumption (dof_values_on_cell));
}



template<int dim, typename VECTOR, class DH>
std::size_t
SolutionTransfer<dim, VECTOR, DH>::Pointerstruct::memory_consumption () const
{
  return sizeof(*this);
}



/*-------------- Explicit Instantiations -------------------------------*/
#ifdef SOLUTION_TRANSFER_INSTANTIATE_PART_TWO
#define DIM_A 3
#define DIM_B 3
#else
#define DIM_A 1
#define DIM_B 2
#endif

// This file compiles the first quarter of the instantiations from solution_transfer.cc
// to reduce the compilation unit (and memory consumption)
#include "solution_transfer.inst"

DEAL_II_NAMESPACE_CLOSE
