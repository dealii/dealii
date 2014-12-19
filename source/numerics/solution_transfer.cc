// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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
      cell_map[std::make_pair(cell->level(),cell->index())]
        = Pointerstruct(&indices_on_cell[i], cell->active_fe_index());
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
          const unsigned int this_fe_index = pointerstruct->second.active_fe_index;
          const unsigned int dofs_per_cell=cell->get_dof_handler().get_fe()[this_fe_index].dofs_per_cell;
          local_values.reinit(dofs_per_cell, true);

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
          cell->set_dof_values_by_interpolation(local_values, out,
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
  template <typename DH>
  void extract_interpolation_matrices (const DH &,
                                       dealii::Table<2,FullMatrix<double> > &)
  {}

  template <int dim, int spacedim>
  void extract_interpolation_matrices (const dealii::hp::DoFHandler<dim,spacedim> &dof,
                                       dealii::Table<2,FullMatrix<double> > &matrices)
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
  Assert (prepared_for!=pure_refinement, ExcAlreadyPrepForRef());
  Assert (prepared_for!=coarsening_and_refinement,
          ExcAlreadyPrepForCoarseAndRef());

  const unsigned int in_size=all_in.size();
  Assert(in_size!=0, ExcNoInVectorsGiven());

  clear();

  const unsigned int n_active_cells = dof_handler->get_tria().n_active_cells();
  n_dofs_old = dof_handler->n_dofs();

  for (unsigned int i=0; i<in_size; ++i)
    {
      Assert(all_in[i].size()==n_dofs_old,
             ExcDimensionMismatch(all_in[i].size(),n_dofs_old));
    }

  // first count the number
  // of cells that will be coarsened
  // and that'll stay or be refined
  unsigned int n_cells_to_coarsen=0;
  unsigned int n_cells_to_stay_or_refine=0;
  for (typename DH::active_cell_iterator act_cell = dof_handler->begin_active();
       act_cell!=dof_handler->end(); ++act_cell)
    {
      if (act_cell->coarsen_flag_set())
        ++n_cells_to_coarsen;
      else
        ++n_cells_to_stay_or_refine;
    }
  Assert((n_cells_to_coarsen+n_cells_to_stay_or_refine)==n_active_cells,
         ExcInternalError());

  unsigned int n_coarsen_fathers=0;
  for (typename DH::cell_iterator cell=dof_handler->begin();
       cell!=dof_handler->end(); ++cell)
    if (!cell->active() && cell->child(0)->coarsen_flag_set())
      ++n_coarsen_fathers;
  Assert(n_cells_to_coarsen>=2*n_coarsen_fathers, ExcInternalError());

  // allocate the needed memory. initialize
  // the following arrays in an efficient
  // way, without copying much
  std::vector<std::vector<types::global_dof_index> >(n_cells_to_stay_or_refine)
  .swap(indices_on_cell);

  std::vector<std::vector<Vector<typename VECTOR::value_type> > >
  (n_coarsen_fathers,
   std::vector<Vector<typename VECTOR::value_type> > (in_size))
  .swap(dof_values_on_cell);

  Table<2,FullMatrix<double> > interpolation_hp;
  std::vector<std::vector<bool> > restriction_is_additive;

  internal::extract_interpolation_matrices (*dof_handler, interpolation_hp);
  internal::restriction_additive (dof_handler->get_fe(), restriction_is_additive);

  // we need counters for
  // the 'to_stay_or_refine' cells 'n_sr' and
  // the 'coarsen_fathers' cells 'n_cf',
  unsigned int n_sr=0, n_cf=0;
  for (typename DH::cell_iterator cell=dof_handler->begin();
       cell!=dof_handler->end(); ++cell)
    {
      // CASE 1: active cell that remains as it is
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
            = Pointerstruct(&indices_on_cell[n_sr], cell->active_fe_index());
          ++n_sr;
        }

      // CASE 2: cell is inactive but will become active
      else if (cell->has_children() && cell->child(0)->coarsen_flag_set())
        {
          // note that if one child has the
          // coarsen flag, then all should
          // have if Tria::prepare_* has
          // worked correctly
          for (unsigned int i=1; i<cell->n_children(); ++i)
            Assert(cell->child(i)->coarsen_flag_set(),
                   ExcMessage("It looks like you didn't call "
                              "Triangulation::prepare_coarsening_and_refinement before "
                              "calling the current function. This can't work."));

          // we will need to interpolate from the children of this cell
          // to the current one. in the hp context, this also means
          // we need to figure out which finite element space to interpolate
          // to since that is not implied by the global FE as in the non-hp
          // case.
          bool different_fe_on_children = false;
          for (unsigned int child=1; child<cell->n_children(); ++child)
            if (cell->child(child)->active_fe_index()
                != cell->child(0)->active_fe_index())
              {
                different_fe_on_children = true;
                break;
              }

          // take FE index from the child with most
          // degrees of freedom locally
          unsigned int most_general_child = 0;
          if (different_fe_on_children == true)
            for (unsigned int child=1; child<cell->n_children(); ++child)
              if (cell->child(child)->get_fe().dofs_per_cell >
                  cell->child(most_general_child)->get_fe().dofs_per_cell)
                most_general_child = child;
          const unsigned int target_fe_index = cell->child(most_general_child)->active_fe_index();

          const unsigned int dofs_per_cell=cell->get_dof_handler().get_fe()[target_fe_index].dofs_per_cell;

          std::vector<Vector<typename VECTOR::value_type> >(in_size,
                                                            Vector<typename VECTOR::value_type>(dofs_per_cell))
          .swap(dof_values_on_cell[n_cf]);


          // store the data of each of the input vectors. get this data
          // as interpolated onto a finite element space that encompasses
          // that of all the children. note that cell->get_interpolated_dof_values
          // already does all of the interpolations between spaces
          for (unsigned int j=0; j<in_size; ++j)
            cell->get_interpolated_dof_values(all_in[j],
                                              dof_values_on_cell[n_cf][j],
                                              target_fe_index);
          cell_map[std::make_pair(cell->level(), cell->index())]
            = Pointerstruct(&dof_values_on_cell[n_cf], target_fe_index);
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

          // cell stayed as it was or was refined
          if (indexptr)
            {
              Assert (valuesptr == 0,
                      ExcInternalError());

              const unsigned int old_fe_index =
                pointerstruct->second.active_fe_index;

              // get the values of
              // each of the input
              // data vectors on this
              // cell and prolong it
              // to its children
              unsigned int in_size = indexptr->size();
              for (unsigned int j=0; j<size; ++j)
                {
                  tmp.reinit (in_size, true);
                  for (unsigned int i=0; i<in_size; ++i)
                    tmp(i) = all_in[j]((*indexptr)[i]);

                  cell->set_dof_values_by_interpolation (tmp, all_out[j],
                                                         old_fe_index);
                }
            }
          else if (valuesptr)
            // the children of this cell were
            // deleted
            {
              Assert (!cell->has_children(), ExcInternalError());
              Assert (indexptr == 0,
                      ExcInternalError());

              const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
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
