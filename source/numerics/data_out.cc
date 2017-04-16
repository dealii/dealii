// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2016 by the deal.II authors
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

#include <deal.II/base/work_stream.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <sstream>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOut
  {
    template <int dim, int spacedim>
    ParallelData<dim,spacedim>::
    ParallelData (const unsigned int n_datasets,
                  const unsigned int n_subdivisions,
                  const std::vector<unsigned int> &n_postprocessor_outputs,
                  const Mapping<dim,spacedim> &mapping,
                  const std::vector<std::shared_ptr<dealii::hp::FECollection<dim,spacedim> > > &finite_elements,
                  const UpdateFlags update_flags,
                  const std::vector<std::vector<unsigned int> > &cell_to_patch_index_map)
      :
      ParallelDataBase<dim,spacedim> (n_datasets,
                                      n_subdivisions,
                                      n_postprocessor_outputs,
                                      mapping,
                                      finite_elements,
                                      update_flags,
                                      false),
      cell_to_patch_index_map (&cell_to_patch_index_map)
    {}
  }
}



template <int dim, typename DoFHandlerType>
void
DataOut<dim,DoFHandlerType>::
build_one_patch
(const std::pair<cell_iterator, unsigned int>                                                *cell_and_index,
 internal::DataOut::ParallelData<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &scratch_data,
 const unsigned int                                                                           n_subdivisions,
 const CurvedCellRegion                                                                       curved_cell_region)
{
  // first create the output object that we will write into
  ::dealii::DataOutBase::Patch<DoFHandlerType::dimension, DoFHandlerType::space_dimension> patch;
  patch.n_subdivisions = n_subdivisions;

  // set the vertices of the patch. if the mapping does not preserve locations
  // (e.g. MappingQEulerian), we need to compute the offset of the vertex for
  // the graphical output. Otherwise, we can just use the vertex info.
  for (unsigned int vertex=0; vertex<GeometryInfo<DoFHandlerType::dimension>::vertices_per_cell; ++vertex)
    if (scratch_data.mapping_collection[0].preserves_vertex_locations())
      patch.vertices[vertex] = cell_and_index->first->vertex(vertex);
    else
      patch.vertices[vertex] = scratch_data.mapping_collection[0].transform_unit_to_real_cell
                               (cell_and_index->first,
                                GeometryInfo<DoFHandlerType::dimension>::unit_cell_vertex (vertex));

  // create DoFHandlerType::active_cell_iterator and initialize FEValues
  scratch_data.reinit_all_fe_values(this->dof_data, cell_and_index->first);

  const FEValuesBase<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &fe_patch_values
    = scratch_data.get_present_fe_values (0);

  const unsigned int n_q_points = fe_patch_values.n_quadrature_points;

  // depending on the requested output of curved cells, if necessary
  // append the quadrature points to the last rows of the patch.data
  // member. This is the case if we want to produce curved cells at the
  // boundary and this cell actually is at the boundary, or else if we
  // want to produce curved cells everywhere
  //
  // note: a cell is *always* at the boundary if dim<spacedim
  if (curved_cell_region==curved_inner_cells
      ||
      (curved_cell_region==curved_boundary
       &&
       (cell_and_index->first->at_boundary()
        ||
        (DoFHandlerType::dimension != DoFHandlerType::space_dimension))))
    {
      Assert(patch.space_dim==DoFHandlerType::space_dimension, ExcInternalError());

      // set the flag indicating that for this cell the points are
      // explicitly given
      patch.points_are_available = true;

      // then resize the patch.data member in order to have enough memory for
      // the quadrature points as well, and copy the quadrature points there
      const std::vector<Point<DoFHandlerType::space_dimension> > &q_points
        = fe_patch_values.get_quadrature_points();

      patch.data.reinit (scratch_data.n_datasets+DoFHandlerType::space_dimension, n_q_points);
      for (unsigned int i=0; i<DoFHandlerType::space_dimension; ++i)
        for (unsigned int q=0; q<n_q_points; ++q)
          patch.data(patch.data.size(0)-DoFHandlerType::space_dimension+i,q) = q_points[q][i];
    }
  else
    {
      patch.data.reinit(scratch_data.n_datasets, n_q_points);
      patch.points_are_available = false;
    }


  if (scratch_data.n_datasets > 0)
    {
      // counter for data records
      unsigned int offset=0;

      // first fill dof_data
      for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
        {
          const FEValuesBase<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &this_fe_patch_values
            = scratch_data.get_present_fe_values (dataset);
          const unsigned int n_components =
            this_fe_patch_values.get_fe().n_components();

          const DataPostprocessor<DoFHandlerType::space_dimension> *postprocessor
            = this->dof_data[dataset]->postprocessor;

          if (postprocessor != nullptr)
            {
              // we have to postprocess the data, so determine, which fields
              // have to be updated
              const UpdateFlags update_flags=postprocessor->get_needed_update_flags();
              if (n_components == 1)
                {
                  // at each point there is only one component of value,
                  // gradient etc.
                  if (update_flags & update_values)
                    this->dof_data[dataset]->get_function_values (this_fe_patch_values,
                                                                  scratch_data.patch_values_scalar.solution_values);
                  if (update_flags & update_gradients)
                    this->dof_data[dataset]->get_function_gradients (this_fe_patch_values,
                                                                     scratch_data.patch_values_scalar.solution_gradients);
                  if (update_flags & update_hessians)
                    this->dof_data[dataset]->get_function_hessians (this_fe_patch_values,
                                                                    scratch_data.patch_values_scalar.solution_hessians);

                  if (update_flags & update_quadrature_points)
                    scratch_data.patch_values_scalar.evaluation_points = this_fe_patch_values.get_quadrature_points();

                  const typename DoFHandlerType::active_cell_iterator dh_cell(&cell_and_index->first->get_triangulation(),
                                                                              cell_and_index->first->level(),
                                                                              cell_and_index->first->index(),
                                                                              this->dof_data[dataset]->dof_handler);
                  scratch_data.patch_values_scalar.template set_cell<DoFHandlerType> (dh_cell);

                  postprocessor->
                  evaluate_scalar_field(scratch_data.patch_values_scalar,
                                        scratch_data.postprocessed_values[dataset]);
                }
              else
                {
                  scratch_data.resize_system_vectors (n_components);

                  // at each point there is a vector valued function and its
                  // derivative...
                  if (update_flags & update_values)
                    this->dof_data[dataset]->get_function_values (this_fe_patch_values,
                                                                  scratch_data.patch_values_system.solution_values);
                  if (update_flags & update_gradients)
                    this->dof_data[dataset]->get_function_gradients (this_fe_patch_values,
                                                                     scratch_data.patch_values_system.solution_gradients);
                  if (update_flags & update_hessians)
                    this->dof_data[dataset]->get_function_hessians (this_fe_patch_values,
                                                                    scratch_data.patch_values_system.solution_hessians);

                  if (update_flags & update_quadrature_points)
                    scratch_data.patch_values_system.evaluation_points = this_fe_patch_values.get_quadrature_points();

                  const typename DoFHandlerType::active_cell_iterator dh_cell(&cell_and_index->first->get_triangulation(),
                                                                              cell_and_index->first->level(),
                                                                              cell_and_index->first->index(),
                                                                              this->dof_data[dataset]->dof_handler);
                  scratch_data.patch_values_system.template set_cell<DoFHandlerType> (dh_cell);

                  postprocessor->
                  evaluate_vector_field(scratch_data.patch_values_system,
                                        scratch_data.postprocessed_values[dataset]);
                }

              for (unsigned int q=0; q<n_q_points; ++q)
                for (unsigned int component=0;
                     component<this->dof_data[dataset]->n_output_variables;
                     ++component)
                  patch.data(offset+component,q)
                    = scratch_data.postprocessed_values[dataset][q](component);
            }
          else
            // now we use the given data vector without modifications. again,
            // we treat single component functions separately for efficiency
            // reasons.
            if (n_components == 1)
              {
                this->dof_data[dataset]->get_function_values (this_fe_patch_values,
                                                              scratch_data.patch_values_scalar.solution_values);
                for (unsigned int q=0; q<n_q_points; ++q)
                  patch.data(offset,q) = scratch_data.patch_values_scalar.solution_values[q];
              }
            else
              {
                scratch_data.resize_system_vectors(n_components);
                this->dof_data[dataset]->get_function_values (this_fe_patch_values,
                                                              scratch_data.patch_values_system.solution_values);
                for (unsigned int component=0; component<n_components;
                     ++component)
                  for (unsigned int q=0; q<n_q_points; ++q)
                    patch.data(offset+component,q) =
                      scratch_data.patch_values_system.solution_values[q](component);
              }
          // increment the counter for the actual data record
          offset+=this->dof_data[dataset]->n_output_variables;
        }

      // then do the cell data. only compute the number of a cell if needed;
      // also make sure that we only access cell data if the
      // first_cell/next_cell functions only return active cells
      if (this->cell_data.size() != 0)
        {
          Assert (!cell_and_index->first->has_children(), ExcNotImplemented());

          for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
            {
              const double value
                = this->cell_data[dataset]->get_cell_data_value (cell_and_index->second);
              for (unsigned int q=0; q<n_q_points; ++q)
                patch.data(offset+dataset,q) = value;
            }
        }
    }


  for (unsigned int f=0; f<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++f)
    {
      // let's look up whether the neighbor behind that face is noted in the
      // table of cells which we treat. this can only happen if the neighbor
      // exists, and is on the same level as this cell, but it may also happen
      // that the neighbor is not a member of the range of cells over which we
      // loop, in which case the respective entry in the
      // cell_to_patch_index_map will have the value no_neighbor. (note that
      // since we allocated only as much space in this array as the maximum
      // index of the cells we loop over, not every neighbor may have its
      // space in it, so we have to assume that it is extended by values
      // no_neighbor)
      if (cell_and_index->first->at_boundary(f)
          ||
          (cell_and_index->first->neighbor(f)->level() != cell_and_index->first->level()))
        {
          patch.neighbors[f] = numbers::invalid_unsigned_int;
          continue;
        }

      const cell_iterator neighbor = cell_and_index->first->neighbor(f);
      Assert (static_cast<unsigned int>(neighbor->level()) <
              scratch_data.cell_to_patch_index_map->size(),
              ExcInternalError());
      if ((static_cast<unsigned int>(neighbor->index()) >=
           (*scratch_data.cell_to_patch_index_map)[neighbor->level()].size())
          ||
          ((*scratch_data.cell_to_patch_index_map)[neighbor->level()][neighbor->index()]
           ==
           dealii::DataOutBase::Patch<DoFHandlerType::dimension>::no_neighbor))
        {
          patch.neighbors[f] = numbers::invalid_unsigned_int;
          continue;
        }

      // now, there is a neighbor, so get its patch number and set it for the
      // neighbor index
      patch.neighbors[f]
        = (*scratch_data.cell_to_patch_index_map)[neighbor->level()][neighbor->index()];
    }

  const unsigned int patch_idx =
    (*scratch_data.cell_to_patch_index_map)[cell_and_index->first->level()][cell_and_index->first->index()];
  // did we mess up the indices?
  Assert(patch_idx < this->patches.size(), ExcInternalError());
  patch.patch_index = patch_idx;

  // Put the patch into the patches vector. instead of copying the data,
  // simply swap the contents to avoid the penalty of writing into another
  // processor's memory
  this->patches[patch_idx].swap (patch);
}



template <int dim, typename DoFHandlerType>
void DataOut<dim,DoFHandlerType>::build_patches (const unsigned int n_subdivisions)
{
  build_patches (StaticMappingQ1<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::mapping,
                 n_subdivisions, no_curved_cells);
}



template <int dim, typename DoFHandlerType>
void DataOut<dim,DoFHandlerType>::build_patches
(const Mapping<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &mapping,
 const unsigned int                                                        n_subdivisions_,
 const CurvedCellRegion                                                    curved_region)
{
  // Check consistency of redundant template parameter
  Assert (dim==DoFHandlerType::dimension, ExcDimensionMismatch(dim, DoFHandlerType::dimension));

  Assert (this->triangulation != nullptr,
          Exceptions::DataOut::ExcNoTriangulationSelected());

  const unsigned int n_subdivisions = (n_subdivisions_ != 0)
                                      ? n_subdivisions_
                                      : this->default_subdivisions;
  Assert (n_subdivisions >= 1,
          Exceptions::DataOut::ExcInvalidNumberOfSubdivisions(n_subdivisions));

  this->validate_dataset_names();

  // First count the cells we want to create patches of. Also fill the object
  // that maps the cell indices to the patch numbers, as this will be needed
  // for generation of neighborship information.
  // Note, there is a confusing mess of different indices here at play:
  // patch_index - the index of a patch in all_cells
  // cell->index - only unique on each level, used in cell_to_patch_index_map
  // active_index - index for a cell when counting from begin_active() using ++cell
  // cell_index - unique index of a cell counted using next_locally_owned_cell()
  //              starting from first_locally_owned_cell()
  //
  // It turns out that we create one patch for each selected cell, so patch_index==cell_index.
  //
  // will be cell_to_patch_index_map[cell->level][cell->index] = patch_index
  std::vector<std::vector<unsigned int> > cell_to_patch_index_map;
  cell_to_patch_index_map.resize (this->triangulation->n_levels());
  for (unsigned int l=0; l<this->triangulation->n_levels(); ++l)
    {
      // max_index is the largest cell->index on level l
      unsigned int max_index = 0;
      for (cell_iterator cell=first_locally_owned_cell(); cell != this->triangulation->end();
           cell = next_locally_owned_cell(cell))
        if (static_cast<unsigned int>(cell->level()) == l)
          max_index = std::max (max_index,
                                static_cast<unsigned int>(cell->index()));

      cell_to_patch_index_map[l].resize (max_index+1,
                                         dealii::DataOutBase::Patch<DoFHandlerType::dimension,
                                         DoFHandlerType::space_dimension>::no_neighbor);
    }

  // will be all_cells[patch_index] = pair(cell, active_index)
  std::vector<std::pair<cell_iterator, unsigned int> > all_cells;
  {
    // important: we need to compute the active_index of the cell in the range
    // 0..n_active_cells() because this is where we need to look up cell
    // data from (cell data vectors do not have the length distance computed by
    // first_locally_owned_cell/next_locally_owned_cell because this might skip
    // some values (FilteredIterator).
    active_cell_iterator active_cell = this->triangulation->begin_active();
    unsigned int active_index = 0;
    cell_iterator cell = first_locally_owned_cell();
    for (; cell != this->triangulation->end();
         cell = next_locally_owned_cell(cell))
      {
        // move forward until active_cell points at the cell (cell) we are looking
        // at to compute the current active_index
        while (active_cell!=this->triangulation->end()
               && cell->active()
               && active_cell_iterator(cell) != active_cell)
          {
            ++active_cell;
            ++active_index;
          }

        Assert (static_cast<unsigned int>(cell->level()) <
                cell_to_patch_index_map.size(),
                ExcInternalError());
        Assert (static_cast<unsigned int>(cell->index()) <
                cell_to_patch_index_map[cell->level()].size(),
                ExcInternalError());
        Assert (active_index < this->triangulation->n_active_cells(),
                ExcInternalError());
        cell_to_patch_index_map[cell->level()][cell->index()] = all_cells.size();

        all_cells.emplace_back (cell, active_index);
      }
  }

  this->patches.clear ();
  this->patches.resize(all_cells.size());

  // now create a default object for the WorkStream object to work with
  unsigned int n_datasets=this->cell_data.size();
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    n_datasets += this->dof_data[i]->n_output_variables;

  std::vector<unsigned int> n_postprocessor_outputs (this->dof_data.size());
  for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
    if (this->dof_data[dataset]->postprocessor)
      n_postprocessor_outputs[dataset] = this->dof_data[dataset]->n_output_variables;
    else
      n_postprocessor_outputs[dataset] = 0;

  const CurvedCellRegion curved_cell_region
    = (n_subdivisions<2 ? no_curved_cells : curved_region);

  UpdateFlags update_flags = update_values;
  if (curved_cell_region != no_curved_cells)
    update_flags |= update_quadrature_points;

  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    if (this->dof_data[i]->postprocessor)
      update_flags |= this->dof_data[i]->postprocessor->get_needed_update_flags();
  // perhaps update_normal_vectors is present, which would only be useful on
  // faces, but we may not use it here.
  Assert (!(update_flags & update_normal_vectors),
          ExcMessage("The update of normal vectors may not be requested for evaluation of "
                     "data on cells via DataPostprocessor."));

  internal::DataOut::ParallelData<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
  thread_data (n_datasets, n_subdivisions,
               n_postprocessor_outputs,
               mapping,
               this->get_finite_elements(),
               update_flags,
               cell_to_patch_index_map);

  // now build the patches in parallel
  if (all_cells.size() > 0)
    WorkStream::run (&all_cells[0],
                     &all_cells[0]+all_cells.size(),
                     std::bind(&DataOut<dim,DoFHandlerType>::build_one_patch,
                               this,
                               std::placeholders::_1,
                               std::placeholders::_2,
                               /* no std::placeholders::_3, since this function doesn't actually need a
                                  copy data object -- it just writes everything right into the
                                  output array */
                               n_subdivisions,
                               curved_cell_region),
                     // no copy-local-to-global function needed here
                     std::function<void (const int &)>(),
                     thread_data,
                     /* dummy CopyData object = */ 0,
                     // experimenting shows that we can make things run a bit
                     // faster if we increase the number of cells we work on
                     // per item (i.e., WorkStream's chunk_size argument,
                     // about 10% improvement) and the items in flight at any
                     // given time (another 5% on the testcase discussed in
                     // @ref workstream_paper, on 32 cores) and if
                     8*MultithreadInfo::n_threads(),
                     64);
}



template <int dim, typename DoFHandlerType>
typename DataOut<dim,DoFHandlerType>::cell_iterator
DataOut<dim,DoFHandlerType>::first_cell ()
{
  return this->triangulation->begin_active ();
}



template <int dim, typename DoFHandlerType>
typename DataOut<dim,DoFHandlerType>::cell_iterator
DataOut<dim,DoFHandlerType>::next_cell
(const typename DataOut<dim,DoFHandlerType>::cell_iterator &cell)
{
  // convert the iterator to an active_iterator and advance this to the next
  // active cell
  typename Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::
  active_cell_iterator active_cell = cell;
  ++active_cell;
  return active_cell;
}



template <int dim, typename DoFHandlerType>
typename DataOut<dim,DoFHandlerType>::cell_iterator
DataOut<dim,DoFHandlerType>::first_locally_owned_cell ()
{
  typename DataOut<dim,DoFHandlerType>::cell_iterator
  cell = first_cell();

  // skip cells if the current one has no children (is active) and is a ghost
  // or artificial cell
  while ((cell != this->triangulation->end()) &&
         (cell->has_children() == false) &&
         !cell->is_locally_owned())
    cell = next_cell(cell);

  return cell;
}



template <int dim, typename DoFHandlerType>
typename DataOut<dim,DoFHandlerType>::cell_iterator
DataOut<dim,DoFHandlerType>::next_locally_owned_cell
(const typename DataOut<dim,DoFHandlerType>::cell_iterator &old_cell)
{
  typename DataOut<dim,DoFHandlerType>::cell_iterator
  cell = next_cell(old_cell);
  while ((cell != this->triangulation->end()) &&
         (cell->has_children() == false) &&
         !cell->is_locally_owned())
    cell = next_cell(cell);
  return cell;
}


// explicit instantiations
#include "data_out.inst"

DEAL_II_NAMESPACE_CLOSE
