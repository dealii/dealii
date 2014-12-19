// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/data_out_rotation.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


//TODO: Update documentation
//TODO: Unify code for dimensions


//TODO: build_some_patches isn't going to work if first_cell/next_cell
//don't iterate over all cells and if cell data is requested. in that
//case, we need to calculate cell_number as in the DataOut class

// Not implemented for 3D


namespace internal
{
  namespace DataOutRotation
  {
    template <int dim, int spacedim>
    ParallelData<dim,spacedim>::
    ParallelData (const unsigned int n_datasets,
                  const unsigned int n_subdivisions,
                  const unsigned int n_patches_per_circle,
                  const std::vector<unsigned int> &n_postprocessor_outputs,
                  const Mapping<dim,spacedim> &mapping,
                  const std::vector<std_cxx11::shared_ptr<dealii::hp::FECollection<dim,spacedim> > > &finite_elements,
                  const UpdateFlags update_flags)
      :
      internal::DataOut::
      ParallelDataBase<dim,spacedim> (n_datasets,
                                      n_subdivisions,
                                      n_postprocessor_outputs,
                                      mapping,
                                      finite_elements,
                                      update_flags,
                                      false),
      n_patches_per_circle (n_patches_per_circle)
    {}



    /**
     * In a WorkStream context, use this function to append the patch computed
     * by the parallel stage to the array of patches.
     */
    template <int dim, int spacedim>
    void
    append_patch_to_list (const std::vector<DataOutBase::Patch<dim+1,spacedim+1> > &new_patches,
                          std::vector<DataOutBase::Patch<dim+1,spacedim+1> > &patches)
    {
      for (unsigned int i=0; i<new_patches.size(); ++i)
        {
          patches.push_back (new_patches[i]);
          patches.back().patch_index = patches.size()-1;
        }
    }
  }
}



template <int dim, class DH>
void
DataOutRotation<dim,DH>::
build_one_patch (const cell_iterator *cell,
                 internal::DataOutRotation::ParallelData<dimension, space_dimension> &data,
                 std::vector<DataOutBase::Patch<dimension+1,space_dimension+1> > &patches)
{
  if (dim == 3)
    {
      // would this function make any sense after all? who would want to
      // output/compute in four space dimensions?
      Assert (false, ExcNotImplemented());
      return;
    }

  Assert ((*cell)->is_locally_owned(),
          ExcNotImplemented());

  const unsigned int n_patches_per_circle = data.n_patches_per_circle;

  // another abbreviation denoting the number of q_points in each direction
  const unsigned int n_points = data.n_subdivisions+1;

  // set up an array that holds the directions in the plane of rotation in
  // which we will put points in the whole domain (not the rotationally
  // reduced one in which the computation took place. for simplicity add the
  // initial direction at the end again
  std::vector<Point<dimension+1> > angle_directions (n_patches_per_circle+1);
  for (unsigned int i=0; i<=n_patches_per_circle; ++i)
    {
      angle_directions[i][dimension-1] = std::cos(2*numbers::PI *
                                                  i/n_patches_per_circle);
      angle_directions[i][dimension] = std::sin(2*numbers::PI *
                                                i/n_patches_per_circle);
    }

  for (unsigned int angle=0; angle<n_patches_per_circle; ++angle)
    {
      // first compute the vertices of the patch. note that they will have to
      // be computed from the vertices of the cell, which has one dimension
      // less, however.
      switch (dimension)
        {
        case 1:
        {
          const double r1 = (*cell)->vertex(0)(0),
                       r2 = (*cell)->vertex(1)(0);
          Assert (r1 >= 0, ExcRadialVariableHasNegativeValues(r1));
          Assert (r2 >= 0, ExcRadialVariableHasNegativeValues(r2));

          patches[angle].vertices[0] = r1*angle_directions[angle];
          patches[angle].vertices[1] = r2*angle_directions[angle];
          patches[angle].vertices[2] = r1*angle_directions[angle+1];
          patches[angle].vertices[3] = r2*angle_directions[angle+1];

          break;
        };

        case 2:
        {
          for (unsigned int vertex=0;
               vertex<GeometryInfo<dimension>::vertices_per_cell;
               ++vertex)
            {
              const Point<dimension> v = (*cell)->vertex(vertex);

              // make sure that the radial variable does attain negative
              // values
              Assert (v(0) >= 0, ExcRadialVariableHasNegativeValues(v(0)));

              // now set the vertices of the patch
              patches[angle].vertices[vertex] = v(0) * angle_directions[angle];
              patches[angle].vertices[vertex][0] = v(1);

              patches[angle].vertices[vertex+GeometryInfo<dimension>::vertices_per_cell]
                = v(0) * angle_directions[angle+1];
              patches[angle].vertices[vertex+GeometryInfo<dimension>::vertices_per_cell][0]
                = v(1);
            };

          break;
        };

        default:
          Assert (false, ExcNotImplemented());
        };

      unsigned int offset=0;

      // then fill in data
      if (data.n_datasets > 0)
        {
          data.reinit_all_fe_values(this->dof_data, *cell);
          // first fill dof_data
          for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
            {
              const FEValuesBase<dimension> &fe_patch_values
                = data.get_present_fe_values(dataset);
              const unsigned int n_components
                = fe_patch_values.get_fe().n_components();
              const DataPostprocessor<dim> *postprocessor=this->dof_data[dataset]->postprocessor;
              if (postprocessor != 0)
                {
                  // we have to postprocess the
                  // data, so determine, which
                  // fields have to be updated
                  const UpdateFlags update_flags=postprocessor->get_needed_update_flags();

                  if (n_components == 1)
                    {
                      // at each point there is
                      // only one component of
                      // value, gradient etc.
                      if (update_flags & update_values)
                        this->dof_data[dataset]->get_function_values (fe_patch_values,
                                                                      data.patch_values);
                      if (update_flags & update_gradients)
                        this->dof_data[dataset]->get_function_gradients (fe_patch_values,
                                                                         data.patch_gradients);
                      if (update_flags & update_hessians)
                        this->dof_data[dataset]->get_function_hessians (fe_patch_values,
                                                                        data.patch_hessians);

                      if (update_flags & update_quadrature_points)
                        data.patch_evaluation_points = fe_patch_values.get_quadrature_points();

                      std::vector<Point<space_dimension> > dummy_normals;
                      postprocessor->
                      compute_derived_quantities_scalar(data.patch_values,
                                                        data.patch_gradients,
                                                        data.patch_hessians,
                                                        dummy_normals,
                                                        data.patch_evaluation_points,
                                                        data.postprocessed_values[dataset]);
                    }
                  else
                    {
                      data.resize_system_vectors(n_components);

                      // at each point there is a vector valued function and
                      // its derivative...
                      if (update_flags & update_values)
                        this->dof_data[dataset]->get_function_values (fe_patch_values,
                                                                      data.patch_values_system);
                      if (update_flags & update_gradients)
                        this->dof_data[dataset]->get_function_gradients (fe_patch_values,
                                                                         data.patch_gradients_system);
                      if (update_flags & update_hessians)
                        this->dof_data[dataset]->get_function_hessians (fe_patch_values,
                                                                        data.patch_hessians_system);

                      if (update_flags & update_quadrature_points)
                        data.patch_evaluation_points = fe_patch_values.get_quadrature_points();

                      std::vector<Point<space_dimension> > dummy_normals;
                      postprocessor->
                      compute_derived_quantities_vector(data.patch_values_system,
                                                        data.patch_gradients_system,
                                                        data.patch_hessians_system,
                                                        dummy_normals,
                                                        data.patch_evaluation_points,
                                                        data.postprocessed_values[dataset]);
                    }

                  for (unsigned int component=0;
                       component<this->dof_data[dataset]->n_output_variables;
                       ++component)
                    {
                      switch (dimension)
                        {
                        case 1:
                          for (unsigned int x=0; x<n_points; ++x)
                            for (unsigned int y=0; y<n_points; ++y)
                              patches[angle].data(offset+component,
                                                  x*n_points + y)
                                = data.postprocessed_values[dataset][x](component);
                          break;

                        case 2:
                          for (unsigned int x=0; x<n_points; ++x)
                            for (unsigned int y=0; y<n_points; ++y)
                              for (unsigned int z=0; z<n_points; ++z)
                                patches[angle].data(offset+component,
                                                    x*n_points*n_points +
                                                    y*n_points +
                                                    z)
                                  = data.postprocessed_values[dataset][x*n_points+z](component);
                          break;

                        default:
                          Assert (false, ExcNotImplemented());
                        }
                    }
                }
              else if (n_components == 1)
                {
                  this->dof_data[dataset]->get_function_values (fe_patch_values,
                                                                data.patch_values);

                  switch (dimension)
                    {
                    case 1:
                      for (unsigned int x=0; x<n_points; ++x)
                        for (unsigned int y=0; y<n_points; ++y)
                          patches[angle].data(offset,
                                              x*n_points + y)
                            = data.patch_values[x];
                      break;

                    case 2:
                      for (unsigned int x=0; x<n_points; ++x)
                        for (unsigned int y=0; y<n_points; ++y)
                          for (unsigned int z=0; z<n_points; ++z)
                            patches[angle].data(offset,
                                                x*n_points*n_points +
                                                y +
                                                z*n_points)
                              = data.patch_values[x*n_points+z];
                      break;

                    default:
                      Assert (false, ExcNotImplemented());
                    }
                }
              else
                // system of components
                {
                  data.resize_system_vectors(n_components);
                  this->dof_data[dataset]->get_function_values (fe_patch_values,
                                                                data.patch_values_system);

                  for (unsigned int component=0; component<n_components;
                       ++component)
                    {
                      switch (dimension)
                        {
                        case 1:
                          for (unsigned int x=0; x<n_points; ++x)
                            for (unsigned int y=0; y<n_points; ++y)
                              patches[angle].data(offset+component,
                                                  x*n_points + y)
                                = data.patch_values_system[x](component);
                          break;

                        case 2:
                          for (unsigned int x=0; x<n_points; ++x)
                            for (unsigned int y=0; y<n_points; ++y)
                              for (unsigned int z=0; z<n_points; ++z)
                                patches[angle].data(offset+component,
                                                    x*n_points*n_points +
                                                    y*n_points +
                                                    z)
                                  = data.patch_values_system[x*n_points+z](component);
                          break;

                        default:
                          Assert (false, ExcNotImplemented());
                        }
                    }
                }
              offset+=this->dof_data[dataset]->n_output_variables;
            }

          // then do the cell data
          for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
            {
              // we need to get at the number of the cell to which this face
              // belongs in order to access the cell data. this is not readily
              // available, so choose the following rather inefficient way:
              Assert ((*cell)->active(),
                      ExcMessage("Cell must be active for cell data"));
              const unsigned int cell_number
                = std::distance (this->triangulation->begin_active(),
                                 typename Triangulation<dimension,space_dimension>::active_cell_iterator(*cell));
              const double value
                = this->cell_data[dataset]->get_cell_data_value (cell_number);
              switch (dimension)
                {
                case 1:
                  for (unsigned int x=0; x<n_points; ++x)
                    for (unsigned int y=0; y<n_points; ++y)
                      patches[angle].data(dataset+offset,
                                          x*n_points +
                                          y)
                        = value;
                  break;

                case 2:
                  for (unsigned int x=0; x<n_points; ++x)
                    for (unsigned int y=0; y<n_points; ++y)
                      for (unsigned int z=0; z<n_points; ++z)
                        patches[angle].data(dataset+offset,
                                            x*n_points*n_points +
                                            y*n_points +
                                            z)
                          = value;
                  break;

                default:
                  Assert (false, ExcNotImplemented());
                }
            }
        }
    }
}



template <int dim, class DH>
void DataOutRotation<dim,DH>::build_patches (const unsigned int n_patches_per_circle,
                                             const unsigned int nnnn_subdivisions)
{
  // Check consistency of redundant
  // template parameter
  Assert (dim==dimension, ExcDimensionMismatch(dim, dimension));
  typedef DataOut_DoFData<DH,dimension+1> BaseClass;
  Assert (this->triangulation != 0,
          typename BaseClass::ExcNoTriangulationSelected());

  const unsigned int n_subdivisions = (nnnn_subdivisions != 0)
                                      ? nnnn_subdivisions
                                      : this->default_subdivisions;
  Assert (n_subdivisions >= 1,
          ExcInvalidNumberOfSubdivisions(n_subdivisions));

  unsigned int n_datasets=this->cell_data.size();
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    n_datasets+= this->dof_data[i]->n_output_variables;

  UpdateFlags update_flags=update_values | update_quadrature_points;
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    if (this->dof_data[i]->postprocessor)
      update_flags |= this->dof_data[i]->postprocessor->get_needed_update_flags();
  // perhaps update_normal_vectors is present,
  // which would only be useful on faces, but
  // we may not use it here.
  Assert (!(update_flags & update_normal_vectors),
          ExcMessage("The update of normal vectors may not be requested for "
                     "evaluation of data on cells via DataPostprocessor."));

  // first count the cells we want to
  // create patches of and make sure
  // there is enough memory for that
  std::vector<cell_iterator> all_cells;
  for (cell_iterator cell=first_cell(); cell != this->triangulation->end();
       cell = next_cell(cell))
    all_cells.push_back (cell);

  // then also take into account that
  // we want more than one patch to
  // come out of every cell, as they
  // are repeated around the axis of
  // rotation
  this->patches.clear();
  this->patches.reserve (all_cells.size() * n_patches_per_circle);


  std::vector<unsigned int> n_postprocessor_outputs (this->dof_data.size());
  for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
    if (this->dof_data[dataset]->postprocessor)
      n_postprocessor_outputs[dataset] = this->dof_data[dataset]->n_output_variables;
    else
      n_postprocessor_outputs[dataset] = 0;

  const MappingQ1<dimension, space_dimension> mapping;
  internal::DataOutRotation::ParallelData<dimension, space_dimension>
  thread_data (n_datasets,
               n_subdivisions, n_patches_per_circle,
               n_postprocessor_outputs, mapping, this->get_finite_elements(),
               update_flags);
  std::vector<DataOutBase::Patch<dimension+1,space_dimension+1> >
  new_patches (n_patches_per_circle);
  for (unsigned int i=0; i<new_patches.size(); ++i)
    {
      new_patches[i].n_subdivisions = n_subdivisions;
      new_patches[i].data.reinit (n_datasets,
                                  Utilities::fixed_power<dimension+1>(n_subdivisions+1));
    }

  // now build the patches in parallel
  WorkStream::run (&all_cells[0],
                   &all_cells[0]+all_cells.size(),
                   std_cxx11::bind(&DataOutRotation<dim,DH>::build_one_patch,
                                   this, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3),
                   std_cxx11::bind(&internal::DataOutRotation
                                   ::append_patch_to_list<dim,space_dimension>,
                                   std_cxx11::_1, std_cxx11::ref(this->patches)),
                   thread_data,
                   new_patches);
}



template <int dim, class DH>
typename DataOutRotation<dim,DH>::cell_iterator
DataOutRotation<dim,DH>::first_cell ()
{
  return this->triangulation->begin_active ();
}


template <int dim, class DH>
typename DataOutRotation<dim,DH>::cell_iterator
DataOutRotation<dim,DH>::next_cell (const cell_iterator &cell)
{
  // convert the iterator to an
  // active_iterator and advance
  // this to the next active cell
  typename Triangulation<dimension,space_dimension>::active_cell_iterator active_cell = cell;
  ++active_cell;
  return active_cell;
}



// explicit instantiations
#include "data_out_rotation.inst"


DEAL_II_NAMESPACE_CLOSE
