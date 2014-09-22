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
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutFaces
  {
    template <int dim, int spacedim>
    ParallelData<dim,spacedim>::
    ParallelData (const unsigned int n_datasets,
                  const unsigned int n_subdivisions,
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
                                      true)
    {}



    /**
     * In a WorkStream context, use this function to append the patch computed
     * by the parallel stage to the array of patches.
     */
    template <int dim, int spacedim>
    void
    append_patch_to_list (const DataOutBase::Patch<dim-1,spacedim> &patch,
                          std::vector<DataOutBase::Patch<dim-1,spacedim> > &patches)
    {
      patches.push_back (patch);
      patches.back().patch_index = patches.size()-1;
    }
  }
}



template <int dim, class DH>
DataOutFaces<dim,DH>::DataOutFaces(const bool so)
  :
  surface_only(so)
{
  Assert (dim == DH::dimension,
          ExcNotImplemented());
}



template <int dim, class DH>
void
DataOutFaces<dim,DH>::
build_one_patch (const FaceDescriptor *cell_and_face,
                 internal::DataOutFaces::ParallelData<dimension, dimension> &data,
                 DataOutBase::Patch<dimension-1,space_dimension>  &patch)
{
  Assert (cell_and_face->first->is_locally_owned(),
          ExcNotImplemented());

  // we use the mapping to transform the vertices. However, the mapping works
  // on cells, not faces, so transform the face vertex to a cell vertex, that
  // to a unit cell vertex and then, finally, that to the mapped vertex. In
  // most cases this complicated procedure will be the identity.
  for (unsigned int vertex=0; vertex<GeometryInfo<dimension-1>::vertices_per_cell; ++vertex)
    patch.vertices[vertex] = data.mapping_collection[0].transform_unit_to_real_cell
                             (cell_and_face->first,
                              GeometryInfo<dimension>::unit_cell_vertex
                              (GeometryInfo<dim>::face_to_cell_vertices
                               (cell_and_face->second,
                                vertex,
                                cell_and_face->first->face_orientation(cell_and_face->second),
                                cell_and_face->first->face_flip(cell_and_face->second),
                                cell_and_face->first->face_rotation(cell_and_face->second))));

  if (data.n_datasets > 0)
    {
      data.reinit_all_fe_values(this->dof_data, cell_and_face->first,
                                cell_and_face->second);
      const FEValuesBase<dimension> &fe_patch_values
        = data.get_present_fe_values (0);

      const unsigned int n_q_points = fe_patch_values.n_quadrature_points;

      // store the intermediate points
      Assert(patch.space_dim==dimension, ExcInternalError());
      const std::vector<Point<dimension> > &q_points=fe_patch_values.get_quadrature_points();
      // resize the patch.data member in order to have enough memory for the
      // quadrature points as well
      patch.data.reinit(data.n_datasets+dimension,
                        patch.data.size(1));
      // set the flag indicating that for this cell the points are explicitly
      // given
      patch.points_are_available=true;
      // copy points to patch.data
      for (unsigned int i=0; i<dimension; ++i)
        for (unsigned int q=0; q<n_q_points; ++q)
          patch.data(patch.data.size(0)-dimension+i,q)=q_points[q][i];

      // counter for data records
      unsigned int offset=0;

      // first fill dof_data
      for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
        {
          const FEValuesBase<dimension> &this_fe_patch_values
            = data.get_present_fe_values (dataset);
          const unsigned int n_components
            = this_fe_patch_values.get_fe().n_components();
          const DataPostprocessor<dim> *postprocessor=this->dof_data[dataset]->postprocessor;
          if (postprocessor != 0)
            {
              // we have to postprocess the data, so determine, which fields
              // have to be updated
              const UpdateFlags update_flags=postprocessor->get_needed_update_flags();

              // get normals, if needed. this is a geometrical information and
              // thus does not depend on the number of components of the data
              // vector
              if (update_flags & update_normal_vectors)
                data.patch_normals=this_fe_patch_values.get_normal_vectors();

              if (n_components == 1)
                {
                  // at each point there is only one component of value,
                  // gradient etc.
                  if (update_flags & update_values)
                    this->dof_data[dataset]->get_function_values (this_fe_patch_values,
                                                                  data.patch_values);
                  if (update_flags & update_gradients)
                    this->dof_data[dataset]->get_function_gradients (this_fe_patch_values,
                                                                     data.patch_gradients);
                  if (update_flags & update_hessians)
                    this->dof_data[dataset]->get_function_hessians (this_fe_patch_values,
                                                                    data.patch_hessians);

                  if (update_flags & update_quadrature_points)
                    data.patch_evaluation_points = this_fe_patch_values.get_quadrature_points();

                  postprocessor->
                  compute_derived_quantities_scalar(data.patch_values,
                                                    data.patch_gradients,
                                                    data.patch_hessians,
                                                    data.patch_normals,
                                                    data.patch_evaluation_points,
                                                    data.postprocessed_values[dataset]);
                }
              else
                {
                  // at each point there is a vector valued function and its
                  // derivative...
                  data.resize_system_vectors(n_components);
                  if (update_flags & update_values)
                    this->dof_data[dataset]->get_function_values (this_fe_patch_values,
                                                                  data.patch_values_system);
                  if (update_flags & update_gradients)
                    this->dof_data[dataset]->get_function_gradients (this_fe_patch_values,
                                                                     data.patch_gradients_system);
                  if (update_flags & update_hessians)
                    this->dof_data[dataset]->get_function_hessians (this_fe_patch_values,
                                                                    data.patch_hessians_system);

                  if (update_flags & update_quadrature_points)
                    data.patch_evaluation_points = this_fe_patch_values.get_quadrature_points();

                  postprocessor->
                  compute_derived_quantities_vector(data.patch_values_system,
                                                    data.patch_gradients_system,
                                                    data.patch_hessians_system,
                                                    data.patch_normals,
                                                    data.patch_evaluation_points,
                                                    data.postprocessed_values[dataset]);
                }

              for (unsigned int q=0; q<n_q_points; ++q)
                for (unsigned int component=0;
                     component<this->dof_data[dataset]->n_output_variables; ++component)
                  patch.data(offset+component,q)
                    = data.postprocessed_values[dataset][q](component);
            }
          else
            // now we use the given data vector without modifications. again,
            // we treat single component functions separately for efficiency
            // reasons.
            if (n_components == 1)
              {
                this->dof_data[dataset]->get_function_values (this_fe_patch_values,
                                                              data.patch_values);
                for (unsigned int q=0; q<n_q_points; ++q)
                  patch.data(offset,q) = data.patch_values[q];
              }
            else
              {
                data.resize_system_vectors(n_components);
                this->dof_data[dataset]->get_function_values (this_fe_patch_values,
                                                              data.patch_values_system);
                for (unsigned int component=0; component<n_components;
                     ++component)
                  for (unsigned int q=0; q<n_q_points; ++q)
                    patch.data(offset+component,q) =
                      data.patch_values_system[q](component);
              }
          // increment the counter for the actual data record
          offset+=this->dof_data[dataset]->n_output_variables;
        }

      // then do the cell data
      for (unsigned int dataset=0; dataset<this->cell_data.size(); ++dataset)
        {
          // we need to get at the number of the cell to which this face
          // belongs in order to access the cell data. this is not readily
          // available, so choose the following rather inefficient way:
          Assert (cell_and_face->first->active(), ExcCellNotActiveForCellData());
          const unsigned int cell_number
            = std::distance (this->triangulation->begin_active(),
                             typename Triangulation<dimension,space_dimension>::active_cell_iterator(cell_and_face->first));

          const double value
            = this->cell_data[dataset]->get_cell_data_value (cell_number);
          for (unsigned int q=0; q<n_q_points; ++q)
            patch.data(dataset+offset,q) = value;
        }
    }
}




template <int dim, class DH>
void DataOutFaces<dim,DH>::build_patches (const unsigned int n_subdivisions_)
{
  build_patches (StaticMappingQ1<dimension>::mapping, n_subdivisions_);
}



template <int dim, class DH>
void DataOutFaces<dim,DH>::build_patches (const Mapping<dimension> &mapping,
                                          const unsigned int n_subdivisions_)
{
  // Check consistency of redundant template parameter
  Assert (dim==dimension, ExcDimensionMismatch(dim, dimension));

  const unsigned int n_subdivisions = (n_subdivisions_ != 0)
                                      ? n_subdivisions_
                                      : this->default_subdivisions;

  Assert (n_subdivisions >= 1,
          ExcInvalidNumberOfSubdivisions(n_subdivisions));

  typedef DataOut_DoFData<DH,dimension+1> BaseClass;
  Assert (this->triangulation != 0,
          typename BaseClass::ExcNoTriangulationSelected());

  unsigned int n_datasets     = this->cell_data.size();
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    n_datasets += this->dof_data[i]->n_output_variables;

  // first count the cells we want to create patches of and make sure there is
  // enough memory for that
  std::vector<FaceDescriptor> all_faces;
  for (FaceDescriptor face=first_face();
       face.first != this->triangulation->end();
       face = next_face(face))
    all_faces.push_back (face);

  // clear the patches array and allocate the right number of elements
  this->patches.clear ();
  this->patches.reserve (all_faces.size());
  Assert (this->patches.size() == 0, ExcInternalError());


  std::vector<unsigned int> n_postprocessor_outputs (this->dof_data.size());
  for (unsigned int dataset=0; dataset<this->dof_data.size(); ++dataset)
    if (this->dof_data[dataset]->postprocessor)
      n_postprocessor_outputs[dataset] = this->dof_data[dataset]->n_output_variables;
    else
      n_postprocessor_outputs[dataset] = 0;

  UpdateFlags update_flags=update_values;
  for (unsigned int i=0; i<this->dof_data.size(); ++i)
    if (this->dof_data[i]->postprocessor)
      update_flags |= this->dof_data[i]->postprocessor->get_needed_update_flags();
  update_flags |= update_quadrature_points;

  internal::DataOutFaces::ParallelData<dimension, space_dimension>
  thread_data (n_datasets,
               n_subdivisions,
               n_postprocessor_outputs,
               mapping,
               this->get_finite_elements(),
               update_flags);
  DataOutBase::Patch<dimension-1,space_dimension> sample_patch;
  sample_patch.n_subdivisions = n_subdivisions;
  sample_patch.data.reinit (n_datasets,
                            Utilities::fixed_power<dimension-1>(n_subdivisions+1));

  // now build the patches in parallel
  WorkStream::run (&all_faces[0],
                   &all_faces[0]+all_faces.size(),
                   std_cxx11::bind(&DataOutFaces<dim,DH>::build_one_patch,
                                   this, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3),
                   std_cxx11::bind(&internal::DataOutFaces::
                                   append_patch_to_list<dim,space_dimension>,
                                   std_cxx11::_1, std_cxx11::ref(this->patches)),
                   thread_data,
                   sample_patch);
}



template <int dim, class DH>
typename DataOutFaces<dim,DH>::FaceDescriptor
DataOutFaces<dim,DH>::first_face ()
{
  // simply find first active cell
  // with a face on the boundary
  typename Triangulation<dimension,space_dimension>::active_cell_iterator cell = this->triangulation->begin_active();
  for (; cell != this->triangulation->end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dimension>::faces_per_cell; ++f)
      if (!surface_only || cell->face(f)->at_boundary())
        return FaceDescriptor(cell, f);

  // ups, triangulation has no
  // boundary? impossible!
  Assert (false, ExcInternalError());

  return FaceDescriptor();
}



template <int dim, class DH>
typename DataOutFaces<dim,DH>::FaceDescriptor
DataOutFaces<dim,DH>::next_face (const FaceDescriptor &old_face)
{
  FaceDescriptor face = old_face;

  // first check whether the present cell has more faces on the boundary
  for (unsigned int f=face.second+1; f<GeometryInfo<dimension>::faces_per_cell; ++f)
    if (!surface_only || face.first->face(f)->at_boundary())
      // yup, that is so, so return it
      {
        face.second = f;
        return face;
      };

  // otherwise find the next active cell that has a face on the boundary

  // convert the iterator to an active_iterator and advance this to the next
  // active cell
  typename Triangulation<dimension,space_dimension>::active_cell_iterator active_cell = face.first;

  // increase face pointer by one
  ++active_cell;

  // while there are active cells
  while (active_cell != this->triangulation->end())
    {
      // check all the faces of this active cell
      for (unsigned int f=0; f<GeometryInfo<dimension>::faces_per_cell; ++f)
        if (!surface_only || active_cell->face(f)->at_boundary())
          {
            face.first  = active_cell;
            face.second = f;
            return face;
          };
      // the present cell had no faces on the boundary, so check next cell
      ++active_cell;
    };

  // we fell off the edge, so return with invalid pointer
  face.first  = this->triangulation->end();
  face.second = 0;
  return face;
}



// explicit instantiations
#include "data_out_faces.inst"

DEAL_II_NAMESPACE_CLOSE
