// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out_faces.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutFacesImplementation
  {
    template <int dim, int spacedim>
    ParallelData<dim, spacedim>::ParallelData(
      const unsigned int               n_datasets,
      const unsigned int               n_subdivisions,
      const std::vector<unsigned int> &n_postprocessor_outputs,
      const Mapping<dim, spacedim>    &mapping,
      const std::vector<
        std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                       &finite_elements,
      const UpdateFlags update_flags)
      : internal::DataOutImplementation::ParallelDataBase<dim, spacedim>(
          n_datasets,
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
    append_patch_to_list(
      const DataOutBase::Patch<DataOutFaces<dim, spacedim>::patch_dim,
                               DataOutFaces<dim, spacedim>::patch_spacedim>
        &patch,
      std::vector<
        DataOutBase::Patch<DataOutFaces<dim, spacedim>::patch_dim,
                           DataOutFaces<dim, spacedim>::patch_spacedim>>
        &patches)
    {
      patches.push_back(patch);
      patches.back().patch_index = patches.size() - 1;
    }
  } // namespace DataOutFacesImplementation
} // namespace internal



template <int dim, int spacedim>
DataOutFaces<dim, spacedim>::DataOutFaces(const bool so)
  : surface_only(so)
{}



template <int dim, int spacedim>
void
DataOutFaces<dim, spacedim>::build_one_patch(
  const FaceDescriptor *cell_and_face,
  internal::DataOutFacesImplementation::ParallelData<dim, spacedim> &data,
  DataOutBase::Patch<patch_dim, patch_spacedim>                     &patch)
{
  const cell_iterator cell        = cell_and_face->first;
  const unsigned int  face_number = cell_and_face->second;

  Assert(cell->is_locally_owned(), ExcNotImplemented());

  // First set the kind of object we are dealing with here in the 'patch'
  // object.
  patch.reference_cell = cell->face(face_number)->reference_cell();

  // We use the mapping to transform the vertices. However, the mapping works
  // on cells, not faces, so transform the face vertex to a cell vertex, that
  // to a unit cell vertex and then, finally, that to the mapped vertex. In
  // most cases this complicated procedure will be the identity.
  for (const unsigned int vertex : cell->face(face_number)->vertex_indices())
    {
      const Point<dim> vertex_reference_coordinates =
        cell->reference_cell().template vertex<dim>(
          cell->reference_cell().face_to_cell_vertices(
            face_number, vertex, cell->combined_face_orientation(face_number)));

      const Point<dim> vertex_real_coordinates =
        data.mapping_collection[0].transform_unit_to_real_cell(
          cell, vertex_reference_coordinates);

      patch.vertices[vertex] = vertex_real_coordinates;
    }


  if (data.n_datasets > 0)
    {
      data.reinit_all_fe_values(this->dof_data, cell, face_number);
      const FEValuesBase<dim> &fe_patch_values = data.get_present_fe_values(0);

      const unsigned int n_q_points = fe_patch_values.n_quadrature_points;

      // store the intermediate points
      Assert(patch.space_dim == dim, ExcInternalError());
      const std::vector<Point<dim>> &q_points =
        fe_patch_values.get_quadrature_points();
      // size the patch.data member in order to have enough memory for the
      // quadrature points as well
      patch.data.reinit(data.n_datasets + dim, q_points.size());
      // set the flag indicating that for this cell the points are explicitly
      // given
      patch.points_are_available = true;
      // copy points to patch.data
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int q = 0; q < n_q_points; ++q)
          patch.data(patch.data.size(0) - dim + i, q) = q_points[q][i];

      // counter for data records
      unsigned int offset = 0;

      // first fill dof_data
      for (unsigned int dataset = 0; dataset < this->dof_data.size(); ++dataset)
        {
          const FEValuesBase<dim> &this_fe_patch_values =
            data.get_present_fe_values(dataset);
          const unsigned int n_components =
            this_fe_patch_values.get_fe().n_components();
          const DataPostprocessor<dim> *postprocessor =
            this->dof_data[dataset]->postprocessor;
          if (postprocessor != nullptr)
            {
              // we have to postprocess the data, so determine, which fields
              // have to be updated
              const UpdateFlags update_flags =
                postprocessor->get_needed_update_flags();

              if (n_components == 1)
                {
                  // at each point there is only one component of value,
                  // gradient etc.
                  if (update_flags & update_values)
                    this->dof_data[dataset]->get_function_values(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      data.patch_values_scalar.solution_values);
                  if (update_flags & update_gradients)
                    this->dof_data[dataset]->get_function_gradients(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      data.patch_values_scalar.solution_gradients);
                  if (update_flags & update_hessians)
                    this->dof_data[dataset]->get_function_hessians(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      data.patch_values_scalar.solution_hessians);

                  if (update_flags & update_quadrature_points)
                    data.patch_values_scalar.evaluation_points =
                      this_fe_patch_values.get_quadrature_points();

                  if (update_flags & update_normal_vectors)
                    data.patch_values_scalar.normals =
                      this_fe_patch_values.get_normal_vectors();

                  const typename DoFHandler<dim, spacedim>::active_cell_iterator
                    dh_cell(&cell->get_triangulation(),
                            cell->level(),
                            cell->index(),
                            this->dof_data[dataset]->dof_handler);
                  data.patch_values_scalar.template set_cell_and_face<dim>(
                    dh_cell, face_number);

                  postprocessor->evaluate_scalar_field(
                    data.patch_values_scalar,
                    data.postprocessed_values[dataset]);
                }
              else
                {
                  // at each point there is a vector valued function and its
                  // derivative...
                  data.resize_system_vectors(n_components);
                  if (update_flags & update_values)
                    this->dof_data[dataset]->get_function_values(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      data.patch_values_system.solution_values);
                  if (update_flags & update_gradients)
                    this->dof_data[dataset]->get_function_gradients(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      data.patch_values_system.solution_gradients);
                  if (update_flags & update_hessians)
                    this->dof_data[dataset]->get_function_hessians(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      data.patch_values_system.solution_hessians);

                  if (update_flags & update_quadrature_points)
                    data.patch_values_system.evaluation_points =
                      this_fe_patch_values.get_quadrature_points();

                  if (update_flags & update_normal_vectors)
                    data.patch_values_system.normals =
                      this_fe_patch_values.get_normal_vectors();

                  const typename DoFHandler<dim, spacedim>::active_cell_iterator
                    dh_cell(&cell->get_triangulation(),
                            cell->level(),
                            cell->index(),
                            this->dof_data[dataset]->dof_handler);
                  data.patch_values_system.template set_cell_and_face<dim>(
                    dh_cell, face_number);

                  postprocessor->evaluate_vector_field(
                    data.patch_values_system,
                    data.postprocessed_values[dataset]);
                }

              for (unsigned int q = 0; q < n_q_points; ++q)
                for (unsigned int component = 0;
                     component < this->dof_data[dataset]->n_output_variables;
                     ++component)
                  patch.data(offset + component, q) =
                    data.postprocessed_values[dataset][q](component);
            }
          else
            // now we use the given data vector without modifications. again,
            // we treat single component functions separately for efficiency
            // reasons.
            if (n_components == 1)
              {
                this->dof_data[dataset]->get_function_values(
                  this_fe_patch_values,
                  internal::DataOutImplementation::ComponentExtractor::
                    real_part,
                  data.patch_values_scalar.solution_values);
                for (unsigned int q = 0; q < n_q_points; ++q)
                  patch.data(offset, q) =
                    data.patch_values_scalar.solution_values[q];
              }
            else
              {
                data.resize_system_vectors(n_components);
                this->dof_data[dataset]->get_function_values(
                  this_fe_patch_values,
                  internal::DataOutImplementation::ComponentExtractor::
                    real_part,
                  data.patch_values_system.solution_values);
                for (unsigned int component = 0; component < n_components;
                     ++component)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    patch.data(offset + component, q) =
                      data.patch_values_system.solution_values[q](component);
              }
          // increment the counter for the actual data record
          offset += this->dof_data[dataset]->n_output_variables;
        }

      // then do the cell data
      for (unsigned int dataset = 0; dataset < this->cell_data.size();
           ++dataset)
        {
          // we need to get at the number of the cell to which this face
          // belongs in order to access the cell data. this is not readily
          // available, so choose the following rather inefficient way:
          Assert(
            cell->is_active(),
            ExcMessage(
              "The current function is trying to generate cell-data output "
              "for a face that does not belong to an active cell. This is "
              "not supported."));
          const unsigned int cell_number = std::distance(
            this->triangulation->begin_active(),
            typename Triangulation<dim, spacedim>::active_cell_iterator(cell));

          const double value = this->cell_data[dataset]->get_cell_data_value(
            cell_number,
            internal::DataOutImplementation::ComponentExtractor::real_part);
          for (unsigned int q = 0; q < n_q_points; ++q)
            patch.data(dataset + offset, q) = value;
        }
    }
}



template <int dim, int spacedim>
void
DataOutFaces<dim, spacedim>::build_patches(const unsigned int n_subdivisions)
{
  if (this->triangulation->get_reference_cells().size() == 1)
    build_patches(this->triangulation->get_reference_cells()[0]
                    .template get_default_linear_mapping<dim, spacedim>(),
                  n_subdivisions);
  else
    Assert(false,
           ExcMessage("The DataOutFaces class can currently not be "
                      "used on meshes that do not have the same cell type "
                      "throughout."));
}



template <int dim, int spacedim>
void
DataOutFaces<dim, spacedim>::build_patches(
  const Mapping<dim, spacedim> &mapping,
  const unsigned int            n_subdivisions_)
{
  const unsigned int n_subdivisions =
    (n_subdivisions_ != 0) ? n_subdivisions_ : this->default_subdivisions;

  Assert(n_subdivisions >= 1,
         Exceptions::DataOutImplementation::ExcInvalidNumberOfSubdivisions(
           n_subdivisions));

  Assert(this->triangulation != nullptr,
         Exceptions::DataOutImplementation::ExcNoTriangulationSelected());

  this->validate_dataset_names();

  unsigned int n_datasets = this->cell_data.size();
  for (unsigned int i = 0; i < this->dof_data.size(); ++i)
    n_datasets += this->dof_data[i]->n_output_variables;

  // first collect the cells we want to create patches of; we will
  // then iterate over them. the end-condition of the loop needs to
  // test that next_face() returns an end iterator, as well as for the
  // case that first_face() returns an invalid FaceDescriptor object
  std::vector<FaceDescriptor> all_faces;
  for (FaceDescriptor face = first_face();
       ((face.first != this->triangulation->end()) &&
        (face != FaceDescriptor()));
       face = next_face(face))
    all_faces.push_back(face);

  // clear the patches array and allocate the right number of elements
  this->patches.clear();
  this->patches.reserve(all_faces.size());
  Assert(this->patches.empty(), ExcInternalError());


  std::vector<unsigned int> n_postprocessor_outputs(this->dof_data.size());
  for (unsigned int dataset = 0; dataset < this->dof_data.size(); ++dataset)
    if (this->dof_data[dataset]->postprocessor)
      n_postprocessor_outputs[dataset] =
        this->dof_data[dataset]->n_output_variables;
    else
      n_postprocessor_outputs[dataset] = 0;

  UpdateFlags update_flags = update_values;
  for (unsigned int i = 0; i < this->dof_data.size(); ++i)
    if (this->dof_data[i]->postprocessor)
      update_flags |=
        this->dof_data[i]->postprocessor->get_needed_update_flags();
  update_flags |= update_quadrature_points;

  internal::DataOutFacesImplementation::ParallelData<dim, spacedim> thread_data(
    n_datasets,
    n_subdivisions,
    n_postprocessor_outputs,
    mapping,
    this->get_fes(),
    update_flags);
  DataOutBase::Patch<patch_dim, patch_spacedim> sample_patch;
  sample_patch.n_subdivisions = n_subdivisions;

  // now build the patches in parallel
  WorkStream::run(
    all_faces.data(),
    all_faces.data() + all_faces.size(),
    [this](
      const FaceDescriptor *cell_and_face,
      internal::DataOutFacesImplementation::ParallelData<dim, spacedim> &data,
      DataOutBase::Patch<patch_dim, patch_spacedim> &patch) {
      this->build_one_patch(cell_and_face, data, patch);
    },
    [this](const DataOutBase::Patch<patch_dim, patch_spacedim> &patch) {
      internal::DataOutFacesImplementation::append_patch_to_list<dim, spacedim>(
        patch, this->patches);
    },
    thread_data,
    sample_patch);
}



template <int dim, int spacedim>
typename DataOutFaces<dim, spacedim>::FaceDescriptor
DataOutFaces<dim, spacedim>::first_face()
{
  // simply find first active cell with a face on the boundary
  for (const auto &cell : this->triangulation->active_cell_iterators())
    if (cell->is_locally_owned())
      for (const unsigned int f : cell->face_indices())
        if ((surface_only == false) || cell->face(f)->at_boundary())
          return FaceDescriptor(cell, f);

  // just return an invalid descriptor if we haven't found a locally
  // owned face. this can happen in parallel where all boundary
  // faces are owned by other processors
  return FaceDescriptor();
}



template <int dim, int spacedim>
typename DataOutFaces<dim, spacedim>::FaceDescriptor
DataOutFaces<dim, spacedim>::next_face(const FaceDescriptor &old_face)
{
  FaceDescriptor face = old_face;

  // first check whether the present cell has more faces on the boundary. since
  // we started with this face, its cell must clearly be locally owned
  Assert(face.first->is_locally_owned(), ExcInternalError());
  for (unsigned int f = face.second + 1; f < face.first->n_faces(); ++f)
    if (!surface_only || face.first->face(f)->at_boundary())
      // yup, that is so, so return it
      {
        face.second = f;
        return face;
      }

  // otherwise find the next active cell that has a face on the boundary

  // convert the iterator to an active_iterator and advance this to the next
  // active cell
  typename Triangulation<dim, spacedim>::active_cell_iterator active_cell =
    face.first;

  // increase face pointer by one
  ++active_cell;

  // while there are active cells
  while (active_cell != this->triangulation->end())
    {
      // check all the faces of this active cell. but skip it altogether
      // if it isn't locally owned
      if (active_cell->is_locally_owned())
        for (const unsigned int f : face.first->face_indices())
          if (!surface_only || active_cell->face(f)->at_boundary())
            {
              face.first  = active_cell;
              face.second = f;
              return face;
            }

      // the present cell had no faces on the boundary (or was not locally
      // owned), so check next cell
      ++active_cell;
    }

  // we fell off the edge, so return with invalid pointer
  face.first  = this->triangulation->end();
  face.second = 0;
  return face;
}



// explicit instantiations
#include "numerics/data_out_faces.inst"

DEAL_II_NAMESPACE_CLOSE
