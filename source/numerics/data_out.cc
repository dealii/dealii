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

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/data_out.h>

#include <sstream>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutImplementation
  {
    template <int dim, int spacedim>
    ParallelData<dim, spacedim>::ParallelData(
      const unsigned int               n_datasets,
      const unsigned int               n_subdivisions,
      const std::vector<unsigned int> &n_postprocessor_outputs,
      const Mapping<dim, spacedim> &   mapping,
      const std::vector<
        std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
        &                                           finite_elements,
      const UpdateFlags                             update_flags,
      const std::vector<std::vector<unsigned int>> &cell_to_patch_index_map)
      : ParallelDataBase<dim, spacedim>(n_datasets,
                                        n_subdivisions,
                                        n_postprocessor_outputs,
                                        mapping,
                                        finite_elements,
                                        update_flags,
                                        false)
      , cell_to_patch_index_map(&cell_to_patch_index_map)
    {}
  } // namespace DataOutImplementation
} // namespace internal



template <int dim, typename DoFHandlerType>
DataOut<dim, DoFHandlerType>::DataOut()
{
  // For the moment, just call the existing virtual functions. This
  // preserves backward compatibility. When these deprecated functions are
  // removed, we'll just inline their functionality into the lambda
  // functions created here.
  set_cell_selection(
    [this](const Triangulation<dim, spacedim> &) {
      return this->first_locally_owned_cell();
    },
    [this](const Triangulation<dim, spacedim> &, const cell_iterator &cell) {
      return this->next_locally_owned_cell(cell);
    });
}



template <int dim, typename DoFHandlerType>
void
DataOut<dim, DoFHandlerType>::build_one_patch(
  const std::pair<cell_iterator, unsigned int> *cell_and_index,
  internal::DataOutImplementation::ParallelData<DoFHandlerType::dimension,
                                                DoFHandlerType::space_dimension>
    &                    scratch_data,
  const unsigned int     n_subdivisions,
  const CurvedCellRegion curved_cell_region)
{
  // first create the output object that we will write into
  ::dealii::DataOutBase::Patch<DoFHandlerType::dimension,
                               DoFHandlerType::space_dimension>
    patch;
  patch.n_subdivisions = n_subdivisions;

  // set the vertices of the patch. if the mapping does not preserve locations
  // (e.g. MappingQEulerian), we need to compute the offset of the vertex for
  // the graphical output. Otherwise, we can just use the vertex info.
  for (const unsigned int vertex :
       GeometryInfo<DoFHandlerType::dimension>::vertex_indices())
    if (scratch_data.mapping_collection[0].preserves_vertex_locations())
      patch.vertices[vertex] = cell_and_index->first->vertex(vertex);
    else
      patch.vertices[vertex] =
        scratch_data.mapping_collection[0].transform_unit_to_real_cell(
          cell_and_index->first,
          GeometryInfo<DoFHandlerType::dimension>::unit_cell_vertex(vertex));

  // initialize FEValues
  scratch_data.reinit_all_fe_values(this->dof_data, cell_and_index->first);

  const FEValuesBase<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
    &fe_patch_values = scratch_data.get_present_fe_values(0);

  const unsigned int n_q_points = fe_patch_values.n_quadrature_points;

  // First fill the geometric information for the patch: Where are the
  // nodes in question located.
  //
  // Depending on the requested output of curved cells, if necessary
  // append the quadrature points to the last rows of the patch.data
  // member. This is the case if we want to produce curved cells at the
  // boundary and this cell actually is at the boundary, or else if we
  // want to produce curved cells everywhere
  //
  // note: a cell is *always* at the boundary if dim<spacedim
  if (curved_cell_region == curved_inner_cells ||
      (curved_cell_region == curved_boundary &&
       (cell_and_index->first->at_boundary() ||
        (DoFHandlerType::dimension != DoFHandlerType::space_dimension))))
    {
      Assert(patch.space_dim == DoFHandlerType::space_dimension,
             ExcInternalError());

      // set the flag indicating that for this cell the points are
      // explicitly given
      patch.points_are_available = true;

      // then resize the patch.data member in order to have enough memory for
      // the quadrature points as well, and copy the quadrature points there
      const std::vector<Point<DoFHandlerType::space_dimension>> &q_points =
        fe_patch_values.get_quadrature_points();

      patch.data.reinit(scratch_data.n_datasets +
                          DoFHandlerType::space_dimension,
                        n_q_points);
      for (unsigned int i = 0; i < DoFHandlerType::space_dimension; ++i)
        for (unsigned int q = 0; q < n_q_points; ++q)
          patch.data(patch.data.size(0) - DoFHandlerType::space_dimension + i,
                     q) = q_points[q][i];
    }
  else
    {
      patch.data.reinit(scratch_data.n_datasets, n_q_points);
      patch.points_are_available = false;
    }


  // Next fill the information we get from DoF data
  if (scratch_data.n_datasets > 0)
    {
      // counter for data records
      unsigned int offset = 0;

      // first fill dof_data
      unsigned int dataset_number = 0;
      for (const auto &dataset : this->dof_data)
        {
          const FEValuesBase<DoFHandlerType::dimension,
                             DoFHandlerType::space_dimension>
            &this_fe_patch_values =
              scratch_data.get_present_fe_values(dataset_number);
          const unsigned int n_components =
            this_fe_patch_values.get_fe().n_components();

          const DataPostprocessor<DoFHandlerType::space_dimension>
            *postprocessor = dataset->postprocessor;

          if (postprocessor != nullptr)
            {
              // We have to postprocess the data, so determine, which fields
              // have to be updated
              const UpdateFlags update_flags =
                postprocessor->get_needed_update_flags();

              if ((n_components == 1) &&
                  (dataset->is_complex_valued() == false))
                {
                  // At each point there is only one component of value,
                  // gradient etc. Based on the 'if' statement above, we
                  // know that the solution is scalar and real-valued, so
                  // we do not need to worry about getting any imaginary
                  // components to the postprocessor, and we can safely
                  // call the function that evaluates a scalar field
                  if (update_flags & update_values)
                    dataset->get_function_values(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      scratch_data.patch_values_scalar.solution_values);

                  if (update_flags & update_gradients)
                    dataset->get_function_gradients(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      scratch_data.patch_values_scalar.solution_gradients);

                  if (update_flags & update_hessians)
                    dataset->get_function_hessians(
                      this_fe_patch_values,
                      internal::DataOutImplementation::ComponentExtractor::
                        real_part,
                      scratch_data.patch_values_scalar.solution_hessians);

                  // Also fill some of the other fields postprocessors may
                  // want to access.
                  if (update_flags & update_quadrature_points)
                    scratch_data.patch_values_scalar.evaluation_points =
                      this_fe_patch_values.get_quadrature_points();

                  const typename DoFHandlerType::cell_iterator dh_cell(
                    &cell_and_index->first->get_triangulation(),
                    cell_and_index->first->level(),
                    cell_and_index->first->index(),
                    dataset->dof_handler);
                  scratch_data.patch_values_scalar
                    .template set_cell<DoFHandlerType>(dh_cell);

                  // Finally call the postprocessor's function that
                  // deals with scalar inputs.
                  postprocessor->evaluate_scalar_field(
                    scratch_data.patch_values_scalar,
                    scratch_data.postprocessed_values[dataset_number]);
                }
              else
                {
                  // At each point we now have to evaluate a vector valued
                  // function and its derivatives. It may be that the solution
                  // is scalar and complex-valued, but we treat this as a vector
                  // field with two components.

                  // At this point, we need to ask how we fill the fields that
                  // we want to pass on to the postprocessor. If the field in
                  // question is real-valued, we'll just extract the (only)
                  // real component from the solution fields
                  if (dataset->is_complex_valued() == false)
                    {
                      scratch_data.resize_system_vectors(n_components);

                      if (update_flags & update_values)
                        dataset->get_function_values(
                          this_fe_patch_values,
                          internal::DataOutImplementation::ComponentExtractor::
                            real_part,
                          scratch_data.patch_values_system.solution_values);

                      if (update_flags & update_gradients)
                        dataset->get_function_gradients(
                          this_fe_patch_values,
                          internal::DataOutImplementation::ComponentExtractor::
                            real_part,
                          scratch_data.patch_values_system.solution_gradients);

                      if (update_flags & update_hessians)
                        dataset->get_function_hessians(
                          this_fe_patch_values,
                          internal::DataOutImplementation::ComponentExtractor::
                            real_part,
                          scratch_data.patch_values_system.solution_hessians);
                    }
                  else
                    {
                      // The solution is complex-valued. Let's cover the scalar
                      // case first (i.e., one scalar but complex-valued field,
                      // which we will have to split into its real and imaginar
                      // parts).
                      if (n_components == 1)
                        {
                          scratch_data.resize_system_vectors(2);

                          // First get the real component of the scalar solution
                          // and copy the data into the
                          // scratch_data.patch_values_system output fields
                          if (update_flags & update_values)
                            {
                              dataset->get_function_values(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::real_part,
                                scratch_data.patch_values_scalar
                                  .solution_values);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_scalar
                                         .solution_values.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_values[i]
                                      .size(),
                                    2);
                                  scratch_data.patch_values_system
                                    .solution_values[i][0] =
                                    scratch_data.patch_values_scalar
                                      .solution_values[i];
                                }
                            }

                          if (update_flags & update_gradients)
                            {
                              dataset->get_function_gradients(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::real_part,
                                scratch_data.patch_values_scalar
                                  .solution_gradients);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_scalar
                                         .solution_gradients.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_values[i]
                                      .size(),
                                    2);
                                  scratch_data.patch_values_system
                                    .solution_gradients[i][0] =
                                    scratch_data.patch_values_scalar
                                      .solution_gradients[i];
                                }
                            }

                          if (update_flags & update_hessians)
                            {
                              dataset->get_function_hessians(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::real_part,
                                scratch_data.patch_values_scalar
                                  .solution_hessians);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_scalar
                                         .solution_hessians.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_hessians[i]
                                      .size(),
                                    2);
                                  scratch_data.patch_values_system
                                    .solution_hessians[i][0] =
                                    scratch_data.patch_values_scalar
                                      .solution_hessians[i];
                                }
                            }

                          // Now we also have to get the imaginary
                          // component of the scalar solution
                          // and copy the data into the
                          // scratch_data.patch_values_system output fields
                          // that follow the real one
                          if (update_flags & update_values)
                            {
                              dataset->get_function_values(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::imaginary_part,
                                scratch_data.patch_values_scalar
                                  .solution_values);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_scalar
                                         .solution_values.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_values[i]
                                      .size(),
                                    2);
                                  scratch_data.patch_values_system
                                    .solution_values[i][1] =
                                    scratch_data.patch_values_scalar
                                      .solution_values[i];
                                }
                            }

                          if (update_flags & update_gradients)
                            {
                              dataset->get_function_gradients(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::imaginary_part,
                                scratch_data.patch_values_scalar
                                  .solution_gradients);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_scalar
                                         .solution_gradients.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_values[i]
                                      .size(),
                                    2);
                                  scratch_data.patch_values_system
                                    .solution_gradients[i][1] =
                                    scratch_data.patch_values_scalar
                                      .solution_gradients[i];
                                }
                            }

                          if (update_flags & update_hessians)
                            {
                              dataset->get_function_hessians(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::imaginary_part,
                                scratch_data.patch_values_scalar
                                  .solution_hessians);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_scalar
                                         .solution_hessians.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_hessians[i]
                                      .size(),
                                    2);
                                  scratch_data.patch_values_system
                                    .solution_hessians[i][1] =
                                    scratch_data.patch_values_scalar
                                      .solution_hessians[i];
                                }
                            }
                        }
                      else
                        {
                          scratch_data.resize_system_vectors(2 * n_components);

                          // This is the vector-valued, complex-valued case. In
                          // essence, we just need to do the same as above,
                          // i.e., call the functions in dataset
                          // to retrieve first the real and then the imaginary
                          // part of the solution, then copy them to the
                          // scratch_data.patch_values_system. The difference to
                          // the scalar case is that there, we could (ab)use the
                          // scratch_data.patch_values_scalar members for first
                          // the real part and then the imaginary part, copying
                          // them into the scratch_data.patch_values_system
                          // variable one after the other. We can't do this
                          // here because the solution is vector-valued, and
                          // so using the single
                          // scratch_data.patch_values_system object doesn't
                          // work.
                          //
                          // Rather, we need to come up with a temporary object
                          // for this (one for each the values, the gradients,
                          // and the hessians).
                          //
                          // Compared to the previous code path, we here
                          // first get the real and imaginary parts of the
                          // values, then the real and imaginary parts of the
                          // gradients, etc. This allows us to scope the
                          // temporary objects better
                          if (update_flags & update_values)
                            {
                              std::vector<Vector<double>> tmp(
                                scratch_data.patch_values_system.solution_values
                                  .size(),
                                Vector<double>(n_components));

                              // First get the real part into the tmp object
                              dataset->get_function_values(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::real_part,
                                tmp);

                              // Then copy these values into the first
                              // n_components slots of the output object.
                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_system
                                         .solution_values.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_values[i]
                                      .size(),
                                    2 * n_components);
                                  for (unsigned int j = 0; j < n_components;
                                       ++j)
                                    scratch_data.patch_values_system
                                      .solution_values[i][j] = tmp[i][j];
                                }

                              // Then do the same with the imaginary part,
                              // copying past the end of the previous set of
                              // values.
                              dataset->get_function_values(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::imaginary_part,
                                tmp);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_system
                                         .solution_values.size();
                                   ++i)
                                {
                                  for (unsigned int j = 0; j < n_components;
                                       ++j)
                                    scratch_data.patch_values_system
                                      .solution_values[i][j + n_components] =
                                      tmp[i][j];
                                }
                            }

                          // Now do the exact same thing for the gradients
                          if (update_flags & update_gradients)
                            {
                              std::vector<std::vector<
                                Tensor<1, DoFHandlerType::space_dimension>>>
                                tmp(
                                  scratch_data.patch_values_system
                                    .solution_gradients.size(),
                                  std::vector<
                                    Tensor<1, DoFHandlerType::space_dimension>>(
                                    n_components));

                              // First the real part
                              dataset->get_function_gradients(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::real_part,
                                tmp);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_system
                                         .solution_gradients.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_gradients[i]
                                      .size(),
                                    2 * n_components);
                                  for (unsigned int j = 0; j < n_components;
                                       ++j)
                                    scratch_data.patch_values_system
                                      .solution_gradients[i][j] = tmp[i][j];
                                }

                              // Then the imaginary part
                              dataset->get_function_gradients(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::imaginary_part,
                                tmp);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_system
                                         .solution_gradients.size();
                                   ++i)
                                {
                                  for (unsigned int j = 0; j < n_components;
                                       ++j)
                                    scratch_data.patch_values_system
                                      .solution_gradients[i][j + n_components] =
                                      tmp[i][j];
                                }
                            }

                          // And finally the Hessians. Same scheme as above.
                          if (update_flags & update_hessians)
                            {
                              std::vector<std::vector<
                                Tensor<2, DoFHandlerType::space_dimension>>>
                                tmp(
                                  scratch_data.patch_values_system
                                    .solution_gradients.size(),
                                  std::vector<
                                    Tensor<2, DoFHandlerType::space_dimension>>(
                                    n_components));

                              // First the real part
                              dataset->get_function_hessians(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::real_part,
                                tmp);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_system
                                         .solution_hessians.size();
                                   ++i)
                                {
                                  AssertDimension(
                                    scratch_data.patch_values_system
                                      .solution_hessians[i]
                                      .size(),
                                    2 * n_components);
                                  for (unsigned int j = 0; j < n_components;
                                       ++j)
                                    scratch_data.patch_values_system
                                      .solution_hessians[i][j] = tmp[i][j];
                                }

                              // Then the imaginary part
                              dataset->get_function_hessians(
                                this_fe_patch_values,
                                internal::DataOutImplementation::
                                  ComponentExtractor::imaginary_part,
                                tmp);

                              for (unsigned int i = 0;
                                   i < scratch_data.patch_values_system
                                         .solution_hessians.size();
                                   ++i)
                                {
                                  for (unsigned int j = 0; j < n_components;
                                       ++j)
                                    scratch_data.patch_values_system
                                      .solution_hessians[i][j + n_components] =
                                      tmp[i][j];
                                }
                            }
                        }
                    }

                  // Now set other fields we may need
                  if (update_flags & update_quadrature_points)
                    scratch_data.patch_values_system.evaluation_points =
                      this_fe_patch_values.get_quadrature_points();

                  const typename DoFHandlerType::cell_iterator dh_cell(
                    &cell_and_index->first->get_triangulation(),
                    cell_and_index->first->level(),
                    cell_and_index->first->index(),
                    dataset->dof_handler);
                  scratch_data.patch_values_system
                    .template set_cell<DoFHandlerType>(dh_cell);

                  // Whether the solution was complex-scalar or
                  // complex-vector-valued doesn't matter -- we took it apart
                  // into several fields and so we have to call the
                  // evaluate_vector_field() function.
                  postprocessor->evaluate_vector_field(
                    scratch_data.patch_values_system,
                    scratch_data.postprocessed_values[dataset_number]);
                }

              // Now we need to copy the result of the postprocessor to
              // the Patch object where it can then be further processed
              // by the functions in DataOutBase
              for (unsigned int q = 0; q < n_q_points; ++q)
                for (unsigned int component = 0;
                     component < dataset->n_output_variables;
                     ++component)
                  patch.data(offset + component, q) =
                    scratch_data.postprocessed_values[dataset_number][q](
                      component);

              // Move the counter for the output location forward as
              // appropriate
              offset += dataset->n_output_variables;
            }
          else
            {
              // use the given data vector directly, without a postprocessor.
              // again, we treat single component functions separately for
              // efficiency reasons.
              if (n_components == 1)
                {
                  Assert(dataset->n_output_variables == 1, ExcInternalError());

                  // First output the real part of the solution vector
                  dataset->get_function_values(
                    this_fe_patch_values,
                    internal::DataOutImplementation::ComponentExtractor::
                      real_part,
                    scratch_data.patch_values_scalar.solution_values);
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    patch.data(offset, q) =
                      scratch_data.patch_values_scalar.solution_values[q];
                  offset += 1;

                  // And if there is one, also output the imaginary part. Note
                  // that the problem is scalar-valued, so we can freely add the
                  // imaginary part after the real part without having to worry
                  // that we are interleaving the real components of a vector
                  // with the imaginary components of the same vector.
                  if (dataset->is_complex_valued() == true)
                    {
                      dataset->get_function_values(
                        this_fe_patch_values,
                        internal::DataOutImplementation::ComponentExtractor::
                          imaginary_part,
                        scratch_data.patch_values_scalar.solution_values);
                      for (unsigned int q = 0; q < n_q_points; ++q)
                        patch.data(offset, q) =
                          scratch_data.patch_values_scalar.solution_values[q];
                      offset += 1;
                    }
                }
              else
                {
                  scratch_data.resize_system_vectors(n_components);

                  // So we have a multi-component DoFHandler here. That's more
                  // complicated. If the vector is real-valued, then we can just
                  // get everything at all quadrature points and copy them into
                  // the output array. In fact, we don't have to worry at all
                  // about the interpretation of the components.
                  if (dataset->is_complex_valued() == false)
                    {
                      dataset->get_function_values(
                        this_fe_patch_values,
                        internal::DataOutImplementation::ComponentExtractor::
                          real_part,
                        scratch_data.patch_values_system.solution_values);
                      for (unsigned int component = 0; component < n_components;
                           ++component)
                        for (unsigned int q = 0; q < n_q_points; ++q)
                          patch.data(offset + component, q) =
                            scratch_data.patch_values_system.solution_values[q](
                              component);

                      // Increment the counter for the actual data record.
                      offset += dataset->n_output_variables;
                    }
                  else
                    // The situation is more complicated if the input vector is
                    // complex-valued. The easiest approach would have been to
                    // just have all real and then all imaginary components.
                    // This would have been conceptually easy, but it has the
                    // annoying downside that if you have a vector-valued
                    // problem (say, [u v]) then the output order would have
                    // been [u_re, v_re, u_im, v_im]. That's tolerable, but not
                    // quite so nice because one typically thinks of real and
                    // imaginary parts as belonging together. We would really
                    // like the output order to be [u_re, u_im, v_re, v_im].
                    // That, too, would have been easy to implement because one
                    // just has to interleave real and imaginary parts.
                    //
                    // But that's also not what we want. That's because if one
                    // were, for example, to solve a complex-valued Stokes
                    // problem (e.g., computing eigenfunctions of the Stokes
                    // operator), then one has solution components
                    // [[u v] p] and the proper output order is
                    //     [[u_re v_re] [u_im v_im] p_re p_im].
                    // In other words, the order in which we want to output
                    // data depends on the *interpretation* of components.
                    //
                    // Doing this requires a bit more code, and also needs to
                    // be in sync with what we do in
                    // DataOut_DoFData::get_dataset_names() and
                    // DataOut_DoFData::get_nonscalar_data_ranges().
                    {
                      // Given this description, first get the real parts of
                      // all components:
                      dataset->get_function_values(
                        this_fe_patch_values,
                        internal::DataOutImplementation::ComponentExtractor::
                          real_part,
                        scratch_data.patch_values_system.solution_values);

                      // Then we need to distribute them to the correct
                      // location. This requires knowledge of the interpretation
                      // of components as discussed above.
                      {
                        Assert(dataset->data_component_interpretation.size() ==
                                 n_components,
                               ExcInternalError());

                        unsigned int destination = offset;
                        for (unsigned int component = 0;
                             component < n_components;
                             /* component is updated below */)
                          {
                            switch (
                              dataset->data_component_interpretation[component])
                              {
                                case DataComponentInterpretation::
                                  component_is_scalar:
                                  {
                                    // OK, a scalar component. Put all of the
                                    // values into the current row
                                    // ('destination'); then move 'component'
                                    // forward by one (so we treat the next
                                    // component) and 'destination' forward by
                                    // two (because we're going to put the
                                    // imaginary part of the current component
                                    // into the next slot).
                                    for (unsigned int q = 0; q < n_q_points;
                                         ++q)
                                      patch.data(destination, q) =
                                        scratch_data.patch_values_system
                                          .solution_values[q](component);

                                    ++component;
                                    destination += 2;

                                    break;
                                  }

                                case DataComponentInterpretation::
                                  component_is_part_of_vector:
                                  {
                                    // A vector component. Put the
                                    // DoFHandlerType::space_dimension
                                    // components into the next set of
                                    // contiguous rows
                                    // ('destination+c'); then move 'component'
                                    // forward by spacedim (so we get to the
                                    // next component after the current vector)
                                    // and 'destination' forward by two*spacedim
                                    // (because we're going to put the imaginary
                                    // part of the vector into the subsequent
                                    // spacedim slots).
                                    const unsigned int size =
                                      DoFHandlerType::space_dimension;
                                    for (unsigned int c = 0; c < size; ++c)
                                      for (unsigned int q = 0; q < n_q_points;
                                           ++q)
                                        patch.data(destination + c, q) =
                                          scratch_data.patch_values_system
                                            .solution_values[q](component + c);

                                    component += size;
                                    destination += 2 * size;

                                    break;
                                  }

                                case DataComponentInterpretation::
                                  component_is_part_of_tensor:
                                  {
                                    // Same approach as for vectors above.
                                    const unsigned int size =
                                      DoFHandlerType::space_dimension *
                                      DoFHandlerType::space_dimension;
                                    for (unsigned int c = 0; c < size; ++c)
                                      for (unsigned int q = 0; q < n_q_points;
                                           ++q)
                                        patch.data(destination + c, q) =
                                          scratch_data.patch_values_system
                                            .solution_values[q](component + c);

                                    component += size;
                                    destination += 2 * size;

                                    break;
                                  }

                                default:
                                  Assert(false, ExcNotImplemented());
                              }
                          }
                      }

                      // And now we need to do the same thing again for the
                      // imaginary parts, starting at the top of the list of
                      // components/destinations again.
                      dataset->get_function_values(
                        this_fe_patch_values,
                        internal::DataOutImplementation::ComponentExtractor::
                          imaginary_part,
                        scratch_data.patch_values_system.solution_values);
                      {
                        unsigned int destination = offset;
                        for (unsigned int component = 0;
                             component < n_components;
                             /* component is updated below */)
                          {
                            switch (
                              dataset->data_component_interpretation[component])
                              {
                                case DataComponentInterpretation::
                                  component_is_scalar:
                                  {
                                    // OK, a scalar component. Put all of the
                                    // values into the row past the current one
                                    // ('destination+1') since 'destination' is
                                    // occupied by the real part.
                                    for (unsigned int q = 0; q < n_q_points;
                                         ++q)
                                      patch.data(destination + 1, q) =
                                        scratch_data.patch_values_system
                                          .solution_values[q](component);

                                    ++component;
                                    destination += 2;

                                    break;
                                  }

                                case DataComponentInterpretation::
                                  component_is_part_of_vector:
                                  {
                                    // A vector component. Put the
                                    // DoFHandlerType::space_dimension
                                    // components into the set of contiguous
                                    // rows that follow the real parts
                                    // ('destination+spacedim+c').
                                    const unsigned int size =
                                      DoFHandlerType::space_dimension;
                                    for (unsigned int c = 0; c < size; ++c)
                                      for (unsigned int q = 0; q < n_q_points;
                                           ++q)
                                        patch.data(destination + size + c, q) =
                                          scratch_data.patch_values_system
                                            .solution_values[q](component + c);

                                    component += size;
                                    destination += 2 * size;

                                    break;
                                  }

                                case DataComponentInterpretation::
                                  component_is_part_of_tensor:
                                  {
                                    // Same as for vectors.
                                    const unsigned int size =
                                      DoFHandlerType::space_dimension *
                                      DoFHandlerType::space_dimension;
                                    for (unsigned int c = 0; c < size; ++c)
                                      for (unsigned int q = 0; q < n_q_points;
                                           ++q)
                                        patch.data(destination + size + c, q) =
                                          scratch_data.patch_values_system
                                            .solution_values[q](component + c);

                                    component += size;
                                    destination += 2 * size;

                                    break;
                                  }

                                default:
                                  Assert(false, ExcNotImplemented());
                              }
                          }
                      }

                      // Increment the counter for the actual data record. We
                      // need to move it forward a number of positions equal to
                      // the number of components of this data set, times two
                      // because we dealt with a complex-valued input vector
                      offset += dataset->n_output_variables * 2;
                    }
                }
            }

          // Also update the dataset_number index that we carry along with the
          // for-loop over all data sets.
          ++dataset_number;
        }

      // Then do the cell data. At least, we don't have to worry about
      // complex-valued vectors/tensors since cell data is always scalar.
      if (this->cell_data.size() != 0)
        {
          Assert(!cell_and_index->first->has_children(), ExcNotImplemented());

          for (const auto &dataset : this->cell_data)
            {
              // as above, first output the real part
              {
                const double value =
                  dataset->get_cell_data_value(cell_and_index->second,
                                               internal::DataOutImplementation::
                                                 ComponentExtractor::real_part);
                for (unsigned int q = 0; q < n_q_points; ++q)
                  patch.data(offset, q) = value;
              }

              // and if there is one, also output the imaginary part
              if (dataset->is_complex_valued() == true)
                {
                  const double value = dataset->get_cell_data_value(
                    cell_and_index->second,
                    internal::DataOutImplementation::ComponentExtractor::
                      imaginary_part);
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    patch.data(offset + 1, q) = value;
                }

              offset += (dataset->is_complex_valued() ? 2 : 1);
            }
        }
    }


  for (const unsigned int f :
       GeometryInfo<DoFHandlerType::dimension>::face_indices())
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
      if (cell_and_index->first->at_boundary(f) ||
          (cell_and_index->first->neighbor(f)->level() !=
           cell_and_index->first->level()))
        {
          patch.neighbors[f] = numbers::invalid_unsigned_int;
          continue;
        }

      const cell_iterator neighbor = cell_and_index->first->neighbor(f);
      Assert(static_cast<unsigned int>(neighbor->level()) <
               scratch_data.cell_to_patch_index_map->size(),
             ExcInternalError());
      if ((static_cast<unsigned int>(neighbor->index()) >=
           (*scratch_data.cell_to_patch_index_map)[neighbor->level()].size()) ||
          ((*scratch_data.cell_to_patch_index_map)[neighbor->level()]
                                                  [neighbor->index()] ==
           dealii::DataOutBase::Patch<DoFHandlerType::dimension>::no_neighbor))
        {
          patch.neighbors[f] = numbers::invalid_unsigned_int;
          continue;
        }

      // now, there is a neighbor, so get its patch number and set it for the
      // neighbor index
      patch.neighbors[f] =
        (*scratch_data
            .cell_to_patch_index_map)[neighbor->level()][neighbor->index()];
    }

  const unsigned int patch_idx =
    (*scratch_data.cell_to_patch_index_map)[cell_and_index->first->level()]
                                           [cell_and_index->first->index()];
  // did we mess up the indices?
  Assert(patch_idx < this->patches.size(), ExcInternalError());
  patch.patch_index = patch_idx;

  // Put the patch into the patches vector. instead of copying the data,
  // simply swap the contents to avoid the penalty of writing into another
  // processor's memory
  this->patches[patch_idx].swap(patch);
}



template <int dim, typename DoFHandlerType>
void
DataOut<dim, DoFHandlerType>::build_patches(const unsigned int n_subdivisions)
{
  build_patches(StaticMappingQ1<DoFHandlerType::dimension,
                                DoFHandlerType::space_dimension>::mapping,
                n_subdivisions,
                no_curved_cells);
}



template <int dim, typename DoFHandlerType>
void
DataOut<dim, DoFHandlerType>::build_patches(
  const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
    &                    mapping,
  const unsigned int     n_subdivisions_,
  const CurvedCellRegion curved_region)
{
  // Check consistency of redundant template parameter
  Assert(dim == DoFHandlerType::dimension,
         ExcDimensionMismatch(dim, DoFHandlerType::dimension));

  Assert(this->triangulation != nullptr,
         Exceptions::DataOutImplementation::ExcNoTriangulationSelected());

  const unsigned int n_subdivisions =
    (n_subdivisions_ != 0) ? n_subdivisions_ : this->default_subdivisions;
  Assert(n_subdivisions >= 1,
         Exceptions::DataOutImplementation::ExcInvalidNumberOfSubdivisions(
           n_subdivisions));

  this->validate_dataset_names();

  // First count the cells we want to create patches of. Also fill the object
  // that maps the cell indices to the patch numbers, as this will be needed
  // for generation of neighborship information.
  // Note, there is a confusing mess of different indices here at play:
  // - patch_index: the index of a patch in all_cells
  // - cell->index: only unique on each level, used in cell_to_patch_index_map
  // - active_index: index for a cell when counting from begin_active() using
  //   ++cell (identical to cell->active_cell_index())
  // - cell_index: unique index of a cell counted using
  //   next_cell_function() starting from first_cell_function()
  //
  // It turns out that we create one patch for each selected cell, so
  // patch_index==cell_index.
  //
  // Now construct the map such that
  // cell_to_patch_index_map[cell->level][cell->index] = patch_index
  std::vector<std::vector<unsigned int>> cell_to_patch_index_map;
  cell_to_patch_index_map.resize(this->triangulation->n_levels());
  for (unsigned int l = 0; l < this->triangulation->n_levels(); ++l)
    {
      // max_index is the largest cell->index on level l
      unsigned int max_index = 0;
      for (cell_iterator cell = first_cell_function(*this->triangulation);
           cell != this->triangulation->end();
           cell = next_cell_function(*this->triangulation, cell))
        if (static_cast<unsigned int>(cell->level()) == l)
          max_index =
            std::max(max_index, static_cast<unsigned int>(cell->index()));

      cell_to_patch_index_map[l].resize(
        max_index + 1,
        dealii::DataOutBase::Patch<
          DoFHandlerType::dimension,
          DoFHandlerType::space_dimension>::no_neighbor);
    }

  // will be all_cells[patch_index] = pair(cell, active_index)
  std::vector<std::pair<cell_iterator, unsigned int>> all_cells;
  {
    // important: we need to compute the active_index of the cell in the range
    // 0..n_active_cells() because this is where we need to look up cell
    // data from (cell data vectors do not have the length distance computed by
    // first_cell_function/next_cell_function because this might skip
    // some values (FilteredIterator).
    auto          active_cell  = this->triangulation->begin_active();
    unsigned int  active_index = 0;
    cell_iterator cell         = first_cell_function(*this->triangulation);
    for (; cell != this->triangulation->end();
         cell = next_cell_function(*this->triangulation, cell))
      {
        // move forward until active_cell points at the cell (cell) we are
        // looking at to compute the current active_index
        while (active_cell != this->triangulation->end() && cell->is_active() &&
               decltype(active_cell)(cell) != active_cell)
          {
            ++active_cell;
            ++active_index;
          }

        Assert(static_cast<unsigned int>(cell->level()) <
                 cell_to_patch_index_map.size(),
               ExcInternalError());
        Assert(static_cast<unsigned int>(cell->index()) <
                 cell_to_patch_index_map[cell->level()].size(),
               ExcInternalError());
        Assert(active_index < this->triangulation->n_active_cells(),
               ExcInternalError());
        cell_to_patch_index_map[cell->level()][cell->index()] =
          all_cells.size();

        all_cells.emplace_back(cell, active_index);
      }
  }

  this->patches.clear();
  this->patches.resize(all_cells.size());

  // Now create a default object for the WorkStream object to work with. The
  // first step is to count how many output data sets there will be. This is,
  // in principle, just the number of components of each data set, but we
  // need to allocate two entries per component if there are
  // complex-valued input data (unless we use a postprocessor on this
  // output -- all postprocessor outputs are real-valued)
  unsigned int n_datasets = 0;
  for (unsigned int i = 0; i < this->cell_data.size(); ++i)
    n_datasets += (this->cell_data[i]->is_complex_valued() &&
                       (this->cell_data[i]->postprocessor == nullptr) ?
                     2 :
                     1);
  for (unsigned int i = 0; i < this->dof_data.size(); ++i)
    n_datasets += (this->dof_data[i]->n_output_variables *
                   (this->dof_data[i]->is_complex_valued() &&
                        (this->dof_data[i]->postprocessor == nullptr) ?
                      2 :
                      1));

  std::vector<unsigned int> n_postprocessor_outputs(this->dof_data.size());
  for (unsigned int dataset = 0; dataset < this->dof_data.size(); ++dataset)
    if (this->dof_data[dataset]->postprocessor)
      n_postprocessor_outputs[dataset] =
        this->dof_data[dataset]->n_output_variables;
    else
      n_postprocessor_outputs[dataset] = 0;

  const CurvedCellRegion curved_cell_region =
    (n_subdivisions < 2 ? no_curved_cells : curved_region);

  UpdateFlags update_flags = update_values;
  if (curved_cell_region != no_curved_cells)
    update_flags |= update_quadrature_points;

  for (unsigned int i = 0; i < this->dof_data.size(); ++i)
    if (this->dof_data[i]->postprocessor)
      update_flags |=
        this->dof_data[i]->postprocessor->get_needed_update_flags();
  // perhaps update_normal_vectors is present, which would only be useful on
  // faces, but we may not use it here.
  Assert(
    !(update_flags & update_normal_vectors),
    ExcMessage(
      "The update of normal vectors may not be requested for evaluation of "
      "data on cells via DataPostprocessor."));

  internal::DataOutImplementation::ParallelData<DoFHandlerType::dimension,
                                                DoFHandlerType::space_dimension>
    thread_data(n_datasets,
                n_subdivisions,
                n_postprocessor_outputs,
                mapping,
                this->get_fes(),
                update_flags,
                cell_to_patch_index_map);

  auto worker = [this, n_subdivisions, curved_cell_region](
                  const std::pair<cell_iterator, unsigned int> *cell_and_index,
                  internal::DataOutImplementation::ParallelData<
                    DoFHandlerType::dimension,
                    DoFHandlerType::space_dimension> &scratch_data,
                  // this function doesn't actually need a copy data object --
                  // it just writes everything right into the output array
                  int) {
    this->build_one_patch(cell_and_index,
                          scratch_data,
                          n_subdivisions,
                          curved_cell_region);
  };

  // now build the patches in parallel
  if (all_cells.size() > 0)
    WorkStream::run(all_cells.data(),
                    all_cells.data() + all_cells.size(),
                    worker,
                    // no copy-local-to-global function needed here
                    std::function<void(const int)>(),
                    thread_data,
                    /* dummy CopyData object = */ 0,
                    // experimenting shows that we can make things run a bit
                    // faster if we increase the number of cells we work on
                    // per item (i.e., WorkStream's chunk_size argument,
                    // about 10% improvement) and the items in flight at any
                    // given time (another 5% on the testcase discussed in
                    // @ref workstream_paper, on 32 cores) and if
                    8 * MultithreadInfo::n_threads(),
                    64);
}



template <int dim, typename DoFHandlerType>
void
DataOut<dim, DoFHandlerType>::set_cell_selection(
  const std::function<cell_iterator(const Triangulation<dim, spacedim> &)>
    &                                                        first_cell,
  const std::function<cell_iterator(const Triangulation<dim, spacedim> &,
                                    const cell_iterator &)> &next_cell)
{
  first_cell_function = first_cell;
  next_cell_function  = next_cell;
}



template <int dim, typename DoFHandlerType>
void
DataOut<dim, DoFHandlerType>::set_cell_selection(
  const FilteredIterator<cell_iterator> &filtered_iterator)
{
  const auto first_cell =
    [filtered_iterator](const Triangulation<dim, spacedim> &triangulation) {
      // Create a copy of the filtered iterator so that we can
      // call a non-const function -- though we are really only
      // interested in the return value of that function, not the
      // state of the object
      FilteredIterator<cell_iterator> x = filtered_iterator;
      return x.set_to_next_positive(triangulation.begin());
    };


  const auto next_cell =
    [filtered_iterator](const Triangulation<dim, spacedim> &,
                        const cell_iterator &cell) {
      // Create a copy of the filtered iterator so that we can
      // call a non-const function -- though we are really only
      // interested in the return value of that function, not the
      // state of the object
      FilteredIterator<cell_iterator> x = filtered_iterator;

      // Set the iterator to 'cell'. Since 'cell' must satisfy the
      // predicate (that's how it was created), set_to_next_positive
      // simply sets the iterator to 'cell'.
      x.set_to_next_positive(cell);

      // Advance by one:
      ++x;

      return x;
    };

  set_cell_selection(first_cell, next_cell);
}



template <int dim, typename DoFHandlerType>
const std::pair<typename DataOut<dim, DoFHandlerType>::FirstCellFunctionType,
                typename DataOut<dim, DoFHandlerType>::NextCellFunctionType>
DataOut<dim, DoFHandlerType>::get_cell_selection() const
{
  return std::make_pair(first_cell_function, next_cell_function);
}



template <int dim, typename DoFHandlerType>
typename DataOut<dim, DoFHandlerType>::cell_iterator
DataOut<dim, DoFHandlerType>::first_cell()
{
  return this->triangulation->begin_active();
}



template <int dim, typename DoFHandlerType>
typename DataOut<dim, DoFHandlerType>::cell_iterator
DataOut<dim, DoFHandlerType>::next_cell(
  const typename DataOut<dim, DoFHandlerType>::cell_iterator &cell)
{
  // convert the iterator to an active_iterator and advance this to the next
  // active cell
  typename Triangulation<DoFHandlerType::dimension,
                         DoFHandlerType::space_dimension>::active_cell_iterator
    active_cell = cell;
  ++active_cell;
  return active_cell;
}



template <int dim, typename DoFHandlerType>
typename DataOut<dim, DoFHandlerType>::cell_iterator
DataOut<dim, DoFHandlerType>::first_locally_owned_cell()
{
  typename DataOut<dim, DoFHandlerType>::cell_iterator cell = first_cell();

  // skip cells if the current one has no children (is active) and is a ghost
  // or artificial cell
  while ((cell != this->triangulation->end()) && cell->is_active() &&
         !cell->is_locally_owned())
    cell = next_cell(cell);

  return cell;
}



template <int dim, typename DoFHandlerType>
typename DataOut<dim, DoFHandlerType>::cell_iterator
DataOut<dim, DoFHandlerType>::next_locally_owned_cell(
  const typename DataOut<dim, DoFHandlerType>::cell_iterator &old_cell)
{
  typename DataOut<dim, DoFHandlerType>::cell_iterator cell =
    next_cell(old_cell);
  while ((cell != this->triangulation->end()) && cell->is_active() &&
         !cell->is_locally_owned())
    cell = next_cell(cell);
  return cell;
}


// explicit instantiations
#include "data_out.inst"

DEAL_II_NAMESPACE_CLOSE
