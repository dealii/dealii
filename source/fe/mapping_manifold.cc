// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q1.h>

#include <cmath>
#include <algorithm>
#include <numeric>


DEAL_II_NAMESPACE_OPEN

template<int dim, int spacedim>
MappingManifold<dim,spacedim>::InternalData::InternalData ()
{}



template<int dim, int spacedim>
std::size_t
MappingManifold<dim,spacedim>::InternalData::memory_consumption () const
{
  return (Mapping<dim,spacedim>::InternalDataBase::memory_consumption() );
  // MemoryConsumption::memory_consumption (shape_values) +
  // MemoryConsumption::memory_consumption (shape_derivatives) +
  // MemoryConsumption::memory_consumption (covariant) +
  // MemoryConsumption::memory_consumption (contravariant) +
  // MemoryConsumption::memory_consumption (unit_tangentials) +
  // MemoryConsumption::memory_consumption (aux) +
  // MemoryConsumption::memory_consumption (mapping_support_points) +
  // MemoryConsumption::memory_consumption (cell_of_current_support_points) +
  // MemoryConsumption::memory_consumption (volume_elements) +
  // MemoryConsumption::memory_consumption (polynomial_degree) +
  // MemoryConsumption::memory_consumption (n_shape_functions));
}


template <int dim, int spacedim>
void
MappingManifold<dim,spacedim>::InternalData::
initialize (const UpdateFlags      update_flags,
            const Quadrature<dim> &q,
            const unsigned int     n_original_q_points)
{
  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values()
  this->update_each = update_flags;

  const unsigned int n_q_points = q.size();

  // Update the weights used in the Manifold Quadrature formulas
  compute_manifold_quadrature_weights(q);

  // see if we need the (transformation) shape function values
  // and/or gradients and resize the necessary arrays
  if (this->update_each & update_quadrature_points)
    cell_manifold_quadratures.resize(q.size());

  // if (this->update_each & (update_covariant_transformation
  //                          | update_contravariant_transformation
  //                          | update_JxW_values
  //                          | update_boundary_forms
  //                          | update_normal_vectors
  //                          | update_jacobians
  //                          | update_jacobian_grads
  //                          | update_inverse_jacobians
  //                          | update_jacobian_pushed_forward_grads
  //                          | update_jacobian_2nd_derivatives
  //                          | update_jacobian_pushed_forward_2nd_derivatives
  //                          | update_jacobian_3rd_derivatives
  //                          | update_jacobian_pushed_forward_3rd_derivatives))
  //   shape_derivatives.resize(n_shape_functions * n_q_points);

  // if (this->update_each & update_covariant_transformation)
  //   covariant.resize(n_original_q_points);

  // if (this->update_each & update_contravariant_transformation)
  //   contravariant.resize(n_original_q_points);

  // if (this->update_each & update_volume_elements)
  //   volume_elements.resize(n_original_q_points);

  // if (this->update_each &
  //     (update_jacobian_grads | update_jacobian_pushed_forward_grads) )
  //   shape_second_derivatives.resize(n_shape_functions * n_q_points);

  // if (this->update_each &
  //     (update_jacobian_2nd_derivatives | update_jacobian_pushed_forward_2nd_derivatives) )
  //   shape_third_derivatives.resize(n_shape_functions * n_q_points);

  // if (this->update_each &
  //     (update_jacobian_3rd_derivatives | update_jacobian_pushed_forward_3rd_derivatives) )
  //   shape_fourth_derivatives.resize(n_shape_functions * n_q_points);

  // // now also fill the various fields with their correct values
  // compute_shape_function_values (q.get_points());
}



template <int dim, int spacedim>
void
MappingManifold<dim,spacedim>::InternalData::
initialize_face (const UpdateFlags      update_flags,
                 const Quadrature<dim> &q,
                 const unsigned int     n_original_q_points)
{
  // initialize (update_flags, q, n_original_q_points);

  // if (dim > 1)
  //   {
  //     if (this->update_each & update_boundary_forms)
  //       {
  //         aux.resize (dim-1, std::vector<Tensor<1,spacedim> > (n_original_q_points));

  //         // Compute tangentials to the
  //         // unit cell.
  //         const unsigned int nfaces = GeometryInfo<dim>::faces_per_cell;
  //         unit_tangentials.resize (nfaces*(dim-1),
  //                                  std::vector<Tensor<1,dim> > (n_original_q_points));
  //         if (dim==2)
  //           {
  //             // ensure a counterclockwise
  //             // orientation of tangentials
  //             static const int tangential_orientation[4]= {-1,1,1,-1};
  //             for (unsigned int i=0; i<nfaces; ++i)
  //               {
  //                 Tensor<1,dim> tang;
  //                 tang[1-i/2]=tangential_orientation[i];
  //                 std::fill (unit_tangentials[i].begin(),
  //                            unit_tangentials[i].end(), tang);
  //               }
  //           }
  //         else if (dim==3)
  //           {
  //             for (unsigned int i=0; i<nfaces; ++i)
  //               {
  //                 Tensor<1,dim> tang1, tang2;

  //                 const unsigned int nd=
  //                   GeometryInfo<dim>::unit_normal_direction[i];

  //                 // first tangential
  //                 // vector in direction
  //                 // of the (nd+1)%3 axis
  //                 // and inverted in case
  //                 // of unit inward normal
  //                 tang1[(nd+1)%dim]=GeometryInfo<dim>::unit_normal_orientation[i];
  //                 // second tangential
  //                 // vector in direction
  //                 // of the (nd+2)%3 axis
  //                 tang2[(nd+2)%dim]=1.;

  //                 // same unit tangents
  //                 // for all quadrature
  //                 // points on this face
  //                 std::fill (unit_tangentials[i].begin(),
  //                            unit_tangentials[i].end(), tang1);
  //                 std::fill (unit_tangentials[nfaces+i].begin(),
  //                            unit_tangentials[nfaces+i].end(), tang2);
  //               }
  //           }
  //       }
  //   }
}


template<int dim, int spacedim>
MappingManifold<dim,spacedim>::MappingManifold ()
  :
  fe_q(1)
  // support_point_weights_on_quad (compute_support_point_weights_on_quad<dim>(this->polynomial_degree)),
  // support_point_weights_on_hex (compute_support_point_weights_on_hex<dim>(this->polynomial_degree)),
{
}



template<int dim, int spacedim>
MappingManifold<dim,spacedim>::MappingManifold (const MappingManifold<dim,spacedim> &mapping)
  :
  fe_q(1)
{}




template<int dim, int spacedim>
Mapping<dim,spacedim> *
MappingManifold<dim,spacedim>::clone () const
{
  return new MappingManifold<dim,spacedim>(*this);
}



template<int dim, int spacedim>
Point<spacedim>
MappingManifold<dim,spacedim>::
transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<dim> &p) const
{
  std::vector<Point<spacedim> > vertices;
  std::vector<Point<spacedim> > weights;
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
    {
      vertices.push_back(cell->vertex(v));
      weights.push_back(fe_q.shape_value(v,p));
    }
  return cell->get_manifold().get_new_point(Quadrature<spacedim>(vertices, weights));
}


// In the code below, GCC tries to instantiate MappingManifold<3,4> when
// seeing which of the overloaded versions of
// do_transform_real_to_unit_cell_internal() to call. This leads to bad
// error messages and, generally, nothing very good. Avoid this by ensuring
// that this class exists, but does not have an inner InternalData
// type, thereby ruling out the codim-1 version of the function
// below when doing overload resolution.
template <>
class MappingManifold<3,4>
{};


template<int dim, int spacedim>
UpdateFlags
MappingManifold<dim,spacedim>::requires_update_flags (const UpdateFlags in) const
{
  // add flags if the respective quantities are necessary to compute
  // what we need. note that some flags appear in both the conditions
  // and in subsequent set operations. this leads to some circular
  // logic. the only way to treat this is to iterate. since there are
  // 5 if-clauses in the loop, it will take at most 5 iterations to
  // converge. do them:
  UpdateFlags out = in;
  for (unsigned int i=0; i<5; ++i)
    {
      // The following is a little incorrect:
      // If not applied on a face,
      // update_boundary_forms does not
      // make sense. On the other hand,
      // it is necessary on a
      // face. Currently,
      // update_boundary_forms is simply
      // ignored for the interior of a
      // cell.
      if (out & (update_JxW_values
                 | update_normal_vectors))
        out |= update_boundary_forms;

      if (out & (update_covariant_transformation
                 | update_JxW_values
                 | update_jacobians
                 | update_jacobian_grads
                 | update_boundary_forms
                 | update_normal_vectors))
        out |= update_contravariant_transformation;

      if (out & (update_inverse_jacobians
                 | update_jacobian_pushed_forward_grads
                 | update_jacobian_pushed_forward_2nd_derivatives
                 | update_jacobian_pushed_forward_3rd_derivatives) )
        out |= update_covariant_transformation;

      // The contravariant transformation
      // used in the Piola transformation, which
      // requires the determinant of the
      // Jacobi matrix of the transformation.
      // Because we have no way of knowing here whether the finite
      // elements wants to use the contravariant of the Piola
      // transforms, we add the JxW values to the list of flags to be
      // updated for each cell.
      if (out & update_contravariant_transformation)
        out |= update_JxW_values;

      if (out & update_normal_vectors)
        out |= update_JxW_values;
    }

  return out;
}



template<int dim, int spacedim>
typename MappingManifold<dim,spacedim>::InternalData *
MappingManifold<dim,spacedim>::get_data (const UpdateFlags update_flags,
                                         const Quadrature<dim> &q) const
{
  InternalData *data = new InternalData();
  data->initialize (this->requires_update_flags(update_flags), q, q.size());

  return data;
}



template<int dim, int spacedim>
typename MappingManifold<dim,spacedim>::InternalData *
MappingManifold<dim,spacedim>::get_face_data (const UpdateFlags        update_flags,
                                              const Quadrature<dim-1> &quadrature) const
{
  InternalData *data = new InternalData();
  data->initialize_face (this->requires_update_flags(update_flags),
                         QProjector<dim>::project_to_all_faces(quadrature),
                         quadrature.size());

  return data;
}



template<int dim, int spacedim>
typename MappingManifold<dim,spacedim>::InternalData *
MappingManifold<dim,spacedim>::get_subface_data (const UpdateFlags update_flags,
                                                 const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData();
  data->initialize_face (this->requires_update_flags(update_flags),
                         QProjector<dim>::project_to_all_subfaces(quadrature),
                         quadrature.size());

  return data;
}



namespace internal
{
  namespace
  {
    /**
     * Compute the locations of quadrature points on the object described by
     * the first argument (and the cell for which the mapping support points
     * have already been set), but only if the update_flags of the @p data
     * argument indicate so.
     */
    template <int dim, int spacedim>
    void
    maybe_compute_q_points (const typename QProjector<dim>::DataSetDescriptor                 data_set,
                            const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data,
                            std::vector<Point<spacedim> >                                     &quadrature_points)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_quadrature_points)
        {
          for (unsigned int point=0; point<quadrature_points.size(); ++point)
            {
              const double *shape = &data.shape(point+data_set,0);
              Point<spacedim> result = (shape[0] *
                                        data.mapping_support_points[0]);
              for (unsigned int k=1; k<data.n_shape_functions; ++k)
                for (unsigned int i=0; i<spacedim; ++i)
                  result[i] += shape[k] * data.mapping_support_points[k][i];
              quadrature_points[point] = result;
            }
        }
    }


    /**
     * Update the co- and contravariant matrices as well as their determinant, for the cell
     * described stored in the data object, but only if the update_flags of the @p data
     * argument indicate so.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_Jacobians (const CellSimilarity::Similarity                                   cell_similarity,
                            const typename dealii::QProjector<dim>::DataSetDescriptor          data_set,
                            const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_contravariant_transformation)
        // if the current cell is just a
        // translation of the previous one, no
        // need to recompute jacobians...
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();

            std::fill(data.contravariant.begin(), data.contravariant.end(),
                      DerivativeForm<1,dim,spacedim>());

            Assert (data.n_shape_functions > 0, ExcInternalError());
            const Tensor<1,spacedim> *supp_pts =
              &data.mapping_support_points[0];

            for (unsigned int point=0; point<n_q_points; ++point)
              {
                const Tensor<1,dim> *data_derv =
                  &data.derivative(point+data_set, 0);

                double result [spacedim][dim];

                // peel away part of sum to avoid zeroing the
                // entries and adding for the first time
                for (unsigned int i=0; i<spacedim; ++i)
                  for (unsigned int j=0; j<dim; ++j)
                    result[i][j] = data_derv[0][j] * supp_pts[0][i];
                for (unsigned int k=1; k<data.n_shape_functions; ++k)
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      result[i][j] += data_derv[k][j] * supp_pts[k][i];

                // write result into contravariant data. for
                // j=dim in the case dim<spacedim, there will
                // never be any nonzero data that arrives in
                // here, so it is ok anyway because it was
                // initialized to zero at the initialization
                for (unsigned int i=0; i<spacedim; ++i)
                  for (unsigned int j=0; j<dim; ++j)
                    data.contravariant[point][i][j] = result[i][j];
              }
          }

      if (update_flags & update_covariant_transformation)
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();
            for (unsigned int point=0; point<n_q_points; ++point)
              {
                data.covariant[point] = (data.contravariant[point]).covariant_form();
              }
          }

      if (update_flags & update_volume_elements)
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();
            for (unsigned int point=0; point<n_q_points; ++point)
              data.volume_elements[point] = data.contravariant[point].determinant();
          }

    }

    /**
     * Update the Hessian of the transformation from unit to real cell, the
     * Jacobian gradients.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_grads (const CellSimilarity::Similarity                                   cell_similarity,
                                 const typename QProjector<dim>::DataSetDescriptor                  data_set,
                                 const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data,
                                 std::vector<DerivativeForm<2,dim,spacedim> >                      &jacobian_grads)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_grads)
        {
          const unsigned int n_q_points = jacobian_grads.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const Tensor<2,dim> *second =
                    &data.second_derivative(point+data_set, 0);
                  double result [spacedim][dim][dim];
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        result[i][j][l] = (second[0][j][l] *
                                           data.mapping_support_points[0][i]);
                  for (unsigned int k=1; k<data.n_shape_functions; ++k)
                    for (unsigned int i=0; i<spacedim; ++i)
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          result[i][j][l]
                          += (second[k][j][l]
                              *
                              data.mapping_support_points[k][i]);

                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        jacobian_grads[point][i][j][l] = result[i][j][l];
                }
            }
        }
    }

    /**
     * Update the Hessian of the transformation from unit to real cell, the
     * Jacobian gradients, pushed forward to the real cell coordinates.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_pushed_forward_grads (const CellSimilarity::Similarity                                   cell_similarity,
                                                const typename QProjector<dim>::DataSetDescriptor                  data_set,
                                                const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data,
                                                std::vector<Tensor<3,spacedim> >                      &jacobian_pushed_forward_grads)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_pushed_forward_grads)
        {
          const unsigned int n_q_points = jacobian_pushed_forward_grads.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              double tmp[spacedim][spacedim][spacedim];
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const Tensor<2,dim> *second =
                    &data.second_derivative(point+data_set, 0);
                  double result [spacedim][dim][dim];
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        result[i][j][l] = (second[0][j][l] *
                                           data.mapping_support_points[0][i]);
                  for (unsigned int k=1; k<data.n_shape_functions; ++k)
                    for (unsigned int i=0; i<spacedim; ++i)
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          result[i][j][l]
                          += (second[k][j][l]
                              *
                              data.mapping_support_points[k][i]);

                  // first push forward the j-components
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        {
                          tmp[i][j][l] = result[i][0][l] *
                                         data.covariant[point][j][0];
                          for (unsigned int jr=1; jr<dim; ++jr)
                            {
                              tmp[i][j][l] += result[i][jr][l] *
                                              data.covariant[point][j][jr];
                            }
                        }

                  // now, pushing forward the l-components
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<spacedim; ++l)
                        {
                          jacobian_pushed_forward_grads[point][i][j][l] = tmp[i][j][0] *
                                                                          data.covariant[point][l][0];
                          for (unsigned int lr=1; lr<dim; ++lr)
                            {
                              jacobian_pushed_forward_grads[point][i][j][l] += tmp[i][j][lr] *
                                                                               data.covariant[point][l][lr];
                            }

                        }
                }
            }
        }
    }

    /**
     * Update the third derivatives of the transformation from unit to real cell, the
     * Jacobian hessians.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_2nd_derivatives (const CellSimilarity::Similarity                              cell_similarity,
                                           const typename QProjector<dim>::DataSetDescriptor             data_set,
                                           const typename dealii::MappingManifold<dim,spacedim>::InternalData &data,
                                           std::vector<DerivativeForm<3,dim,spacedim> >                 &jacobian_2nd_derivatives)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_2nd_derivatives)
        {
          const unsigned int n_q_points = jacobian_2nd_derivatives.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const Tensor<3,dim> *third =
                    &data.third_derivative(point+data_set, 0);
                  double result [spacedim][dim][dim][dim];
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          result[i][j][l][m] = (third[0][j][l][m] *
                                                data.mapping_support_points[0][i]);
                  for (unsigned int k=1; k<data.n_shape_functions; ++k)
                    for (unsigned int i=0; i<spacedim; ++i)
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          for (unsigned int m=0; m<dim; ++m)
                            result[i][j][l][m]
                            += (third[k][j][l][m]
                                *
                                data.mapping_support_points[k][i]);

                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          jacobian_2nd_derivatives[point][i][j][l][m] = result[i][j][l][m];
                }
            }
        }
    }

    /**
     * Update the Hessian of the Hessian of the transformation from unit
     * to real cell, the Jacobian Hessian gradients, pushed forward to the
     * real cell coordinates.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_pushed_forward_2nd_derivatives (const CellSimilarity::Similarity                                   cell_similarity,
                                                          const typename QProjector<dim>::DataSetDescriptor                  data_set,
                                                          const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data,
                                                          std::vector<Tensor<4,spacedim> >                      &jacobian_pushed_forward_2nd_derivatives)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_pushed_forward_2nd_derivatives)
        {
          const unsigned int n_q_points = jacobian_pushed_forward_2nd_derivatives.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              double tmp[spacedim][spacedim][spacedim][spacedim];
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const Tensor<3,dim> *third =
                    &data.third_derivative(point+data_set, 0);
                  double result [spacedim][dim][dim][dim];
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          result[i][j][l][m] = (third[0][j][l][m] *
                                                data.mapping_support_points[0][i]);
                  for (unsigned int k=1; k<data.n_shape_functions; ++k)
                    for (unsigned int i=0; i<spacedim; ++i)
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          for (unsigned int m=0; m<dim; ++m)
                            result[i][j][l][m]
                            += (third[k][j][l][m]
                                *
                                data.mapping_support_points[k][i]);

                  // push forward the j-coordinate
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          {
                            jacobian_pushed_forward_2nd_derivatives[point][i][j][l][m]
                              = result[i][0][l][m]*
                                data.covariant[point][j][0];
                            for (unsigned int jr=1; jr<dim; ++jr)
                              jacobian_pushed_forward_2nd_derivatives[point][i][j][l][m]
                              += result[i][jr][l][m]*
                                 data.covariant[point][j][jr];
                          }

                  // push forward the l-coordinate
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<spacedim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          {
                            tmp[i][j][l][m]
                              = jacobian_pushed_forward_2nd_derivatives[point][i][j][0][m]*
                                data.covariant[point][l][0];
                            for (unsigned int lr=1; lr<dim; ++lr)
                              tmp[i][j][l][m]
                              += jacobian_pushed_forward_2nd_derivatives[point][i][j][lr][m]*
                                 data.covariant[point][l][lr];
                          }

                  // push forward the m-coordinate
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<spacedim; ++l)
                        for (unsigned int m=0; m<spacedim; ++m)
                          {
                            jacobian_pushed_forward_2nd_derivatives[point][i][j][l][m]
                              = tmp[i][j][l][0]*
                                data.covariant[point][m][0];
                            for (unsigned int mr=1; mr<dim; ++mr)
                              jacobian_pushed_forward_2nd_derivatives[point][i][j][l][m]
                              += tmp[i][j][l][mr]*
                                 data.covariant[point][m][mr];
                          }
                }
            }
        }
    }

    /**
         * Update the fourth derivatives of the transformation from unit to real cell, the
         * Jacobian hessian gradients.
         *
         * Skip the computation if possible as indicated by the first argument.
         */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_3rd_derivatives (const CellSimilarity::Similarity                              cell_similarity,
                                           const typename QProjector<dim>::DataSetDescriptor             data_set,
                                           const typename dealii::MappingManifold<dim,spacedim>::InternalData &data,
                                           std::vector<DerivativeForm<4,dim,spacedim> >                 &jacobian_3rd_derivatives)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_3rd_derivatives)
        {
          const unsigned int n_q_points = jacobian_3rd_derivatives.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const Tensor<4,dim> *fourth =
                    &data.fourth_derivative(point+data_set, 0);
                  double result [spacedim][dim][dim][dim][dim];
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          for (unsigned int n=0; n<dim; ++n)
                            result[i][j][l][m][n] = (fourth[0][j][l][m][n] *
                                                     data.mapping_support_points[0][i]);
                  for (unsigned int k=1; k<data.n_shape_functions; ++k)
                    for (unsigned int i=0; i<spacedim; ++i)
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          for (unsigned int m=0; m<dim; ++m)
                            for (unsigned int n=0; n<dim; ++n)
                              result[i][j][l][m][n]
                              += (fourth[k][j][l][m][n]
                                  *
                                  data.mapping_support_points[k][i]);

                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          for (unsigned int n=0; n<dim; ++n)
                            jacobian_3rd_derivatives[point][i][j][l][m][n] = result[i][j][l][m][n];
                }
            }
        }
    }

    /**
     * Update the Hessian gradient of the transformation from unit to real cell, the
     * Jacobian Hessians, pushed forward to the real cell coordinates.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_pushed_forward_3rd_derivatives (const CellSimilarity::Similarity                                   cell_similarity,
                                                          const typename QProjector<dim>::DataSetDescriptor                  data_set,
                                                          const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data,
                                                          std::vector<Tensor<5,spacedim> >                      &jacobian_pushed_forward_3rd_derivatives)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_pushed_forward_3rd_derivatives)
        {
          const unsigned int n_q_points = jacobian_pushed_forward_3rd_derivatives.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              double tmp[spacedim][spacedim][spacedim][spacedim][spacedim];
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const Tensor<4,dim> *fourth =
                    &data.fourth_derivative(point+data_set, 0);
                  double result [spacedim][dim][dim][dim][dim];
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          for (unsigned int n=0; n<dim; ++n)
                            result[i][j][l][m][n] = (fourth[0][j][l][m][n] *
                                                     data.mapping_support_points[0][i]);
                  for (unsigned int k=1; k<data.n_shape_functions; ++k)
                    for (unsigned int i=0; i<spacedim; ++i)
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          for (unsigned int m=0; m<dim; ++m)
                            for (unsigned int n=0; n<dim; ++n)
                              result[i][j][l][m][n]
                              += (fourth[k][j][l][m][n]
                                  *
                                  data.mapping_support_points[k][i]);

                  // push-forward the j-coordinate
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          for (unsigned int n=0; n<dim; ++n)
                            {
                              tmp[i][j][l][m][n] = result[i][0][l][m][n] *
                                                   data.covariant[point][j][0];
                              for (unsigned int jr=1; jr<dim; ++jr)
                                tmp[i][j][l][m][n] += result[i][jr][l][m][n] *
                                                      data.covariant[point][j][jr];
                            }

                  // push-forward the l-coordinate
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<spacedim; ++l)
                        for (unsigned int m=0; m<dim; ++m)
                          for (unsigned int n=0; n<dim; ++n)
                            {
                              jacobian_pushed_forward_3rd_derivatives[point][i][j][l][m][n]
                                = tmp[i][j][0][m][n] *
                                  data.covariant[point][l][0];
                              for (unsigned int lr=1; lr<dim; ++lr)
                                jacobian_pushed_forward_3rd_derivatives[point][i][j][l][m][n]
                                += tmp[i][j][lr][m][n] *
                                   data.covariant[point][l][lr];
                            }

                  // push-forward the m-coordinate
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<spacedim; ++l)
                        for (unsigned int m=0; m<spacedim; ++m)
                          for (unsigned int n=0; n<dim; ++n)
                            {
                              tmp[i][j][l][m][n]
                                = jacobian_pushed_forward_3rd_derivatives[point][i][j][l][0][n] *
                                  data.covariant[point][m][0];
                              for (unsigned int mr=1; mr<dim; ++mr)
                                tmp[i][j][l][m][n]
                                += jacobian_pushed_forward_3rd_derivatives[point][i][j][l][mr][n] *
                                   data.covariant[point][m][mr];
                            }

                  // push-forward the n-coordinate
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<spacedim; ++j)
                      for (unsigned int l=0; l<spacedim; ++l)
                        for (unsigned int m=0; m<spacedim; ++m)
                          for (unsigned int n=0; n<spacedim; ++n)
                            {
                              jacobian_pushed_forward_3rd_derivatives[point][i][j][l][m][n]
                                = tmp[i][j][l][m][0] *
                                  data.covariant[point][n][0];
                              for (unsigned int nr=1; nr<dim; ++nr)
                                jacobian_pushed_forward_3rd_derivatives[point][i][j][l][m][n]
                                += tmp[i][j][l][m][nr] *
                                   data.covariant[point][n][nr];
                            }
                }
            }
        }
    }
  }
}




template<int dim, int spacedim>
CellSimilarity::Similarity
MappingManifold<dim,spacedim>::
fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                const CellSimilarity::Similarity                           cell_similarity,
                const Quadrature<dim>                                     &quadrature,
                const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // ensure that the following static_cast is really correct:
  Assert (dynamic_cast<const InternalData *>(&internal_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  const unsigned int n_q_points=quadrature.size();

  // if necessary, recompute the support points of the transformation of this cell
  // (note that we need to first check the triangulation pointer, since otherwise
  // the second test might trigger an exception if the triangulations are not the
  // same)
  if ((data.mapping_support_points.size() == 0)
      ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation())
      ||
      (cell != data.cell_of_current_support_points))
    {
      data.mapping_support_points = this->compute_mapping_support_points(cell);
      data.cell_of_current_support_points = cell;
    }

  internal::maybe_compute_q_points<dim,spacedim> (QProjector<dim>::DataSetDescriptor::cell (),
                                                  data,
                                                  output_data.quadrature_points);
  internal::maybe_update_Jacobians<dim,spacedim> (cell_similarity,
                                                  QProjector<dim>::DataSetDescriptor::cell (),
                                                  data);

  const UpdateFlags update_flags = data.update_each;
  const std::vector<double> &weights=quadrature.get_weights();

  // Multiply quadrature weights by absolute value of Jacobian determinants or
  // the area element g=sqrt(DX^t DX) in case of codim > 0

  if (update_flags & (update_normal_vectors
                      | update_JxW_values))
    {
      AssertDimension (output_data.JxW_values.size(), n_q_points);

      Assert( !(update_flags & update_normal_vectors ) ||
              (output_data.normal_vectors.size() == n_q_points),
              ExcDimensionMismatch(output_data.normal_vectors.size(), n_q_points));


      if (cell_similarity != CellSimilarity::translation)
        for (unsigned int point=0; point<n_q_points; ++point)
          {

            if (dim == spacedim)
              {
                const double det = data.contravariant[point].determinant();

                // check for distorted cells.

                // TODO: this allows for anisotropies of up to 1e6 in 3D and
                // 1e12 in 2D. might want to find a finer
                // (dimension-independent) criterion
                Assert (det > 1e-12*Utilities::fixed_power<dim>(cell->diameter()/
                                                                std::sqrt(double(dim))),
                        (typename Mapping<dim,spacedim>::ExcDistortedMappedCell(cell->center(), det, point)));

                output_data.JxW_values[point] = weights[point] * det;
              }
            // if dim==spacedim, then there is no cell normal to
            // compute. since this is for FEValues (and not FEFaceValues),
            // there are also no face normals to compute
            else //codim>0 case
              {
                Tensor<1, spacedim> DX_t [dim];
                for (unsigned int i=0; i<spacedim; ++i)
                  for (unsigned int j=0; j<dim; ++j)
                    DX_t[j][i] = data.contravariant[point][i][j];

                Tensor<2, dim> G; //First fundamental form
                for (unsigned int i=0; i<dim; ++i)
                  for (unsigned int j=0; j<dim; ++j)
                    G[i][j] = DX_t[i] * DX_t[j];

                output_data.JxW_values[point]
                  = sqrt(determinant(G)) * weights[point];

                if (cell_similarity == CellSimilarity::inverted_translation)
                  {
                    // we only need to flip the normal
                    if (update_flags & update_normal_vectors)
                      output_data.normal_vectors[point] *= -1.;
                  }
                else
                  {
                    const unsigned int codim = spacedim-dim;
                    (void)codim;

                    if (update_flags & update_normal_vectors)
                      {
                        Assert( codim==1 , ExcMessage("There is no cell normal in codim 2."));

                        if (dim==1)
                          output_data.normal_vectors[point] =
                            cross_product_2d(-DX_t[0]);
                        else //dim == 2
                          output_data.normal_vectors[point] =
                            cross_product_3d(DX_t[0], DX_t[1]);

                        output_data.normal_vectors[point] /= output_data.normal_vectors[point].norm();

                        if (cell->direction_flag() == false)
                          output_data.normal_vectors[point] *= -1.;
                      }

                  }
              } //codim>0 case

          }
    }



  // copy values from InternalData to vector given by reference
  if (update_flags & update_jacobians)
    {
      AssertDimension (output_data.jacobians.size(), n_q_points);
      if (cell_similarity != CellSimilarity::translation)
        for (unsigned int point=0; point<n_q_points; ++point)
          output_data.jacobians[point] = data.contravariant[point];
    }

  // copy values from InternalData to vector given by reference
  if (update_flags & update_inverse_jacobians)
    {
      AssertDimension (output_data.inverse_jacobians.size(), n_q_points);
      if (cell_similarity != CellSimilarity::translation)
        for (unsigned int point=0; point<n_q_points; ++point)
          output_data.inverse_jacobians[point] = data.covariant[point].transpose();
    }

  internal::maybe_update_jacobian_grads<dim,spacedim> (cell_similarity,
                                                       QProjector<dim>::DataSetDescriptor::cell (),
                                                       data,
                                                       output_data.jacobian_grads);

  internal::maybe_update_jacobian_pushed_forward_grads<dim,spacedim> (cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell (),
      data,
      output_data.jacobian_pushed_forward_grads);

  internal::maybe_update_jacobian_2nd_derivatives<dim,spacedim> (cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell (),
      data,
      output_data.jacobian_2nd_derivatives);

  internal::maybe_update_jacobian_pushed_forward_2nd_derivatives<dim,spacedim> (cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell (),
      data,
      output_data.jacobian_pushed_forward_2nd_derivatives);

  internal::maybe_update_jacobian_3rd_derivatives<dim,spacedim> (cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell (),
      data,
      output_data.jacobian_3rd_derivatives);

  internal::maybe_update_jacobian_pushed_forward_3rd_derivatives<dim,spacedim> (cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell (),
      data,
      output_data.jacobian_pushed_forward_3rd_derivatives);

  return cell_similarity;
}






namespace internal
{
  namespace
  {
    /**
     * Depending on what information is called for in the update flags of the
     * @p data object, compute the various pieces of information that is required
     * by the fill_fe_face_values() and fill_fe_subface_values() functions.
     * This function simply unifies the work that would be done by
     * those two functions.
     *
     * The resulting data is put into the @p output_data argument.
     */
    template <int dim, int spacedim>
    void
    maybe_compute_face_data (const dealii::MappingManifold<dim,spacedim> &mapping,
                             const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                             const unsigned int               face_no,
                             const unsigned int               subface_no,
                             const unsigned int               n_q_points,
                             const std::vector<double>        &weights,
                             const typename dealii::MappingManifold<dim,spacedim>::InternalData &data,
                             internal::FEValues::MappingRelatedData<dim,spacedim>         &output_data)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_boundary_forms)
        {
          AssertDimension (output_data.boundary_forms.size(), n_q_points);
          if (update_flags & update_normal_vectors)
            AssertDimension (output_data.normal_vectors.size(), n_q_points);
          if (update_flags & update_JxW_values)
            AssertDimension (output_data.JxW_values.size(), n_q_points);

          // map the unit tangentials to the real cell. checking for d!=dim-1
          // eliminates compiler warnings regarding unsigned int expressions <
          // 0.
          for (unsigned int d=0; d!=dim-1; ++d)
            {
              Assert (face_no+GeometryInfo<dim>::faces_per_cell*d <
                      data.unit_tangentials.size(),
                      ExcInternalError());
              Assert (data.aux[d].size() <=
                      data.unit_tangentials[face_no+GeometryInfo<dim>::faces_per_cell*d].size(),
                      ExcInternalError());

              mapping.transform (make_array_view(data.unit_tangentials[face_no+GeometryInfo<dim>::faces_per_cell*d]),
                                 mapping_contravariant,
                                 data,
                                 make_array_view(data.aux[d]));
            }

          // if dim==spacedim, we can use the unit tangentials to compute the
          // boundary form by simply taking the cross product
          if (dim == spacedim)
            {
              for (unsigned int i=0; i<n_q_points; ++i)
                switch (dim)
                  {
                  case 1:
                    // in 1d, we don't have access to any of the data.aux
                    // fields (because it has only dim-1 components), but we
                    // can still compute the boundary form by simply
                    // looking at the number of the face
                    output_data.boundary_forms[i][0] = (face_no == 0 ?
                                                        -1 : +1);
                    break;
                  case 2:
                    output_data.boundary_forms[i] =
                      cross_product_2d(data.aux[0][i]);
                    break;
                  case 3:
                    output_data.boundary_forms[i] =
                      cross_product_3d(data.aux[0][i], data.aux[1][i]);
                    break;
                  default:
                    Assert(false, ExcNotImplemented());
                  }
            }
          else //(dim < spacedim)
            {
              // in the codim-one case, the boundary form results from the
              // cross product of all the face tangential vectors and the cell
              // normal vector
              //
              // to compute the cell normal, use the same method used in
              // fill_fe_values for cells above
              AssertDimension (data.contravariant.size(), n_q_points);

              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  if (dim==1)
                    {
                      // J is a tangent vector
                      output_data.boundary_forms[point] = data.contravariant[point].transpose()[0];
                      output_data.boundary_forms[point] /=
                        (face_no == 0 ? -1. : +1.) * output_data.boundary_forms[point].norm();
                    }

                  if (dim==2)
                    {
                      const DerivativeForm<1,spacedim,dim> DX_t =
                        data.contravariant[point].transpose();

                      Tensor<1, spacedim> cell_normal =
                        cross_product_3d(DX_t[0], DX_t[1]);
                      cell_normal /= cell_normal.norm();

                      // then compute the face normal from the face tangent
                      // and the cell normal:
                      output_data.boundary_forms[point] =
                        cross_product_3d(data.aux[0][point], cell_normal);
                    }
                }
            }

          if (update_flags & (update_normal_vectors
                              | update_JxW_values))
            for (unsigned int i=0; i<output_data.boundary_forms.size(); ++i)
              {
                if (update_flags & update_JxW_values)
                  {
                    output_data.JxW_values[i] = output_data.boundary_forms[i].norm() * weights[i];

                    if (subface_no!=numbers::invalid_unsigned_int)
                      {
                        const double area_ratio=GeometryInfo<dim>::subface_ratio(
                                                  cell->subface_case(face_no), subface_no);
                        output_data.JxW_values[i] *= area_ratio;
                      }
                  }

                if (update_flags & update_normal_vectors)
                  output_data.normal_vectors[i] = Point<spacedim>(output_data.boundary_forms[i] /
                                                                  output_data.boundary_forms[i].norm());
              }

          if (update_flags & update_jacobians)
            for (unsigned int point=0; point<n_q_points; ++point)
              output_data.jacobians[point] = data.contravariant[point];

          if (update_flags & update_inverse_jacobians)
            for (unsigned int point=0; point<n_q_points; ++point)
              output_data.inverse_jacobians[point] = data.covariant[point].transpose();
        }
    }


    /**
     * Do the work of MappingManifold::fill_fe_face_values() and
     * MappingManifold::fill_fe_subface_values() in a generic way,
     * using the 'data_set' to differentiate whether we will
     * work on a face (and if so, which one) or subface.
     */
    template<int dim, int spacedim>
    void
    do_fill_fe_face_values (const dealii::MappingManifold<dim,spacedim>                             &mapping,
                            const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                            const unsigned int                                                 face_no,
                            const unsigned int                                                 subface_no,
                            const typename QProjector<dim>::DataSetDescriptor                  data_set,
                            const Quadrature<dim-1>                                           &quadrature,
                            const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data,
                            internal::FEValues::MappingRelatedData<dim,spacedim>              &output_data)
    {
      maybe_compute_q_points<dim,spacedim> (data_set,
                                            data,
                                            output_data.quadrature_points);
      maybe_update_Jacobians<dim,spacedim> (CellSimilarity::none,
                                            data_set,
                                            data);
      maybe_update_jacobian_grads<dim,spacedim> (CellSimilarity::none,
                                                 data_set,
                                                 data,
                                                 output_data.jacobian_grads);
      maybe_update_jacobian_pushed_forward_grads<dim,spacedim> (CellSimilarity::none,
                                                                data_set,
                                                                data,
                                                                output_data.jacobian_pushed_forward_grads);
      maybe_update_jacobian_2nd_derivatives<dim,spacedim> (CellSimilarity::none,
                                                           data_set,
                                                           data,
                                                           output_data.jacobian_2nd_derivatives);
      maybe_update_jacobian_pushed_forward_2nd_derivatives<dim,spacedim> (CellSimilarity::none,
          data_set,
          data,
          output_data.jacobian_pushed_forward_2nd_derivatives);
      maybe_update_jacobian_3rd_derivatives<dim,spacedim> (CellSimilarity::none,
                                                           data_set,
                                                           data,
                                                           output_data.jacobian_3rd_derivatives);
      maybe_update_jacobian_pushed_forward_3rd_derivatives<dim,spacedim> (CellSimilarity::none,
          data_set,
          data,
          output_data.jacobian_pushed_forward_3rd_derivatives);

      maybe_compute_face_data (mapping,
                               cell, face_no, subface_no, quadrature.size(),
                               quadrature.get_weights(), data,
                               output_data);
    }
  }
}



template<int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int                                         face_no,
                     const Quadrature<dim-1>                                   &quadrature,
                     const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                     internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // ensure that the following cast is really correct:
  Assert ((dynamic_cast<const InternalData *>(&internal_data) != 0),
          ExcInternalError());
  const InternalData &data
    = static_cast<const InternalData &>(internal_data);

  // if necessary, recompute the support points of the transformation of this cell
  // (note that we need to first check the triangulation pointer, since otherwise
  // the second test might trigger an exception if the triangulations are not the
  // same)
  if ((data.mapping_support_points.size() == 0)
      ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation())
      ||
      (cell != data.cell_of_current_support_points))
    {
      data.mapping_support_points = this->compute_mapping_support_points(cell);
      data.cell_of_current_support_points = cell;
    }

  internal::do_fill_fe_face_values (*this,
                                    cell, face_no, numbers::invalid_unsigned_int,
                                    QProjector<dim>::DataSetDescriptor::face (face_no,
                                        cell->face_orientation(face_no),
                                        cell->face_flip(face_no),
                                        cell->face_rotation(face_no),
                                        quadrature.size()),
                                    quadrature,
                                    data,
                                    output_data);
}



template<int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                        const unsigned int                                         face_no,
                        const unsigned int                                         subface_no,
                        const Quadrature<dim-1>                                   &quadrature,
                        const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                        internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // ensure that the following cast is really correct:
  Assert ((dynamic_cast<const InternalData *>(&internal_data) != 0),
          ExcInternalError());
  const InternalData &data
    = static_cast<const InternalData &>(internal_data);

  // if necessary, recompute the support points of the transformation of this cell
  // (note that we need to first check the triangulation pointer, since otherwise
  // the second test might trigger an exception if the triangulations are not the
  // same)
  if ((data.mapping_support_points.size() == 0)
      ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation())
      ||
      (cell != data.cell_of_current_support_points))
    {
      data.mapping_support_points = this->compute_mapping_support_points(cell);
      data.cell_of_current_support_points = cell;
    }

  internal::do_fill_fe_face_values (*this,
                                    cell, face_no, subface_no,
                                    QProjector<dim>::DataSetDescriptor::subface (face_no, subface_no,
                                        cell->face_orientation(face_no),
                                        cell->face_flip(face_no),
                                        cell->face_rotation(face_no),
                                        quadrature.size(),
                                        cell->subface_case(face_no)),
                                    quadrature,
                                    data,
                                    output_data);
}



namespace
{
  template <int dim, int spacedim, int rank>
  void
  transform_fields(const ArrayView<const Tensor<rank,dim> >               &input,
                   const MappingType                                       mapping_type,
                   const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                   const ArrayView<Tensor<rank,spacedim> >                &output)
  {
    AssertDimension (input.size(), output.size());
    Assert ((dynamic_cast<const typename MappingManifold<dim,spacedim>::InternalData *>(&mapping_data) != 0),
            ExcInternalError());
    const typename MappingManifold<dim,spacedim>::InternalData
    &data = static_cast<const typename MappingManifold<dim,spacedim>::InternalData &>(mapping_data);

    switch (mapping_type)
      {
      case mapping_contravariant:
      {
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

        for (unsigned int i=0; i<output.size(); ++i)
          output[i] = apply_transformation(data.contravariant[i], input[i]);

        return;
      }

      case mapping_piola:
      {
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
        Assert (data.update_each & update_volume_elements,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));
        Assert (rank==1, ExcMessage("Only for rank 1"));
        if (rank!=1)
          return;

        for (unsigned int i=0; i<output.size(); ++i)
          {
            output[i] = apply_transformation(data.contravariant[i], input[i]);
            output[i] /= data.volume_elements[i];
          }
        return;
      }
      //We still allow this operation as in the
      //reference cell Derivatives are Tensor
      //rather than DerivativeForm
      case mapping_covariant:
      {
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

        for (unsigned int i=0; i<output.size(); ++i)
          output[i] = apply_transformation(data.covariant[i], input[i]);

        return;
      }

      default:
        Assert(false, ExcNotImplemented());
      }
  }


  template <int dim, int spacedim, int rank>
  void
  transform_gradients(const ArrayView<const Tensor<rank,dim> >                &input,
                      const MappingType                                        mapping_type,
                      const typename Mapping<dim,spacedim>::InternalDataBase  &mapping_data,
                      const ArrayView<Tensor<rank,spacedim> >                 &output)
  {
    AssertDimension (input.size(), output.size());
    Assert ((dynamic_cast<const typename MappingManifold<dim,spacedim>::InternalData *>(&mapping_data) != 0),
            ExcInternalError());
    const typename MappingManifold<dim,spacedim>::InternalData
    &data = static_cast<const typename MappingManifold<dim,spacedim>::InternalData &>(mapping_data);

    switch (mapping_type)
      {
      case mapping_contravariant_gradient:
      {
        Assert (data.update_each & update_covariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
        Assert (rank==2, ExcMessage("Only for rank 2"));

        for (unsigned int i=0; i<output.size(); ++i)
          {
            DerivativeForm<1,spacedim,dim> A =
              apply_transformation(data.contravariant[i], transpose(input[i]) );
            output[i] = apply_transformation(data.covariant[i], A.transpose() );
          }

        return;
      }

      case mapping_covariant_gradient:
      {
        Assert (data.update_each & update_covariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
        Assert (rank==2, ExcMessage("Only for rank 2"));

        for (unsigned int i=0; i<output.size(); ++i)
          {
            DerivativeForm<1,spacedim,dim> A =
              apply_transformation(data.covariant[i], transpose(input[i]) );
            output[i] = apply_transformation(data.covariant[i], A.transpose() );
          }

        return;
      }

      case mapping_piola_gradient:
      {
        Assert (data.update_each & update_covariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
        Assert (data.update_each & update_volume_elements,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));
        Assert (rank==2, ExcMessage("Only for rank 2"));

        for (unsigned int i=0; i<output.size(); ++i)
          {
            DerivativeForm<1,spacedim,dim> A =
              apply_transformation(data.covariant[i], input[i] );
            Tensor<2,spacedim> T =
              apply_transformation(data.contravariant[i], A.transpose() );

            output[i] = transpose(T);
            output[i] /= data.volume_elements[i];
          }

        return;
      }

      default:
        Assert(false, ExcNotImplemented());
      }
  }




  template <int dim, int spacedim>
  void
  transform_hessians(const ArrayView<const Tensor<3,dim> >                  &input,
                     const MappingType                                       mapping_type,
                     const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                     const ArrayView<Tensor<3,spacedim> >                   &output)
  {
    AssertDimension (input.size(), output.size());
    Assert ((dynamic_cast<const typename MappingManifold<dim,spacedim>::InternalData *>(&mapping_data) != 0),
            ExcInternalError());
    const typename MappingManifold<dim,spacedim>::InternalData
    &data = static_cast<const typename MappingManifold<dim,spacedim>::InternalData &>(mapping_data);

    switch (mapping_type)
      {
      case mapping_contravariant_hessian:
      {
        Assert (data.update_each & update_covariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

        for (unsigned int q=0; q<output.size(); ++q)
          for (unsigned int i=0; i<spacedim; ++i)
            {
              double tmp1[dim][dim];
              for (unsigned int J=0; J<dim; ++J)
                for (unsigned int K=0; K<dim; ++K)
                  {
                    tmp1[J][K] = data.contravariant[q][i][0] * input[q][0][J][K];
                    for (unsigned int I=1; I<dim; ++I)
                      tmp1[J][K] += data.contravariant[q][i][I] * input[q][I][J][K];
                  }
              for (unsigned int j=0; j<spacedim; ++j)
                {
                  double tmp2[dim];
                  for (unsigned int K=0; K<dim; ++K)
                    {
                      tmp2[K] = data.covariant[q][j][0] * tmp1[0][K];
                      for (unsigned int J=1; J<dim; ++J)
                        tmp2[K] += data.covariant[q][j][J] * tmp1[J][K];
                    }
                  for (unsigned int k=0; k<spacedim; ++k)
                    {
                      output[q][i][j][k] = data.covariant[q][k][0] * tmp2[0];
                      for (unsigned int K=1; K<dim; ++K)
                        output[q][i][j][k] += data.covariant[q][k][K] * tmp2[K];
                    }
                }
            }
        return;
      }

      case mapping_covariant_hessian:
      {
        Assert (data.update_each & update_covariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

        for (unsigned int q=0; q<output.size(); ++q)
          for (unsigned int i=0; i<spacedim; ++i)
            {
              double tmp1[dim][dim];
              for (unsigned int J=0; J<dim; ++J)
                for (unsigned int K=0; K<dim; ++K)
                  {
                    tmp1[J][K] = data.covariant[q][i][0] * input[q][0][J][K];
                    for (unsigned int I=1; I<dim; ++I)
                      tmp1[J][K] += data.covariant[q][i][I] * input[q][I][J][K];
                  }
              for (unsigned int j=0; j<spacedim; ++j)
                {
                  double tmp2[dim];
                  for (unsigned int K=0; K<dim; ++K)
                    {
                      tmp2[K] = data.covariant[q][j][0] * tmp1[0][K];
                      for (unsigned int J=1; J<dim; ++J)
                        tmp2[K] += data.covariant[q][j][J] * tmp1[J][K];
                    }
                  for (unsigned int k=0; k<spacedim; ++k)
                    {
                      output[q][i][j][k] = data.covariant[q][k][0] * tmp2[0];
                      for (unsigned int K=1; K<dim; ++K)
                        output[q][i][j][k] += data.covariant[q][k][K] * tmp2[K];
                    }
                }
            }

        return;
      }

      case mapping_piola_hessian:
      {
        Assert (data.update_each & update_covariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
        Assert (data.update_each & update_volume_elements,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));

        for (unsigned int q=0; q<output.size(); ++q)
          for (unsigned int i=0; i<spacedim; ++i)
            {
              double factor[dim];
              for (unsigned int I=0; I<dim; ++I)
                factor[I] = data.contravariant[q][i][I] / data.volume_elements[q];
              double tmp1[dim][dim];
              for (unsigned int J=0; J<dim; ++J)
                for (unsigned int K=0; K<dim; ++K)
                  {
                    tmp1[J][K] = factor[0] * input[q][0][J][K];
                    for (unsigned int I=1; I<dim; ++I)
                      tmp1[J][K] += factor[I] * input[q][I][J][K];
                  }
              for (unsigned int j=0; j<spacedim; ++j)
                {
                  double tmp2[dim];
                  for (unsigned int K=0; K<dim; ++K)
                    {
                      tmp2[K] = data.covariant[q][j][0] * tmp1[0][K];
                      for (unsigned int J=1; J<dim; ++J)
                        tmp2[K] += data.covariant[q][j][J] * tmp1[J][K];
                    }
                  for (unsigned int k=0; k<spacedim; ++k)
                    {
                      output[q][i][j][k] = data.covariant[q][k][0] * tmp2[0];
                      for (unsigned int K=1; K<dim; ++K)
                        output[q][i][j][k] += data.covariant[q][k][K] * tmp2[K];
                    }
                }
            }

        return;
      }

      default:
        Assert(false, ExcNotImplemented());
      }
  }




  template<int dim, int spacedim, int rank>
  void
  transform_differential_forms(const ArrayView<const DerivativeForm<rank, dim,spacedim> >   &input,
                               const MappingType                                             mapping_type,
                               const typename Mapping<dim,spacedim>::InternalDataBase       &mapping_data,
                               const ArrayView<Tensor<rank+1, spacedim> >                   &output)
  {
    AssertDimension (input.size(), output.size());
    Assert ((dynamic_cast<const typename MappingManifold<dim,spacedim>::InternalData *>(&mapping_data) != 0),
            ExcInternalError());
    const typename MappingManifold<dim,spacedim>::InternalData
    &data = static_cast<const typename MappingManifold<dim,spacedim>::InternalData &>(mapping_data);

    switch (mapping_type)
      {
      case mapping_covariant:
      {
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

        for (unsigned int i=0; i<output.size(); ++i)
          output[i] = apply_transformation(data.covariant[i], input[i]);

        return;
      }
      default:
        Assert(false, ExcNotImplemented());
      }
  }
}



template<int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
transform (const ArrayView<const Tensor<1, dim> >                  &input,
           const MappingType                                        mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase  &mapping_data,
           const ArrayView<Tensor<1, spacedim> >                   &output) const
{
  transform_fields(input, mapping_type, mapping_data, output);
}



template<int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
transform (const ArrayView<const DerivativeForm<1, dim,spacedim> >  &input,
           const MappingType                                         mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase   &mapping_data,
           const ArrayView<Tensor<2, spacedim> >                    &output) const
{
  transform_differential_forms(input, mapping_type, mapping_data, output);
}



template<int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
transform (const ArrayView<const Tensor<2, dim> >                  &input,
           const MappingType                                        mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase  &mapping_data,
           const ArrayView<Tensor<2, spacedim> >                   &output) const
{
  switch (mapping_type)
    {
    case mapping_contravariant:
      transform_fields(input, mapping_type, mapping_data, output);
      return;

    case mapping_piola_gradient:
    case mapping_contravariant_gradient:
    case mapping_covariant_gradient:
      transform_gradients(input, mapping_type, mapping_data, output);
      return;
    default:
      Assert(false, ExcNotImplemented());
    }
}



template<int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
transform (const ArrayView<const  DerivativeForm<2, dim, spacedim> > &input,
           const MappingType                                          mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase    &mapping_data,
           const ArrayView<Tensor<3,spacedim> >                      &output) const
{

  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_covariant_gradient:
    {
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

      for (unsigned int q=0; q<output.size(); ++q)
        for (unsigned int i=0; i<spacedim; ++i)
          for (unsigned int j=0; j<spacedim; ++j)
            {
              double tmp[dim];
              for (unsigned int K=0; K<dim; ++K)
                {
                  tmp[K] = data.covariant[q][j][0] * input[q][i][0][K];
                  for (unsigned int J=1; J<dim; ++J)
                    tmp[K] += data.covariant[q][j][J] * input[q][i][J][K];
                }
              for (unsigned int k=0; k<spacedim; ++k)
                {
                  output[q][i][j][k] = data.covariant[q][k][0] * tmp[0];
                  for (unsigned int K=1; K<dim; ++K)
                    output[q][i][j][k] += data.covariant[q][k][K] * tmp[K];
                }
            }
      return;
    }

    default:
      Assert(false, ExcNotImplemented());
    }
}



template<int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
transform (const ArrayView<const  Tensor<3,dim> >                  &input,
           const MappingType                                        mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase  &mapping_data,
           const ArrayView<Tensor<3,spacedim> >                    &output) const
{
  switch (mapping_type)
    {
    case mapping_piola_hessian:
    case mapping_contravariant_hessian:
    case mapping_covariant_hessian:
      transform_hessians(input, mapping_type, mapping_data, output);
      return;
    default:
      Assert(false, ExcNotImplemented());
    }
}



namespace
{
  /**
   * Ask the manifold descriptor to return intermediate points on lines or
   * faces. The function needs to return one or multiple points (depending on
   * the number of elements in the output vector @p points that lie inside a
   * line, quad or hex). Whether it is a line, quad or hex doesn't really
   * matter to this function but it can be inferred from the number of input
   * points in the @p surrounding_points vector.
   */
  template<int dim, int spacedim>
  void
  get_intermediate_points (const Manifold<dim, spacedim> &manifold,
                           const QGaussLobatto<1>        &line_support_points,
                           const std::vector<Point<spacedim> > &surrounding_points,
                           std::vector<Point<spacedim> > &points)
  {
    Assert(surrounding_points.size() >= 2, ExcMessage("At least 2 surrounding points are required"));
    const unsigned int n=points.size();
    Assert(n>0, ExcMessage("You can't ask for 0 intermediate points."));
    std::vector<double> w(surrounding_points.size());

    switch (surrounding_points.size())
      {
      case 2:
      {
        // If two points are passed, these are the two vertices, and
        // we can only compute degree-1 intermediate points.
        for (unsigned int i=0; i<n; ++i)
          {
            const double x = line_support_points.point(i+1)[0];
            w[1] = x;
            w[0] = (1-x);
            Quadrature<spacedim> quadrature(surrounding_points, w);
            points[i] = manifold.get_new_point(quadrature);
          }
        break;
      }

      case 4:
      {
        Assert(spacedim >= 2, ExcImpossibleInDim(spacedim));
        const unsigned m=
          static_cast<unsigned int>(std::sqrt(static_cast<double>(n)));
        // is n a square number
        Assert(m*m==n, ExcInternalError());

        // If four points are passed, these are the two vertices, and
        // we can only compute (degree-1)*(degree-1) intermediate
        // points.
        for (unsigned int i=0; i<m; ++i)
          {
            const double y=line_support_points.point(1+i)[0];
            for (unsigned int j=0; j<m; ++j)
              {
                const double x=line_support_points.point(1+j)[0];

                w[0] = (1-x)*(1-y);
                w[1] =     x*(1-y);
                w[2] = (1-x)*y    ;
                w[3] =     x*y    ;
                Quadrature<spacedim> quadrature(surrounding_points, w);
                points[i*m+j]=manifold.get_new_point(quadrature);
              }
          }
        break;
      }

      case 8:
        Assert(false, ExcNotImplemented());
        break;
      default:
        Assert(false, ExcInternalError());
        break;
      }
  }




  /**
   * Ask the manifold descriptor to return intermediate points on the object
   * pointed to by the TriaIterator @p iter. This function tries to be
   * backward compatible with respect to the differences between
   * Boundary<dim,spacedim> and Manifold<dim,spacedim>, querying the first
   * whenever the passed @p manifold can be upgraded to a
   * Boundary<dim,spacedim>.
   */
  template <int dim, int spacedim, class TriaIterator>
  void get_intermediate_points_on_object(const Manifold<dim, spacedim> &manifold,
                                         const QGaussLobatto<1>        &line_support_points,
                                         const TriaIterator &iter,
                                         std::vector<Point<spacedim> > &points)
  {
    const unsigned int structdim = TriaIterator::AccessorType::structure_dimension;

    // Try backward compatibility option.
    if (const Boundary<dim,spacedim> *boundary
        = dynamic_cast<const Boundary<dim,spacedim> *>(&manifold))
      // This is actually a boundary. Call old methods.
      {
        switch (structdim)
          {
          case 1:
          {
            const typename Triangulation<dim,spacedim>::line_iterator line = iter;
            boundary->get_intermediate_points_on_line(line, points);
            return;
          }
          case 2:
          {
            const typename Triangulation<dim,spacedim>::quad_iterator quad = iter;
            boundary->get_intermediate_points_on_quad(quad, points);
            return;
          }
          default:
            Assert(false, ExcInternalError());
            return;
          }
      }
    else
      {
        std::vector<Point<spacedim> > sp(GeometryInfo<structdim>::vertices_per_cell);
        for (unsigned int i=0; i<sp.size(); ++i)
          sp[i] = iter->vertex(i);
        get_intermediate_points(manifold, line_support_points, sp, points);
      }
  }


  /**
   * Take a <tt>support_point_weights_on_hex(quad)</tt> and apply it to the vector
   * @p a to compute the inner support points as a linear combination of the
   * exterior points.
   *
   * The vector @p a initially contains the locations of the @p n_outer
   * points, the @p n_inner computed inner points are appended.
   *
   * See equation (7) of the `mapping' report.
   */
  template <int spacedim>
  void add_weighted_interior_points(const Table<2,double>   &lvs,
                                    std::vector<Point<spacedim> > &a)
  {
    const unsigned int n_inner_apply=lvs.n_rows();
    const unsigned int n_outer_apply=lvs.n_cols();
    Assert(a.size()==n_outer_apply,
           ExcDimensionMismatch(a.size(), n_outer_apply));

    // compute each inner point as linear combination of the outer points. the
    // weights are given by the lvs entries, the outer points are the first
    // (existing) elements of a
    for (unsigned int unit_point=0; unit_point<n_inner_apply; ++unit_point)
      {
        Assert(lvs.n_cols()==n_outer_apply, ExcInternalError());
        Point<spacedim> p;
        for (unsigned int k=0; k<n_outer_apply; ++k)
          p+=lvs[unit_point][k]*a[k];

        a.push_back(p);
      }
  }
}


template <int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
add_line_support_points (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                         std::vector<Point<spacedim> > &a) const
{
  // if we only need the midpoint, then ask for it.
  if (this->polynomial_degree==2)
    {
      for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
        {
          const typename Triangulation<dim,spacedim>::line_iterator line =
            (dim == 1  ?
             static_cast<typename Triangulation<dim,spacedim>::line_iterator>(cell) :
             cell->line(line_no));

          const Manifold<dim,spacedim> &manifold =
            ( ( line->manifold_id() == numbers::invalid_manifold_id ) &&
              ( dim < spacedim )
              ?
              cell->get_manifold()
              :
              line->get_manifold() );
          a.push_back(manifold.get_new_point_on_line(line));
        }
    }
  else
    // otherwise call the more complicated functions and ask for inner points
    // from the boundary description
    {
      std::vector<Point<spacedim> > line_points (this->polynomial_degree-1);
      // loop over each of the lines, and if it is at the boundary, then first
      // get the boundary description and second compute the points on it
      for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
        {
          const typename Triangulation<dim,spacedim>::line_iterator
          line = (dim == 1
                  ?
                  static_cast<typename Triangulation<dim,spacedim>::line_iterator>(cell)
                  :
                  cell->line(line_no));

          const Manifold<dim,spacedim> &manifold =
            ( ( line->manifold_id() == numbers::invalid_manifold_id ) &&
              ( dim < spacedim )
              ?
              cell->get_manifold() :
              line->get_manifold() );

          //          get_intermediate_points_on_object (manifold, line_support_points, line, line_points);

          if (dim==3)
            {
              // in 3D, lines might be in wrong orientation. if so, reverse
              // the vector
              if (cell->line_orientation(line_no))
                a.insert (a.end(), line_points.begin(), line_points.end());
              else
                a.insert (a.end(), line_points.rbegin(), line_points.rend());
            }
          else
            // in 2D, lines always have the correct orientation. simply append
            // all points
            a.insert (a.end(), line_points.begin(), line_points.end());
        }
    }
}



template <>
void
MappingManifold<3,3>::
add_quad_support_points(const Triangulation<3,3>::cell_iterator &cell,
                        std::vector<Point<3> >                &a) const
{
//   const unsigned int faces_per_cell    = GeometryInfo<3>::faces_per_cell,
//                      vertices_per_face = GeometryInfo<3>::vertices_per_face,
//                      lines_per_face    = GeometryInfo<3>::lines_per_face,
//                      vertices_per_cell = GeometryInfo<3>::vertices_per_cell;

//   static const StraightBoundary<3> straight_boundary;
//   // used if face quad at boundary or entirely in the interior of the domain
//   std::vector<Point<3> > quad_points ((polynomial_degree-1)*(polynomial_degree-1));
//   // used if only one line of face quad is at boundary
//   std::vector<Point<3> > b(4*polynomial_degree);

//   // Used by the new Manifold interface. This vector collects the
//   // vertices used to compute the intermediate points.
//   std::vector<Point<3> > vertices(4);

//   // loop over all faces and collect points on them
//   for (unsigned int face_no=0; face_no<faces_per_cell; ++face_no)
//     {
//       const Triangulation<3>::face_iterator face = cell->face(face_no);

//       // select the correct mappings for the present face
//       const bool face_orientation = cell->face_orientation(face_no),
//                  face_flip        = cell->face_flip       (face_no),
//                  face_rotation    = cell->face_rotation   (face_no);

// #ifdef DEBUG
//       // some sanity checks up front
//       for (unsigned int i=0; i<vertices_per_face; ++i)
//         Assert(face->vertex_index(i)==cell->vertex_index(
//                  GeometryInfo<3>::face_to_cell_vertices(face_no, i,
//                                                         face_orientation,
//                                                         face_flip,
//                                                         face_rotation)),
//                ExcInternalError());

//       // indices of the lines that bound a face are given by GeometryInfo<3>::
//       // face_to_cell_lines
//       for (unsigned int i=0; i<lines_per_face; ++i)
//         Assert(face->line(i)==cell->line(GeometryInfo<3>::face_to_cell_lines(
//                                            face_no, i, face_orientation, face_flip, face_rotation)),
//                ExcInternalError());
// #endif

//       // if face at boundary, then ask boundary object to return intermediate
//       // points on it
//       if (face->at_boundary())
//         {
//           get_intermediate_points_on_object(face->get_manifold(), line_support_points, face, quad_points);

//           // in 3D, the orientation, flip and rotation of the face might not
//           // match what we expect here, namely the standard orientation. thus
//           // reorder points accordingly. since a Mapping uses the same shape
//           // function as an FE_Q, we can ask a FE_Q to do the reordering for us.
//           for (unsigned int i=0; i<quad_points.size(); ++i)
//             a.push_back(quad_points[fe_q->adjust_quad_dof_index_for_face_orientation(i,
//                                     face_orientation,
//                                     face_flip,
//                                     face_rotation)]);
//         }
//       else
//         {
//           // face is not at boundary, but maybe some of its lines are. count
//           // them
//           unsigned int lines_at_boundary=0;
//           for (unsigned int i=0; i<lines_per_face; ++i)
//             if (face->line(i)->at_boundary())
//               ++lines_at_boundary;

//           Assert(lines_at_boundary<=lines_per_face, ExcInternalError());

//           // if at least one of the lines bounding this quad is at the
//           // boundary, then collect points separately
//           if (lines_at_boundary>0)
//             {
//               // call of function add_weighted_interior_points increases size of b
//               // about 1. There resize b for the case the mentioned function
//               // was already called.
//               b.resize(4*polynomial_degree);

//               // b is of size 4*degree, make sure that this is the right size
//               Assert(b.size()==vertices_per_face+lines_per_face*(polynomial_degree-1),
//                      ExcDimensionMismatch(b.size(),
//                                           vertices_per_face+lines_per_face*(polynomial_degree-1)));

//               // sort the points into b. We used access from the cell (not
//               // from the face) to fill b, so we can assume a standard face
//               // orientation. Doing so, the calculated points will be in
//               // standard orientation as well.
//               for (unsigned int i=0; i<vertices_per_face; ++i)
//                 b[i]=a[GeometryInfo<3>::face_to_cell_vertices(face_no, i)];

//               for (unsigned int i=0; i<lines_per_face; ++i)
//                 for (unsigned int j=0; j<polynomial_degree-1; ++j)
//                   b[vertices_per_face+i*(polynomial_degree-1)+j]=
//                     a[vertices_per_cell + GeometryInfo<3>::face_to_cell_lines(
//                         face_no, i)*(polynomial_degree-1)+j];

//               // Now b includes the support points on the quad and we can
//               // apply the laplace vector
//               add_weighted_interior_points (support_point_weights_on_quad, b);
//               AssertDimension (b.size(),
//                                4*this->polynomial_degree +
//                                (this->polynomial_degree-1)*(this->polynomial_degree-1));

//               for (unsigned int i=0; i<(polynomial_degree-1)*(polynomial_degree-1); ++i)
//                 a.push_back(b[4*polynomial_degree+i]);
//             }
//           else
//             {
//               // face is entirely in the interior. get intermediate
//               // points from the relevant manifold object.
//               vertices.resize(4);
//               for (unsigned int i=0; i<4; ++i)
//                 vertices[i] = face->vertex(i);
//               get_intermediate_points (face->get_manifold(), line_support_points, vertices, quad_points);
//               // in 3D, the orientation, flip and rotation of the face might
//               // not match what we expect here, namely the standard
//               // orientation. thus reorder points accordingly. since a Mapping
//               // uses the same shape function as an FE_Q, we can ask a FE_Q to
//               // do the reordering for us.
//               for (unsigned int i=0; i<quad_points.size(); ++i)
//                 a.push_back(quad_points[fe_q->adjust_quad_dof_index_for_face_orientation(i,
//                                         face_orientation,
//                                         face_flip,
//                                         face_rotation)]);
//             }
//         }
//     }
}



template <>
void
MappingManifold<2,3>::
add_quad_support_points(const Triangulation<2,3>::cell_iterator &cell,
                        std::vector<Point<3> >                &a) const
{
  // std::vector<Point<3> > quad_points ((polynomial_degree-1)*(polynomial_degree-1));
  // get_intermediate_points_on_object (cell->get_manifold(), line_support_points,
  //                                    cell, quad_points);
  // for (unsigned int i=0; i<quad_points.size(); ++i)
  //   a.push_back(quad_points[i]);
}



template <int dim, int spacedim>
void
MappingManifold<dim,spacedim>::
add_quad_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &,
                        std::vector<Point<spacedim> > &) const
{
  Assert (false, ExcInternalError());
}



template<int dim, int spacedim>
std::vector<Point<spacedim> >
MappingManifold<dim,spacedim>::
compute_mapping_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &cell) const
{
  // get the vertices first
  std::vector<Point<spacedim> > a(GeometryInfo<dim>::vertices_per_cell);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a[i] = cell->vertex(i);

  if (this->polynomial_degree>1)
    switch (dim)
      {
      case 1:
        add_line_support_points(cell, a);
        break;
      case 2:
        // in 2d, add the points on the four bounding lines to the exterior
        // (outer) points
        add_line_support_points(cell, a);

        // then get the support points on the quad if we are on a
        // manifold, otherwise compute them from the points around it
        if (dim != spacedim)
          add_quad_support_points(cell, a);
        else
          add_weighted_interior_points (support_point_weights_on_quad, a);
        break;

      case 3:
      {
        // in 3d also add the points located on the boundary faces
        add_line_support_points (cell, a);
        add_quad_support_points (cell, a);

        // then compute the interior points
        add_weighted_interior_points (support_point_weights_on_hex, a);
        break;
      }

      default:
        Assert(false, ExcNotImplemented());
        break;
      }

  return a;
}



//--------------------------- Explicit instantiations -----------------------
#include "mapping_q_generic.inst"


DEAL_II_NAMESPACE_CLOSE

template class dealii::MappingManifold<2,2>;
