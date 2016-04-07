// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2016 by the deal.II authors
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

#include <deal.II/base/utilities.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/numerics/vector_tools.h>

#include <numeric>
#include <fstream>



DEAL_II_NAMESPACE_OPEN


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData::InternalData
(const FiniteElement<dim,spacedim> &fe,
 const ComponentMask                mask)
  :
  n_shape_functions (fe.dofs_per_cell),
  mask (mask),
  local_dof_indices(fe.dofs_per_cell),
  local_dof_values(fe.dofs_per_cell)
{}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
std::size_t
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData::memory_consumption () const
{
  Assert (false, ExcNotImplemented());
  return 0;
}



template<int dim, int spacedim, typename DoFHandlerType, typename VectorType>
double &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::shape
(const unsigned int qpoint,
 const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}


template<int dim, int spacedim, typename DoFHandlerType, typename VectorType>
const Tensor<1,dim> &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::derivative
(const unsigned int qpoint,
 const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}



template<int dim, int spacedim, typename DoFHandlerType, typename VectorType>
Tensor<1,dim> &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::derivative
(const unsigned int qpoint,
 const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
const Tensor<2,dim> &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::second_derivative
(const unsigned int qpoint,
 const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_second_derivatives.size()));
  return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
}



template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
Tensor<2,dim> &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::second_derivative
(const unsigned int qpoint,
 const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_second_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_second_derivatives.size()));
  return shape_second_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
const Tensor<3,dim> &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::third_derivative
(const unsigned int qpoint,
 const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_third_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_third_derivatives.size()));
  return shape_third_derivatives [qpoint*n_shape_functions + shape_nr];
}



template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
Tensor<3,dim> &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::third_derivative
(const unsigned int qpoint,
 const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_third_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_third_derivatives.size()));
  return shape_third_derivatives [qpoint*n_shape_functions + shape_nr];
}


template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
const Tensor<4,dim> &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::fourth_derivative
(const unsigned int qpoint,
 const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_fourth_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_fourth_derivatives.size()));
  return shape_fourth_derivatives [qpoint*n_shape_functions + shape_nr];
}



template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
Tensor<4,dim> &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::fourth_derivative
(const unsigned int qpoint,
 const unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_fourth_derivatives.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_fourth_derivatives.size()));
  return shape_fourth_derivatives [qpoint*n_shape_functions + shape_nr];
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::MappingFEField
(const DoFHandlerType            &euler_dof_handler,
 const VectorType    &euler_vector,
 const ComponentMask  mask)
  :
  euler_vector(&euler_vector),
  fe(&euler_dof_handler.get_fe()),
  euler_dof_handler(&euler_dof_handler),
  fe_mask(mask.size() ? mask :
          ComponentMask(fe->get_nonzero_components(0).size(), true)),
  fe_to_real(fe_mask.size(), numbers::invalid_unsigned_int)
{
  unsigned int size = 0;
  for (unsigned int i=0; i<fe_mask.size(); ++i)
    {
      if (fe_mask[i])
        fe_to_real[i] = size++;
    }
  AssertDimension(size,spacedim);
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::MappingFEField
(const MappingFEField<dim,spacedim,VectorType,DoFHandlerType> &mapping)
  :
  euler_vector(mapping.euler_vector),
  fe(mapping.fe),
  euler_dof_handler(mapping.euler_dof_handler),
  fe_mask(mapping.fe_mask),
  fe_to_real(mapping.fe_to_real)
{}



template<int dim, int spacedim, typename DoFHandlerType, typename VectorType>
inline
const double &
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::InternalData::shape
(const unsigned int qpoint,
 const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
         ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0,
                       shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}



template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
bool
MappingFEField<dim,spacedim,DoFHandlerType,VectorType>::preserves_vertex_locations () const
{
  return false;
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
compute_shapes_virtual (const std::vector<Point<dim> >                       &unit_points,
                        typename MappingFEField<dim, spacedim>::InternalData &data) const
{
  const unsigned int n_points=unit_points.size();

  for (unsigned int point=0; point<n_points; ++point)
    {
      if (data.shape_values.size()!=0)
        for (unsigned int i=0; i<data.n_shape_functions; ++i)
          data.shape(point, i) = fe->shape_value(i, unit_points[point]);

      if (data.shape_derivatives.size()!=0)
        for (unsigned int i=0; i<data.n_shape_functions; ++i)
          data.derivative(point, i) = fe->shape_grad(i, unit_points[point]);

      if (data.shape_second_derivatives.size()!=0)
        for (unsigned int i=0; i<data.n_shape_functions; ++i)
          data.second_derivative(point, i) = fe->shape_grad_grad(i, unit_points[point]);

      if (data.shape_third_derivatives.size()!=0)
        for (unsigned int i=0; i<data.n_shape_functions; ++i)
          data.third_derivative(point, i) = fe->shape_3rd_derivative(i, unit_points[point]);

      if (data.shape_fourth_derivatives.size()!=0)
        for (unsigned int i=0; i<data.n_shape_functions; ++i)
          data.fourth_derivative(point, i) = fe->shape_4th_derivative(i, unit_points[point]);
    }
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
UpdateFlags
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::requires_update_flags (const UpdateFlags in) const
{
  // add flags if the respective quantities are necessary to compute
  // what we need. note that some flags appear in both conditions and
  // in subsequent set operations. this leads to some circular
  // logic. the only way to treat this is to iterate. since there are
  // 5 if-clauses in the loop, it will take at most 4 iterations to
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
      // is a Piola transformation, which
      // requires the determinant of the
      // Jacobi matrix of the transformation.
      // Therefore these values have to be
      // updated for each cell.
      if (out & update_contravariant_transformation)
        out |= update_JxW_values;

      if (out & update_normal_vectors)
        out |= update_JxW_values;
    }

  return out;
}




template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::compute_data
(const UpdateFlags      update_flags,
 const Quadrature<dim> &q,
 const unsigned int     n_original_q_points,
 InternalData          &data) const
{
  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values(). use the transitive hull of the required
  // flags
  data.update_each = requires_update_flags(update_flags);

  const unsigned int n_q_points = q.size();

  // see if we need the (transformation) shape function values
  // and/or gradients and resize the necessary arrays
  if (data.update_each & update_quadrature_points)
    data.shape_values.resize(data.n_shape_functions * n_q_points);

  if (data.update_each & (update_covariant_transformation
                          | update_contravariant_transformation
                          | update_JxW_values
                          | update_boundary_forms
                          | update_normal_vectors
                          | update_jacobians
                          | update_jacobian_grads
                          | update_inverse_jacobians))
    data.shape_derivatives.resize(data.n_shape_functions * n_q_points);

  if (data.update_each & update_covariant_transformation)
    data.covariant.resize(n_original_q_points);

  if (data.update_each & update_contravariant_transformation)
    data.contravariant.resize(n_original_q_points);

  if (data.update_each & update_volume_elements)
    data.volume_elements.resize(n_original_q_points);

  if (data.update_each & (update_jacobian_grads | update_jacobian_pushed_forward_grads) )
    data.shape_second_derivatives.resize(data.n_shape_functions * n_q_points);

  if (data.update_each & (update_jacobian_2nd_derivatives | update_jacobian_pushed_forward_2nd_derivatives) )
    data.shape_third_derivatives.resize(data.n_shape_functions * n_q_points);

  if (data.update_each & (update_jacobian_3rd_derivatives | update_jacobian_pushed_forward_3rd_derivatives) )
    data.shape_fourth_derivatives.resize(data.n_shape_functions * n_q_points);

  compute_shapes_virtual (q.get_points(), data);
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::compute_face_data
(const UpdateFlags      update_flags,
 const Quadrature<dim> &q,
 const unsigned int     n_original_q_points,
 InternalData          &data) const
{
  compute_data (update_flags, q, n_original_q_points, data);

  if (dim > 1)
    {
      if (data.update_each & update_boundary_forms)
        {
          data.aux.resize (dim-1, std::vector<Tensor<1,spacedim> > (n_original_q_points));

          // Compute tangentials to the unit cell.
          for (unsigned int i=0; i<data.unit_tangentials.size(); ++i)
            data.unit_tangentials[i].resize (n_original_q_points);

          if (dim==2)
            {
              // ensure a counterclockwise
              // orientation of tangentials
              static const int tangential_orientation[4]= {-1,1,1,-1};
              for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
                {
                  Tensor<1,dim> tang;
                  tang[1-i/2]=tangential_orientation[i];
                  std::fill (data.unit_tangentials[i].begin(),
                             data.unit_tangentials[i].end(),
                             tang);
                }
            }
          else if (dim==3)
            {
              for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
                {
                  Tensor<1,dim> tang1, tang2;

                  const unsigned int nd=
                    GeometryInfo<dim>::unit_normal_direction[i];

                  // first tangential
                  // vector in direction
                  // of the (nd+1)%3 axis
                  // and inverted in case
                  // of unit inward normal
                  tang1[(nd+1)%dim]=GeometryInfo<dim>::unit_normal_orientation[i];
                  // second tangential
                  // vector in direction
                  // of the (nd+2)%3 axis
                  tang2[(nd+2)%dim]=1.;

                  // same unit tangents
                  // for all quadrature
                  // points on this face
                  std::fill (data.unit_tangentials[i].begin(),
                             data.unit_tangentials[i].end(),
                             tang1);
                  std::fill (data.unit_tangentials[GeometryInfo<dim>::faces_per_cell+i].begin(),
                             data.unit_tangentials[GeometryInfo<dim>::faces_per_cell+i].end(),
                             tang2);
                }
            }
        }
    }
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
typename
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData *
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::get_data (const UpdateFlags      update_flags,
    const Quadrature<dim> &quadrature) const
{
  InternalData *data = new InternalData(*fe, fe_mask);
  this->compute_data (update_flags, quadrature,
                      quadrature.size(), *data);
  return data;
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::get_face_data
(const UpdateFlags        update_flags,
 const Quadrature<dim-1> &quadrature) const
{
  InternalData *data = new InternalData(*fe, fe_mask);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_faces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.size(), *data);

  return data;
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::get_subface_data
(const UpdateFlags        update_flags,
 const Quadrature<dim-1> &quadrature) const
{
  InternalData *data = new InternalData(*fe, fe_mask);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_subfaces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.size(), *data);

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
    maybe_compute_q_points (const typename dealii::QProjector<dim>::DataSetDescriptor data_set,
                            const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
                            const FiniteElement<dim, spacedim>   &fe,
                            const ComponentMask                  &fe_mask,
                            const std::vector<unsigned int>      &fe_to_real,
                            std::vector<Point<spacedim> >        &quadrature_points)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_quadrature_points)
        {
          for (unsigned int point=0; point<quadrature_points.size(); ++point)
            {
              Point<spacedim> result;
              const double *shape = &data.shape(point+data_set,0);

              for (unsigned int k=0; k<data.n_shape_functions; ++k)
                {
                  unsigned int comp_k = fe.system_to_component_index(k).first;
                  if (fe_mask[comp_k])
                    result[fe_to_real[comp_k]] += data.local_dof_values[k] * shape[k];
                }

              quadrature_points[point] = result;
            }
        }
    }

    /**
     * Update the co- and contravariant matrices as well as their determinant,
     * for the cell described stored in the data object, but only if the
     * update_flags of the @p data argument indicate so.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_Jacobians (const CellSimilarity::Similarity    cell_similarity,
                            const typename dealii::QProjector<dim>::DataSetDescriptor  data_set,
                            const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
                            const FiniteElement<dim, spacedim> &fe,
                            const ComponentMask                &fe_mask,
                            const std::vector<unsigned int>    &fe_to_real)
    {
      const UpdateFlags update_flags = data.update_each;

      // then Jacobians
      if (update_flags & update_contravariant_transformation)
        {

          // if the current cell is just a translation of the previous one, no
          // need to recompute jacobians...
          if (cell_similarity != CellSimilarity::translation)
            {
              const unsigned int n_q_points = data.contravariant.size();

              Assert (data.n_shape_functions > 0, ExcInternalError());

              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const Tensor<1,dim> *data_derv =
                    &data.derivative(point+data_set, 0);

                  Tensor<1, dim> result[spacedim];

                  for (unsigned int k=0; k<data.n_shape_functions; ++k)
                    {
                      unsigned int comp_k = fe.system_to_component_index(k).first;
                      if (fe_mask[comp_k])
                        result[fe_to_real[comp_k]] += data.local_dof_values[k] * data_derv[k];
                    }

                  // write result into contravariant data
                  for (unsigned int i=0; i<spacedim; ++i)
                    {
                      data.contravariant[point][i] = result[i];
                    }
                }
            }
        }

      if (update_flags & update_covariant_transformation)
        {
          AssertDimension(data.covariant.size(), data.contravariant.size());
          if (cell_similarity != CellSimilarity::translation)
            for (unsigned int point=0; point<data.contravariant.size(); ++point)
              data.covariant[point] = (data.contravariant[point]).covariant_form();
        }

      if (update_flags & update_volume_elements)
        {
          AssertDimension(data.covariant.size(), data.volume_elements.size());
          if (cell_similarity != CellSimilarity::translation)
            for (unsigned int point=0; point<data.contravariant.size(); ++point)
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
    maybe_update_jacobian_grads (const CellSimilarity::Similarity              cell_similarity,
                                 const typename dealii::QProjector<dim>::DataSetDescriptor data_set,
                                 const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
                                 const FiniteElement<dim, spacedim>           &fe,
                                 const ComponentMask                          &fe_mask,
                                 const std::vector<unsigned int>              &fe_to_real,
                                 std::vector<DerivativeForm<2,dim,spacedim> > &jacobian_grads)
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

                  DerivativeForm<2,dim,spacedim> result;

                  for (unsigned int k=0; k<data.n_shape_functions; ++k)
                    {
                      unsigned int comp_k = fe.system_to_component_index(k).first;
                      if (fe_mask[comp_k])
                        for (unsigned int j=0; j<dim; ++j)
                          for (unsigned int l=0; l<dim; ++l)
                            result[fe_to_real[comp_k]][j][l] += (second[k][j][l]
                                                                 * data.local_dof_values[k]);
                    }

                  // never touch any data for j=dim in case dim<spacedim, so
                  // it will always be zero as it was initialized
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
    maybe_update_jacobian_pushed_forward_grads (
      const CellSimilarity::Similarity              cell_similarity,
      const typename dealii::QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
      const FiniteElement<dim, spacedim>           &fe,
      const ComponentMask                          &fe_mask,
      const std::vector<unsigned int>              &fe_to_real,
      std::vector<Tensor<3,spacedim> >             &jacobian_pushed_forward_grads )
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

                  DerivativeForm<2,dim,spacedim> result;

                  for (unsigned int k=0; k<data.n_shape_functions; ++k)
                    {
                      unsigned int comp_k = fe.system_to_component_index(k).first;
                      if (fe_mask[comp_k])
                        for (unsigned int j=0; j<dim; ++j)
                          for (unsigned int l=0; l<dim; ++l)
                            result[fe_to_real[comp_k]][j][l] += (second[k][j][l]
                                                                 * data.local_dof_values[k]);
                    }

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
     * Update the third derivative of the transformation from unit to real
     * cell, the Jacobian hessians.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_2nd_derivatives (const CellSimilarity::Similarity              cell_similarity,
                                           const typename dealii::QProjector<dim>::DataSetDescriptor data_set,
                                           const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
                                           const FiniteElement<dim, spacedim>           &fe,
                                           const ComponentMask                          &fe_mask,
                                           const std::vector<unsigned int>              &fe_to_real,
                                           std::vector<DerivativeForm<3,dim,spacedim> > &jacobian_2nd_derivatives)
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

                  DerivativeForm<3,dim,spacedim> result;

                  for (unsigned int k=0; k<data.n_shape_functions; ++k)
                    {
                      unsigned int comp_k = fe.system_to_component_index(k).first;
                      if (fe_mask[comp_k])
                        for (unsigned int j=0; j<dim; ++j)
                          for (unsigned int l=0; l<dim; ++l)
                            for (unsigned int m=0; m<dim; ++m)
                              result[fe_to_real[comp_k]][j][l][m] += (third[k][j][l][m]
                                                                      * data.local_dof_values[k]);
                    }

                  // never touch any data for j=dim in case dim<spacedim, so
                  // it will always be zero as it was initialized
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
     * Update the third derivative of the transformation from unit to real cell,
     * the Jacobian hessians, pushed forward to the real cell coordinates.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_pushed_forward_2nd_derivatives (
      const CellSimilarity::Similarity              cell_similarity,
      const typename dealii::QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
      const FiniteElement<dim, spacedim>           &fe,
      const ComponentMask                          &fe_mask,
      const std::vector<unsigned int>              &fe_to_real,
      std::vector<Tensor<4,spacedim> >             &jacobian_pushed_forward_2nd_derivatives )
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

                  DerivativeForm<3,dim,spacedim> result;

                  for (unsigned int k=0; k<data.n_shape_functions; ++k)
                    {
                      unsigned int comp_k = fe.system_to_component_index(k).first;
                      if (fe_mask[comp_k])
                        for (unsigned int j=0; j<dim; ++j)
                          for (unsigned int l=0; l<dim; ++l)
                            for (unsigned int m=0; m<dim; ++m)
                              result[fe_to_real[comp_k]][j][l][m] += (third[k][j][l][m]
                                                                      * data.local_dof_values[k]);
                    }

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
  }

  /**
   * Update the fourth derivative of the transformation from unit to real
   * cell, the Jacobian hessian gradients.
   *
   * Skip the computation if possible as indicated by the first argument.
   */
  template <int dim, int spacedim>
  void
  maybe_update_jacobian_3rd_derivatives (const CellSimilarity::Similarity              cell_similarity,
                                         const typename dealii::QProjector<dim>::DataSetDescriptor data_set,
                                         const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
                                         const FiniteElement<dim, spacedim>           &fe,
                                         const ComponentMask                          &fe_mask,
                                         const std::vector<unsigned int>              &fe_to_real,
                                         std::vector<DerivativeForm<4,dim,spacedim> > &jacobian_3rd_derivatives)
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

                DerivativeForm<4,dim,spacedim> result;

                for (unsigned int k=0; k<data.n_shape_functions; ++k)
                  {
                    unsigned int comp_k = fe.system_to_component_index(k).first;
                    if (fe_mask[comp_k])
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          for (unsigned int m=0; m<dim; ++m)
                            for (unsigned int n=0; n<dim; ++n)
                              result[fe_to_real[comp_k]][j][l][m][n] += (fourth[k][j][l][m][n]
                                                                         * data.local_dof_values[k]);
                  }

                // never touch any data for j,l,m,n=dim in case dim<spacedim, so
                // it will always be zero as it was initialized
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
   * Update the fourth derivative of the transformation from unit to real cell,
   * the Jacobian hessian gradients, pushed forward to the real cell
   * coordinates.
   *
   * Skip the computation if possible as indicated by the first argument.
   */
  template <int dim, int spacedim>
  void
  maybe_update_jacobian_pushed_forward_3rd_derivatives (
    const CellSimilarity::Similarity              cell_similarity,
    const typename dealii::QProjector<dim>::DataSetDescriptor data_set,
    const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
    const FiniteElement<dim, spacedim>           &fe,
    const ComponentMask                          &fe_mask,
    const std::vector<unsigned int>              &fe_to_real,
    std::vector<Tensor<5,spacedim> >             &jacobian_pushed_forward_3rd_derivatives )
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

                DerivativeForm<4,dim,spacedim> result;

                for (unsigned int k=0; k<data.n_shape_functions; ++k)
                  {
                    unsigned int comp_k = fe.system_to_component_index(k).first;
                    if (fe_mask[comp_k])
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          for (unsigned int m=0; m<dim; ++m)
                            for (unsigned int n=0; n<dim; ++n)
                              result[fe_to_real[comp_k]][j][l][m][n]
                              += (fourth[k][j][l][m][n]
                                  * data.local_dof_values[k]);
                  }

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


  /**
   * Depending on what information is called for in the update flags of the
   * @p data object, compute the various pieces of information that is
   * required by the fill_fe_face_values() and fill_fe_subface_values()
   * functions.  This function simply unifies the work that would be done by
   * those two functions.
   *
   * The resulting data is put into the @p output_data argument.
   */
  template <int dim, int spacedim>
  void
  maybe_compute_face_data (const dealii::MappingFEField<dim,spacedim> &mapping,
                           const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                           const unsigned int               face_no,
                           const unsigned int               subface_no,
                           const std::vector<double>        &weights,
                           const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
                           internal::FEValues::MappingRelatedData<dim,spacedim>         &output_data)
  {
    const UpdateFlags update_flags = data.update_each;

    if (update_flags & update_boundary_forms)
      {
        const unsigned int n_q_points = output_data.boundary_forms.size();
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
                  // can still compute the boundary form by simply looking
                  // at the number of the face
                  output_data.boundary_forms[i][0] = (face_no == 0 ?
                                                      -1 : +1);
                  break;
                case 2:
                  output_data.boundary_forms[i] = cross_product_2d(data.aux[0][i]);
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

        if (update_flags & (update_normal_vectors | update_JxW_values))
          for (unsigned int i=0; i<output_data.boundary_forms.size(); ++i)
            {
              if (update_flags & update_JxW_values)
                {
                  output_data.JxW_values[i] = output_data.boundary_forms[i].norm() * weights[i];

                  if (subface_no != numbers::invalid_unsigned_int)
                    {
                      const double area_ratio=GeometryInfo<dim>::subface_ratio(
                                                cell->subface_case(face_no), subface_no);
                      output_data.JxW_values[i] *= area_ratio;
                    }
                }

              if (update_flags & update_normal_vectors)
                output_data.normal_vectors[i] = Point<spacedim>(output_data.boundary_forms[i] / output_data.boundary_forms[i].norm());
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
   * Do the work of MappingFEField::fill_fe_face_values() and
   * MappingFEField::fill_fe_subface_values() in a generic way, using the
   * 'data_set' to differentiate whether we will work on a face (and if so,
   * which one) or subface.
   */
  template<int dim, int spacedim>
  void
  do_fill_fe_face_values (const dealii::MappingFEField<dim,spacedim>                        &mapping,
                          const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                          const unsigned int                                                 face_no,
                          const unsigned int                                                 subface_no,
                          const typename dealii::QProjector<dim>::DataSetDescriptor          data_set,
                          const Quadrature<dim-1>                                           &quadrature,
                          const typename dealii::MappingFEField<dim,spacedim>::InternalData &data,
                          const FiniteElement<dim, spacedim>                                &fe,
                          const ComponentMask                                               &fe_mask,
                          const std::vector<unsigned int>                                   &fe_to_real,
                          internal::FEValues::MappingRelatedData<dim,spacedim>              &output_data)
  {
    maybe_compute_q_points<dim,spacedim> (data_set,
                                          data,
                                          fe, fe_mask, fe_to_real,
                                          output_data.quadrature_points);

    maybe_update_Jacobians<dim,spacedim> (CellSimilarity::none,
                                          data_set,
                                          data,
                                          fe, fe_mask, fe_to_real);

    maybe_update_jacobian_grads<dim,spacedim> (CellSimilarity::none,
                                               data_set,
                                               data,
                                               fe, fe_mask, fe_to_real,
                                               output_data.jacobian_grads);

    maybe_update_jacobian_pushed_forward_grads<dim,spacedim> (CellSimilarity::none,
                                                              data_set,
                                                              data,
                                                              fe, fe_mask, fe_to_real,
                                                              output_data.jacobian_pushed_forward_grads);

    maybe_update_jacobian_2nd_derivatives<dim,spacedim> (CellSimilarity::none,
                                                         data_set,
                                                         data,
                                                         fe, fe_mask, fe_to_real,
                                                         output_data.jacobian_2nd_derivatives);

    maybe_update_jacobian_pushed_forward_2nd_derivatives<dim,spacedim> (CellSimilarity::none,
        data_set,
        data,
        fe, fe_mask, fe_to_real,
        output_data.jacobian_pushed_forward_2nd_derivatives);

    maybe_update_jacobian_3rd_derivatives<dim,spacedim> (CellSimilarity::none,
                                                         data_set,
                                                         data,
                                                         fe, fe_mask, fe_to_real,
                                                         output_data.jacobian_3rd_derivatives);

    maybe_update_jacobian_pushed_forward_3rd_derivatives<dim,spacedim> (CellSimilarity::none,
        data_set,
        data,
        fe, fe_mask, fe_to_real,
        output_data.jacobian_pushed_forward_3rd_derivatives);

    maybe_compute_face_data (mapping,
                             cell, face_no, subface_no,
                             quadrature.get_weights(), data,
                             output_data);
  }
}


// Note that the CellSimilarity flag is modifiable, since MappingFEField can need to
// recalculate data even when cells are similar.
template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
CellSimilarity::Similarity
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                const CellSimilarity::Similarity                           cell_similarity,
                const Quadrature<dim>                                     &quadrature,
                const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0, ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  const unsigned int n_q_points=quadrature.size();
  const  CellSimilarity::Similarity updated_cell_similarity
    = (get_degree() == 1
       ?
       cell_similarity
       :
       CellSimilarity::invalid_next_cell);

  update_internal_dofs(cell, data);

  internal::maybe_compute_q_points(QProjector<dim>::DataSetDescriptor::cell (),
                                   data, *fe, fe_mask, fe_to_real,
                                   output_data.quadrature_points);

  internal::maybe_update_Jacobians(cell_similarity,
                                   QProjector<dim>::DataSetDescriptor::cell (),
                                   data, *fe, fe_mask, fe_to_real);

  const UpdateFlags update_flags = data.update_each;
  const std::vector<double> &weights=quadrature.get_weights();

  // Multiply quadrature weights by absolute value of Jacobian determinants or
  // the area element g=sqrt(DX^t DX) in case of codim > 0

  if (update_flags & (update_normal_vectors | update_JxW_values))
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

                output_data.JxW_values[point] = sqrt(determinant(G)) * weights[point];

                if (cell_similarity == CellSimilarity::inverted_translation)
                  {
                    // we only need to flip the normal
                    if (update_flags & update_normal_vectors)
                      output_data.normal_vectors[point] *= -1.;
                  }
                else
                  {
                    if (update_flags & update_normal_vectors)
                      {
                        Assert (spacedim - dim == 1,
                                ExcMessage("There is no cell normal in codim 2."));

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

  // calculate derivatives of the Jacobians
  internal::maybe_update_jacobian_grads(cell_similarity,
                                        QProjector<dim>::DataSetDescriptor::cell(),
                                        data, *fe, fe_mask, fe_to_real,
                                        output_data.jacobian_grads);

  // calculate derivatives of the Jacobians pushed forward to real cell coordinates
  internal::maybe_update_jacobian_pushed_forward_grads(cell_similarity,
                                                       QProjector<dim>::DataSetDescriptor::cell(),
                                                       data, *fe, fe_mask, fe_to_real,
                                                       output_data.jacobian_pushed_forward_grads);

  // calculate hessians of the Jacobians
  internal::maybe_update_jacobian_2nd_derivatives(cell_similarity,
                                                  QProjector<dim>::DataSetDescriptor::cell(),
                                                  data, *fe, fe_mask, fe_to_real,
                                                  output_data.jacobian_2nd_derivatives);

  // calculate hessians of the Jacobians pushed forward to real cell coordinates
  internal::maybe_update_jacobian_pushed_forward_2nd_derivatives(cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell(),
      data, *fe, fe_mask, fe_to_real,
      output_data.jacobian_pushed_forward_2nd_derivatives);

  // calculate gradients of the hessians of the Jacobians
  internal::maybe_update_jacobian_3rd_derivatives(cell_similarity,
                                                  QProjector<dim>::DataSetDescriptor::cell(),
                                                  data, *fe, fe_mask, fe_to_real,
                                                  output_data.jacobian_3rd_derivatives);

  // calculate gradients of the hessians of the Jacobians pushed forward to real
  // cell coordinates
  internal::maybe_update_jacobian_pushed_forward_3rd_derivatives(cell_similarity,
      QProjector<dim>::DataSetDescriptor::cell(),
      data, *fe, fe_mask, fe_to_real,
      output_data.jacobian_pushed_forward_3rd_derivatives);

  return updated_cell_similarity;
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int                                         face_no,
                     const Quadrature<dim-1>                                   &quadrature,
                     const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                     internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  update_internal_dofs(cell, data);

  internal::do_fill_fe_face_values (*this,
                                    cell, face_no, numbers::invalid_unsigned_int,
                                    QProjector<dim>::DataSetDescriptor::
                                    face (face_no,
                                          cell->face_orientation(face_no),
                                          cell->face_flip(face_no),
                                          cell->face_rotation(face_no),
                                          quadrature.size()),
                                    quadrature,
                                    data,
                                    *fe, fe_mask, fe_to_real,
                                    output_data);
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                        const unsigned int                                         face_no,
                        const unsigned int                                         subface_no,
                        const Quadrature<dim-1>                                   &quadrature,
                        const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                        internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  update_internal_dofs(cell, data);

  internal::do_fill_fe_face_values (*this,
                                    cell, face_no, numbers::invalid_unsigned_int,
                                    QProjector<dim>::DataSetDescriptor::
                                    subface (face_no, subface_no,
                                             cell->face_orientation(face_no),
                                             cell->face_flip(face_no),
                                             cell->face_rotation(face_no),
                                             quadrature.size(),
                                             cell->subface_case(face_no)),
                                    quadrature,
                                    data,
                                    *fe, fe_mask, fe_to_real,
                                    output_data);
}


namespace
{
  template<int dim, int spacedim, int rank, typename VectorType, typename DoFHandlerType>
  void
  transform_fields(const ArrayView<const Tensor<rank,dim> >                &input,
                   const MappingType                                        mapping_type,
                   const typename Mapping<dim,spacedim>::InternalDataBase  &mapping_data,
                   const ArrayView<Tensor<rank,spacedim> >                 &output)
  {
    AssertDimension (input.size(), output.size());
    Assert ((dynamic_cast<const typename MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData *>(&mapping_data) != 0),
            ExcInternalError());
    const typename MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData
    &data = static_cast<const typename MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData &>(mapping_data);

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
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

        for (unsigned int i=0; i<output.size(); ++i)
          output[i] = apply_transformation(data.covariant[i], input[i]);

        return;
      }

      default:
        Assert(false, ExcNotImplemented());
      }
  }


  template<int dim, int spacedim, int rank, typename VectorType, typename DoFHandlerType>
  void
  transform_differential_forms
  (const ArrayView<const DerivativeForm<rank, dim,spacedim> >  &input,
   const MappingType                                            mapping_type,
   const typename Mapping<dim,spacedim>::InternalDataBase      &mapping_data,
   const ArrayView<Tensor<rank+1, spacedim> >                  &output)
  {

    AssertDimension (input.size(), output.size());
    Assert ((dynamic_cast<const typename MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData *>(&mapping_data) != 0),
            ExcInternalError());
    const typename MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData
    &data = static_cast<const typename MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::InternalData &>(mapping_data);

    switch (mapping_type)
      {
      case mapping_covariant:
      {
        Assert (data.update_each & update_contravariant_transformation,
                typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

        for (unsigned int i=0; i<output.size(); ++i)
          output[i] = apply_transformation(data.covariant[i], input[i]);

        return;
      }
      default:
        Assert(false, ExcNotImplemented());
      }

  }
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
transform (const ArrayView<const Tensor<1,dim> >                  &input,
           const MappingType                                       mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           const ArrayView<Tensor<1,spacedim> >                   &output) const
{
  AssertDimension (input.size(), output.size());

  transform_fields<dim,spacedim,1,VectorType,DoFHandlerType>(input, mapping_type, mapping_data, output);
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
transform (const ArrayView<const DerivativeForm<1, dim ,spacedim> > &input,
           const MappingType                                         mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase   &mapping_data,
           const ArrayView<Tensor<2,spacedim> >                     &output) const
{
  AssertDimension (input.size(), output.size());

  transform_differential_forms<dim,spacedim,1,VectorType,DoFHandlerType>(input, mapping_type, mapping_data, output);
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
transform (const ArrayView<const Tensor<2, dim> >                 &input,
           const MappingType                                       ,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           const ArrayView<Tensor<2,spacedim> >                   &output) const
{
  (void)input;
  (void)output;
  (void)mapping_data;
  AssertDimension (input.size(), output.size());

  AssertThrow(false, ExcNotImplemented());
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
transform (const ArrayView<const DerivativeForm<2, dim, spacedim> >  &input,
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
            for (unsigned int k=0; k<spacedim; ++k)
              {
                output[q][i][j][k] = data.covariant[q][j][0]
                                     * data.covariant[q][k][0]
                                     * input[q][i][0][0];
                for (unsigned int J=0; J<dim; ++J)
                  {
                    const unsigned int K0 = (0==J)? 1 : 0;
                    for (unsigned int K=K0; K<dim; ++K)
                      output[q][i][j][k] += data.covariant[q][j][J]
                                            * data.covariant[q][k][K]
                                            * input[q][i][J][K];
                  }

              }
      return;
    }

    default:
      Assert(false, ExcNotImplemented());
    }

}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
transform (const ArrayView<const Tensor<3,dim> >                  &input,
           const MappingType                                     /*mapping_type*/,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           const ArrayView<Tensor<3,spacedim> >                   &output) const
{

  (void)input;
  (void)output;
  (void)mapping_data;
  AssertDimension (input.size(), output.size());

  AssertThrow(false, ExcNotImplemented());

}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
Point<spacedim>
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<dim>                                          &p) const
{
//  Use the get_data function to create an InternalData with data vectors of
//  the right size and transformation shape values already computed at point
//  p.
  const Quadrature<dim> point_quadrature(p);
  std_cxx11::unique_ptr<InternalData> mdata (get_data(update_quadrature_points | update_jacobians,
                                                      point_quadrature));

  update_internal_dofs(cell, *mdata);

  return do_transform_unit_to_real_cell(*mdata);
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
Point<spacedim>
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
do_transform_unit_to_real_cell (const InternalData &data) const
{
  Point<spacedim> p_real;

  for (unsigned int i=0; i<data.n_shape_functions; ++i)
    {
      unsigned int comp_i = fe->system_to_component_index(i).first;
      if (fe_mask[comp_i])
        p_real[fe_to_real[comp_i]] += data.local_dof_values[i] * data.shape(0,i);
    }

  return p_real;
}



template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
Point<dim>
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<spacedim>                                     &p) const
{
  // first a Newton iteration based on the real mapping. It uses the center
  // point of the cell as a starting point
  Point<dim> initial_p_unit;
  try
    {
      initial_p_unit
        = StaticMappingQ1<dim,spacedim>::mapping.transform_real_to_unit_cell(cell, p);
    }
  catch (const typename Mapping<dim,spacedim>::ExcTransformationFailed &)
    {
      // mirror the conditions of the code below to determine if we need to
      // use an arbitrary starting point or if we just need to rethrow the
      // exception
      for (unsigned int d=0; d<dim; ++d)
        initial_p_unit[d] = 0.5;
    }

  initial_p_unit = GeometryInfo<dim>::project_to_unit_cell(initial_p_unit);

  // for (unsigned int d=0; d<dim; ++d)
  //   initial_p_unit[d] = 0.;

  const Quadrature<dim> point_quadrature(initial_p_unit);

  UpdateFlags update_flags = update_quadrature_points | update_jacobians;
  if (spacedim>dim)
    update_flags |= update_jacobian_grads;
  std_cxx11::unique_ptr<InternalData>
  mdata (get_data(update_flags,point_quadrature));

  update_internal_dofs(cell, *mdata);

  return do_transform_real_to_unit_cell(cell, p, initial_p_unit, *mdata);

}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
Point<dim>
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::
do_transform_real_to_unit_cell
(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
 const Point<spacedim>                                     &p,
 const Point<dim>                                          &initial_p_unit,
 InternalData                                              &mdata) const
{
  const unsigned int n_shapes=mdata.shape_values.size();
  (void)n_shapes;
  Assert(n_shapes!=0, ExcInternalError());
  AssertDimension (mdata.shape_derivatives.size(), n_shapes);


  // Newton iteration to solve
  // f(x)=p(x)-p=0
  // x_{n+1}=x_n-[f'(x)]^{-1}f(x)
  // The start value was set to be the
  // linear approximation to the cell
  // The shape values and derivatives
  // of the mapping at this point are
  // previously computed.
  // f(x)
  Point<dim> p_unit = initial_p_unit;
  Point<dim> f;
  compute_shapes_virtual(std::vector<Point<dim> > (1, p_unit), mdata);
  Point<spacedim> p_real(do_transform_unit_to_real_cell(mdata));
  Tensor<1,spacedim> p_minus_F = p - p_real;
  const double eps = 1.e-12*cell->diameter();
  const unsigned int newton_iteration_limit = 20;
  unsigned int newton_iteration=0;
  while (p_minus_F.norm_square() > eps*eps)
    {
      // f'(x)
      Point<spacedim>  DF[dim];
      Tensor<2,dim>  df;
      for (unsigned int k=0; k<mdata.n_shape_functions; ++k)
        {
          const Tensor<1,dim> &grad_k = mdata.derivative(0,k);
          unsigned int comp_k = fe->system_to_component_index(k).first;
          if (fe_mask[comp_k])
            for (unsigned int j=0; j<dim; ++j)
              DF[j][fe_to_real[comp_k]] += mdata.local_dof_values[k] * grad_k[j];
        }
      for (unsigned int j=0; j<dim; ++j)
        {
          f[j] = DF[j] * p_minus_F;
          for (unsigned int l=0; l<dim; ++l)
            df[j][l] = -DF[j] * DF[l];
        }
      // Solve  [f'(x)]d=f(x)
      const Tensor<1, dim> delta =
        invert(df) * static_cast<const Tensor<1, dim> &>(f);
      // do a line search
      double step_length = 1;
      do
        {
          // update of p_unit. The
          // spacedimth component of
          // transformed point is simply
          // ignored in codimension one
          // case. When this component is
          // not zero, then we are
          // projecting the point to the
          // surface or curve identified
          // by the cell.
          Point<dim> p_unit_trial = p_unit;
          for (unsigned int i=0; i<dim; ++i)
            p_unit_trial[i] -= step_length * delta[i];
          // shape values and derivatives
          // at new p_unit point
          compute_shapes_virtual(std::vector<Point<dim> > (1, p_unit_trial), mdata);
          // f(x)
          Point<spacedim> p_real_trial = do_transform_unit_to_real_cell(mdata);
          const Tensor<1,spacedim> f_trial = p - p_real_trial;
          // see if we are making progress with the current step length
          // and if not, reduce it by a factor of two and try again
          if (f_trial.norm() < p_minus_F.norm())
            {
              p_real = p_real_trial;
              p_unit = p_unit_trial;
              p_minus_F = f_trial;
              break;
            }
          else if (step_length > 0.05)
            step_length /= 2;
          else
            goto failure;
        }
      while (true);
      ++newton_iteration;
      if (newton_iteration > newton_iteration_limit)
        goto failure;
    }
  return p_unit;
  // if we get to the following label, then we have either run out
  // of Newton iterations, or the line search has not converged.
  // in either case, we need to give up, so throw an exception that
  // can then be caught
failure:
  AssertThrow (false, (typename Mapping<dim,spacedim>::ExcTransformationFailed()));
  // ...the compiler wants us to return something, though we can
  // of course never get here...
  return Point<dim>();
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
unsigned int
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::get_degree() const
{
  return fe->degree;
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
ComponentMask
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::get_component_mask() const
{
  return this->fe_mask;
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
Mapping<dim,spacedim> *
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::clone () const
{
  return new MappingFEField<dim,spacedim,VectorType,DoFHandlerType>(*this);
}


template<int dim, int spacedim, typename VectorType, typename DoFHandlerType>
void
MappingFEField<dim,spacedim,VectorType,DoFHandlerType>::update_internal_dofs
(const typename Triangulation<dim,spacedim>::cell_iterator  &cell,
 const typename MappingFEField<dim, spacedim>::InternalData &data) const
{
  Assert(euler_dof_handler != 0, ExcMessage("euler_dof_handler is empty"));

  typename DoFHandlerType::cell_iterator dof_cell(*cell, euler_dof_handler);
  Assert (dof_cell->active() == true, ExcInactiveCell());

  dof_cell->get_dof_indices(data.local_dof_indices);

  for (unsigned int i=0; i<data.local_dof_values.size(); ++i)
    data.local_dof_values[i] = (*euler_vector)(data.local_dof_indices[i]);
}

// explicit instantiations
#include "mapping_fe_field.inst"


DEAL_II_NAMESPACE_CLOSE
