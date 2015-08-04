// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
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


template<int dim, int spacedim, class VECTOR, class DH>
MappingFEField<dim,spacedim,VECTOR,DH>::InternalData::InternalData (const FiniteElement<dim,spacedim> &fe,
    const ComponentMask mask)
  :
  n_shape_functions (fe.dofs_per_cell),
  mask (mask),
  local_dof_indices(fe.dofs_per_cell),
  local_dof_values(fe.dofs_per_cell)
{}



template<int dim, int spacedim, class VECTOR, class DH>
std::size_t
MappingFEField<dim,spacedim,VECTOR,DH>::InternalData::memory_consumption () const
{
  Assert (false, ExcNotImplemented());
  return 0;
}


template<int dim, int spacedim, class VECTOR, class DH>
MappingFEField<dim,spacedim,VECTOR,DH>::MappingFEField (const DH      &euler_dof_handler,
                                                        const VECTOR  &euler_vector,
                                                        const ComponentMask mask)
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


template<int dim, int spacedim, class VECTOR, class DH>
MappingFEField<dim,spacedim,VECTOR,DH>::MappingFEField (const MappingFEField<dim,spacedim,VECTOR,DH> &mapping)
  :
  euler_vector(mapping.euler_vector),
  fe(mapping.fe),
  euler_dof_handler(mapping.euler_dof_handler),
  fe_mask(mapping.fe_mask),
  fe_to_real(mapping.fe_to_real)
{}


template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::compute_shapes_virtual (
  const std::vector<Point<dim> > &unit_points,
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
    }
}


template<int dim, int spacedim, class VECTOR, class DH>
UpdateFlags
MappingFEField<dim,spacedim,VECTOR,DH>::update_once (const UpdateFlags in) const
{
  UpdateFlags out = UpdateFlags(in & (update_transformation_values
                                      | update_transformation_gradients));

  // Shape function values
  if (in & update_quadrature_points)
    out |= update_transformation_values;

  // Shape function gradients
  if (in & (update_covariant_transformation
            | update_contravariant_transformation
            | update_JxW_values
            | update_boundary_forms
            | update_normal_vectors
            | update_jacobians
            | update_jacobian_grads
            | update_inverse_jacobians))
    out |= update_transformation_gradients;

  return out;
}



template<int dim, int spacedim, class VECTOR, class DH>
UpdateFlags
MappingFEField<dim,spacedim,VECTOR,DH>::update_each (const UpdateFlags in) const
{
  // Select flags of concern for the
  // transformation.
  UpdateFlags out = UpdateFlags(in & (update_quadrature_points
                                      | update_covariant_transformation
                                      | update_contravariant_transformation
                                      | update_JxW_values
                                      | update_boundary_forms
                                      | update_normal_vectors
                                      | update_volume_elements
                                      | update_jacobians
                                      | update_jacobian_grads
                                      | update_inverse_jacobians));

  // add flags if the respective
  // quantities are necessary to
  // compute what we need. note that
  // some flags appear in both
  // conditions and in subsequents
  // set operations. this leads to
  // some circular logic. the only
  // way to treat this is to
  // iterate. since there are 5
  // if-clauses in the loop, it will
  // take at most 4 iterations to
  // converge. do them:
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

      if (out & (update_inverse_jacobians))
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




template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::compute_data (const UpdateFlags      update_flags,
                                                      const Quadrature<dim>  &q,
                                                      const unsigned int     n_original_q_points,
                                                      InternalData           &data) const
{
  const unsigned int n_q_points = q.size();

  data.update_once = update_once(update_flags);
  data.update_each = update_each(update_flags);
  data.update_flags = data.update_once | data.update_each;

  const UpdateFlags flags(data.update_flags);

  if (flags & update_transformation_values)
    data.shape_values.resize(data.n_shape_functions * n_q_points);

  if (flags & update_transformation_gradients)
    data.shape_derivatives.resize(data.n_shape_functions * n_q_points);

  if (flags & update_covariant_transformation)
    data.covariant.resize(n_original_q_points);

  if (flags & update_contravariant_transformation)
    data.contravariant.resize(n_original_q_points);

  if (flags & update_volume_elements)
    data.volume_elements.resize(n_original_q_points);

  if (flags & update_jacobian_grads)
    data.shape_second_derivatives.resize(data.n_shape_functions * n_q_points);

  compute_shapes_virtual (q.get_points(), data);
}


template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::compute_face_data (const UpdateFlags update_flags,
                                                           const Quadrature<dim> &q,
                                                           const unsigned int n_original_q_points,
                                                           InternalData &data) const
{
  compute_data (update_flags, q, n_original_q_points, data);

  if (dim > 1)
    {
      if (data.update_flags & update_boundary_forms)
        {
          data.aux.resize (dim-1, std::vector<Tensor<1,spacedim> > (n_original_q_points));

          // Compute tangentials to the
          // unit cell.
          const unsigned int nfaces = GeometryInfo<dim>::faces_per_cell;
          data.unit_tangentials.resize (nfaces*(dim-1),
                                        std::vector<Tensor<1,dim> > (n_original_q_points));
          if (dim==2)
            {
              // ensure a counterclock wise
              // orientation of tangentials
              static const int tangential_orientation[4]= {-1,1,1,-1};
              for (unsigned int i=0; i<nfaces; ++i)
                {
                  Tensor<1,dim> tang;
                  tang[1-i/2]=tangential_orientation[i];
                  std::fill (data.unit_tangentials[i].begin(),
                             data.unit_tangentials[i].end(), tang);
                }
            }
          else if (dim==3)
            {
              for (unsigned int i=0; i<nfaces; ++i)
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
                             data.unit_tangentials[i].end(), tang1);
                  std::fill (data.unit_tangentials[nfaces+i].begin(),
                             data.unit_tangentials[nfaces+i].end(), tang2);
                }
            }
        }
    }
}


template<int dim, int spacedim, class VECTOR, class DH>
typename
MappingFEField<dim,spacedim,VECTOR,DH>::InternalData *
MappingFEField<dim,spacedim,VECTOR,DH>::get_data (const UpdateFlags update_flags,
                                                  const Quadrature<dim> &quadrature) const
{
  InternalData *data = new InternalData(*fe, fe_mask);
  this->compute_data (update_flags, quadrature,
                      quadrature.size(), *data);
  return data;
}



template<int dim, int spacedim, class VECTOR, class DH>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingFEField<dim,spacedim,VECTOR,DH>::get_face_data (const UpdateFlags update_flags,
                                                       const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(*fe, fe_mask);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_faces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.size(), *data);

  return data;
}


template<int dim, int spacedim, class VECTOR, class DH>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingFEField<dim,spacedim,VECTOR,DH>::get_subface_data (const UpdateFlags update_flags,
                                                          const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(*fe, fe_mask);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_subfaces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.size(), *data);

  return data;
}


// Note that the CellSimilarity flag is modifyable, since MappingFEField can need to
// recalculate data even when cells are similar.
template<int dim, int spacedim, class VECTOR, class DH>
CellSimilarity::Similarity
MappingFEField<dim,spacedim,VECTOR,DH>::
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

  compute_fill (cell, n_q_points, QProjector<dim>::DataSetDescriptor::cell (),
                updated_cell_similarity,
                data,
                output_data.quadrature_points);

  const UpdateFlags update_flags(data.current_update_flags());
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
                          cross_product(output_data.normal_vectors[point], -DX_t[0]);
                        else //dim == 2
                          cross_product(output_data.normal_vectors[point],DX_t[0],DX_t[1]);

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


  // calculate values of the derivatives of the Jacobians. do it here, since
  // we only do it for cells, not faces.
  if (update_flags & update_jacobian_grads)
    {
      AssertDimension (output_data.jacobian_grads.size(), n_q_points);

      if (cell_similarity != CellSimilarity::translation)
        {
          std::fill(output_data.jacobian_grads.begin(),
                    output_data.jacobian_grads.end(),
                    DerivativeForm<2,dim,spacedim>());

          const unsigned int data_set = QProjector<dim>::DataSetDescriptor::cell();

          for (unsigned int point=0; point<n_q_points; ++point)
            {
              const Tensor<2,dim> *second =
                &data.second_derivative(point+data_set, 0);

              DerivativeForm<2,dim,spacedim> result;

              for (unsigned int k=0; k<data.n_shape_functions; ++k)
                {
                  unsigned int comp_k = fe->system_to_component_index(k).first;
                  if (fe_mask[comp_k])
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        result[fe_to_real[comp_k]][j][l] += (second[k][j][l]
                                                             * data.local_dof_values[k]);
                }

              // never touch any data for j=dim in case dim<spacedim, so it
              // will always be zero as it was initialized
              for (unsigned int i=0; i<spacedim; ++i)
                for (unsigned int j=0; j<dim; ++j)
                  for (unsigned int l=0; l<dim; ++l)
                    output_data.jacobian_grads[point][i][j][l] = result[i][j][l];
            }
        }
    }


  // copy values from InternalData to vector given by reference
  if (update_flags & update_inverse_jacobians)
    {
      AssertDimension (output_data.inverse_jacobians.size(), n_q_points);
      if (cell_similarity != CellSimilarity::translation)
        for (unsigned int point=0; point<n_q_points; ++point)
          output_data.inverse_jacobians[point] = data.covariant[point].transpose();
    }

  return updated_cell_similarity;
}



template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::
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

  const unsigned int n_q_points=quadrature.size();
  this->compute_fill_face (cell, face_no, numbers::invalid_unsigned_int,
                           n_q_points,
                           QProjector<dim>::DataSetDescriptor::
                           face (face_no,
                                 cell->face_orientation(face_no),
                                 cell->face_flip(face_no),
                                 cell->face_rotation(face_no),
                                 n_q_points),
                           quadrature.get_weights(),
                           data,
                           output_data.quadrature_points,
                           output_data.JxW_values,
                           output_data.boundary_forms,
                           output_data.normal_vectors,
                           output_data.jacobians,
                           output_data.inverse_jacobians);
}


template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::
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

  const unsigned int n_q_points=quadrature.size();
  this->compute_fill_face (cell, face_no, subface_no,
                           n_q_points,
                           QProjector<dim>::DataSetDescriptor::
                           subface (face_no, subface_no,
                                    cell->face_orientation(face_no),
                                    cell->face_flip(face_no),
                                    cell->face_rotation(face_no),
                                    n_q_points,
                                    cell->subface_case(face_no)),
                           quadrature.get_weights(),
                           data,
                           output_data.quadrature_points,
                           output_data.JxW_values,
                           output_data.boundary_forms,
                           output_data.normal_vectors,
                           output_data.jacobians,
                           output_data.inverse_jacobians);
}



template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::transform (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  VectorSlice<std::vector<Tensor<1,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());

  transform_fields(input, output, mapping_data, mapping_type);
}



template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::transform (
  const VectorSlice<const std::vector<DerivativeForm<1, dim ,spacedim>  > >  input,
  VectorSlice<std::vector<Tensor<2,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());

  std::vector<DerivativeForm<1, spacedim,spacedim> > aux_output1(output.size());
  VectorSlice< std::vector<DerivativeForm<1, spacedim,spacedim> > >  aux_output( aux_output1);

  transform_differential_forms(input, aux_output, mapping_data, mapping_type);

  for (unsigned int i=0; i<output.size(); i++)
    output[i] = aux_output[i];
}


template<int dim, int spacedim, class VECTOR, class DH>
void MappingFEField<dim,spacedim,VECTOR,DH>::transform
(const VectorSlice<const std::vector<Tensor<2, dim> > >     input,
 VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
 const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
 const MappingType) const
{
  (void)input;
  (void)output;
  (void)mapping_data;
  AssertDimension (input.size(), output.size());

  AssertThrow(false, ExcNotImplemented());
}


template<int dim, int spacedim, class VECTOR, class DH>
template < int rank >
void MappingFEField<dim,spacedim,VECTOR,DH>::transform_fields(
  const VectorSlice<const std::vector<Tensor<rank,dim> > > input,
  VectorSlice<std::vector<Tensor<rank,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_contravariant:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        output[i] = apply_transformation(data.contravariant[i], input[i]);

      return;
    }

    case mapping_piola:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
      Assert (data.update_flags & update_volume_elements,
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
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        output[i] = apply_transformation(data.covariant[i], input[i]);

      return;
    }

    default:
      Assert(false, ExcNotImplemented());
    }
}


template<int dim, int spacedim, class VECTOR, class DH>
template < int rank >
void MappingFEField<dim,spacedim,VECTOR,DH>::transform_differential_forms(
  const VectorSlice<const std::vector<DerivativeForm<rank, dim,spacedim> > >    input,
  VectorSlice<std::vector<DerivativeForm<rank, spacedim,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{

  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_covariant:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        output[i] = apply_transformation(data.covariant[i], input[i]);

      return;
    }
    default:
      Assert(false, ExcNotImplemented());
    }

}



template<int dim, int spacedim, class VECTOR, class DH>
Point<spacedim>
MappingFEField<dim,spacedim,VECTOR,DH>::
transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<dim>                                 &p) const
{
//  Use the get_data function to create an InternalData with data vectors of
//  the right size and transformation shape values already computed at point
//  p.
  const Quadrature<dim> point_quadrature(p);
  std_cxx11::unique_ptr<InternalData>
  mdata (get_data(update_transformation_values, point_quadrature));

  update_internal_dofs(cell, *mdata);

  return this->transform_unit_to_real_cell_internal(*mdata);
}


template<int dim, int spacedim, class VECTOR, class DH>
Point<spacedim>
MappingFEField<dim,spacedim,VECTOR,DH>::
transform_unit_to_real_cell_internal (const InternalData &data) const
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



template<int dim, int spacedim, class VECTOR, class DH>
Point<dim>
MappingFEField<dim,spacedim,VECTOR,DH>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<spacedim>                            &p) const
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

  UpdateFlags update_flags = update_transformation_values|update_transformation_gradients;
  if (spacedim>dim)
    update_flags |= update_jacobian_grads;
  std_cxx11::unique_ptr<InternalData>
  mdata (get_data(update_flags,point_quadrature));

  update_internal_dofs(cell, *mdata);

  return this->transform_real_to_unit_cell_internal(cell, p, initial_p_unit, *mdata);

}


template<int dim, int spacedim, class VECTOR, class DH>
Point<dim>
MappingFEField<dim,spacedim,VECTOR,DH>::
transform_real_to_unit_cell_internal
(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
 const Point<spacedim>                            &p,
 const Point<dim>                                 &initial_p_unit,
 InternalData                                     &mdata) const
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
  Point<spacedim> p_real(transform_unit_to_real_cell_internal(mdata));
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
      Tensor<1, dim> delta;
      contract (delta, invert(df), static_cast<const Tensor<1,dim>&>(f));
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
          Point<spacedim> p_real_trial = transform_unit_to_real_cell_internal(mdata);
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


template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::compute_fill (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const unsigned int  n_q_points,
  const typename QProjector<dim>::DataSetDescriptor  data_set,
  const CellSimilarity::Similarity cell_similarity,
  const InternalData  &data,
  std::vector<Point<spacedim> > &quadrature_points) const
{
  const UpdateFlags update_flags(data.current_update_flags());
  update_internal_dofs(cell, data);

  // first compute quadrature points
  if (update_flags & update_quadrature_points)
    {
      AssertDimension (quadrature_points.size(), n_q_points);

      for (unsigned int point=0; point<n_q_points; ++point)
        {
          Point<spacedim> result;
          const double *shape = &data.shape(point+data_set,0);

          for (unsigned int k=0; k<data.n_shape_functions; ++k)
            {
              unsigned int comp_k = fe->system_to_component_index(k).first;
              if (fe_mask[comp_k])
                result[fe_to_real[comp_k]] += data.local_dof_values[k] * shape[k];
            }

          quadrature_points[point] = result;
        }
    }


  // then Jacobians
  if (update_flags & update_contravariant_transformation)
    {
      AssertDimension (data.contravariant.size(), n_q_points);

      // if the current cell is just a
      // translation of the previous one, no
      // need to recompute jacobians...
      if (cell_similarity != CellSimilarity::translation)
        {
          std::fill(data.contravariant.begin(), data.contravariant.end(),
                    DerivativeForm<1,dim,spacedim>());

          Assert (data.n_shape_functions > 0, ExcInternalError());

          for (unsigned int point=0; point<n_q_points; ++point)
            {
              const Tensor<1,dim> *data_derv =
                &data.derivative(point+data_set, 0);

              Tensor<1, dim> result[spacedim];

              for (unsigned int k=0; k<data.n_shape_functions; ++k)
                {
                  unsigned int comp_k = fe->system_to_component_index(k).first;
                  if (fe_mask[comp_k])
                    result[fe_to_real[comp_k]] += data.local_dof_values[k] * data_derv[k];
                }

              // write result into contravariant data. for
              // j=dim in the case dim<spacedim, there will
              // never be any nonzero data that arrives in
              // here, so it is ok anyway because it was
              // initialized to zero at the initialization
              for (unsigned int i=0; i<spacedim; ++i)
                {
                  data.contravariant[point][i] = result[i];
                }

            }
        }
    }


  if (update_flags & update_covariant_transformation)
    {
      AssertDimension (data.covariant.size(), n_q_points);
      if (cell_similarity != CellSimilarity::translation)
        for (unsigned int point=0; point<n_q_points; ++point)
          data.covariant[point] = (data.contravariant[point]).covariant_form();
    }

  if (update_flags & update_volume_elements)
    if (cell_similarity != CellSimilarity::translation)
      for (unsigned int point=0; point<n_q_points; ++point)
        data.volume_elements[point] = data.contravariant[point].determinant();

}


template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::compute_fill_face (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const unsigned int      face_no,
  const unsigned int      subface_no,
  const unsigned int      n_q_points,//npts
  const typename QProjector<dim>::DataSetDescriptor data_set,
  const std::vector<double>   &weights,
  const InternalData           &data,
  std::vector<Point<spacedim> >    &quadrature_points,
  std::vector<double>         &JxW_values,
  std::vector<Tensor<1,spacedim> > &boundary_forms,
  std::vector<Point<spacedim> > &normal_vectors,
  std::vector<DerivativeForm<1,dim,spacedim> > &/*jacobians*/,
  std::vector<DerivativeForm<1,spacedim,dim> > &/*inverse_jacobians*/) const
{
  compute_fill (cell, n_q_points, data_set, CellSimilarity::none,
                data, quadrature_points);


  const UpdateFlags update_flags(data.current_update_flags());

  if (update_flags & update_boundary_forms)
    {
      AssertDimension (boundary_forms.size(), n_q_points);
      if (update_flags & update_normal_vectors)
        AssertDimension (normal_vectors.size(), n_q_points);
      if (update_flags & update_JxW_values)
        AssertDimension (JxW_values.size(), n_q_points);

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

          transform (data.unit_tangentials[face_no+GeometryInfo<dim>::faces_per_cell*d],
                     data.aux[d],
                     data,
                     mapping_contravariant);
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
                boundary_forms[i][0] = (face_no == 0 ?
                                        -1 : +1);
                break;
              case 2:
                cross_product (boundary_forms[i], data.aux[0][i]);
                break;
              case 3:
                cross_product (boundary_forms[i], data.aux[0][i], data.aux[1][i]);
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
                  boundary_forms[point] = data.contravariant[point].transpose()[0];
                  boundary_forms[point] /=
                    (face_no == 0 ? -1. : +1.) * boundary_forms[point].norm();

                }

              if (dim==2)
                {
                  Tensor<1,spacedim> cell_normal;
                  const DerivativeForm<1,spacedim,dim> DX_t =
                    data.contravariant[point].transpose();
                  cross_product(cell_normal,DX_t[0],DX_t[1]);
                  cell_normal /= cell_normal.norm();

                  // then compute the face normal from the face tangent
                  // and the cell normal:
                  cross_product (boundary_forms[point],
                                 data.aux[0][point], cell_normal);

                }

            }
        }



      if (update_flags & (update_normal_vectors
                          | update_JxW_values))
        for (unsigned int i=0; i<boundary_forms.size(); ++i)
          {
            if (update_flags & update_JxW_values)
              {
                JxW_values[i] = boundary_forms[i].norm() * weights[i];

                if (subface_no != numbers::invalid_unsigned_int)
                  {
                    const double area_ratio=GeometryInfo<dim>::subface_ratio(
                                              cell->subface_case(face_no), subface_no);
                    JxW_values[i] *= area_ratio;
                  }
              }

            if (update_flags & update_normal_vectors)
              normal_vectors[i] = Point<spacedim>(boundary_forms[i] / boundary_forms[i].norm());
          }
    }

}


template<int dim, int spacedim, class VECTOR, class DH>
unsigned int
MappingFEField<dim,spacedim,VECTOR,DH>::get_degree() const
{
  return fe->degree;
}


template<int dim, int spacedim, class VECTOR, class DH>
ComponentMask
MappingFEField<dim,spacedim,VECTOR,DH>::get_component_mask() const
{
  return this->fe_mask;
}


template<int dim, int spacedim, class VECTOR, class DH>
Mapping<dim,spacedim> *
MappingFEField<dim,spacedim,VECTOR,DH>::clone () const
{
  return new MappingFEField<dim,spacedim,VECTOR,DH>(*this);
}


template<int dim, int spacedim, class VECTOR, class DH>
void
MappingFEField<dim,spacedim,VECTOR,DH>::update_internal_dofs (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const typename MappingFEField<dim, spacedim>::InternalData &data) const
{
  Assert(euler_dof_handler != 0, ExcMessage("euler_dof_handler is empty"));

  typename DH::cell_iterator dof_cell(*cell, euler_dof_handler);
  Assert (dof_cell->active() == true, ExcInactiveCell());

  dof_cell->get_dof_indices(data.local_dof_indices);

  for (unsigned int i=0; i<data.local_dof_values.size(); ++i)
    data.local_dof_values[i] = (*euler_vector)(data.local_dof_indices[i]);
}

// explicit instantiations
#include "mapping_fe_field.inst"


DEAL_II_NAMESPACE_CLOSE
