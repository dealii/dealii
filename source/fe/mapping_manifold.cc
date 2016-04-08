// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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
#include <deal.II/grid/manifold.h>
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
  return (Mapping<dim,spacedim>::InternalDataBase::memory_consumption() +
          MemoryConsumption::memory_consumption (vertices) +
          MemoryConsumption::memory_consumption (cell) +
          MemoryConsumption::memory_consumption (quad) +
          MemoryConsumption::memory_consumption (cell_manifold_quadrature_weights) +
          MemoryConsumption::memory_consumption (vertex_weights) +
          MemoryConsumption::memory_consumption (unit_tangentials) +
          MemoryConsumption::memory_consumption (covariant) +
          MemoryConsumption::memory_consumption (contravariant) +
          MemoryConsumption::memory_consumption (aux) +
          MemoryConsumption::memory_consumption (volume_elements) +
          MemoryConsumption::memory_consumption (manifold) );
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

  // Store the quadrature
  this->quad = q;

  // Resize the weights
  this->vertex_weights.resize(GeometryInfo<dim>::vertices_per_cell);

  // see if we need the (transformation) shape function values
  // and/or gradients and resize the necessary arrays
  if (this->update_each & (update_quadrature_points | update_contravariant_transformation) )
    compute_manifold_quadrature_weights(q);

  if (this->update_each & update_covariant_transformation)
    covariant.resize(n_original_q_points);

  if (this->update_each & update_contravariant_transformation)
    contravariant.resize(n_original_q_points);

  if (this->update_each & update_volume_elements)
    volume_elements.resize(n_original_q_points);
}


template <int dim, int spacedim>
void
MappingManifold<dim,spacedim>::InternalData::
initialize_face (const UpdateFlags      update_flags,
                 const Quadrature<dim> &q,
                 const unsigned int     n_original_q_points)
{
  initialize (update_flags, q, n_original_q_points);

  if (dim > 1)
    {
      if (this->update_each & update_boundary_forms)
        {
          aux.resize (dim-1, std::vector<Tensor<1,spacedim> > (n_original_q_points));

          // Compute tangentials to the unit cell.
          for (unsigned int i=0; i<unit_tangentials.size(); ++i)
            unit_tangentials[i].resize (n_original_q_points);
          switch (dim)
            {
            case 2:
            {
              // ensure a counterclockwise
              // orientation of tangentials
              static const int tangential_orientation[4]= {-1,1,1,-1};
              for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
                {
                  Tensor<1,dim> tang;
                  tang[1-i/2]=tangential_orientation[i];
                  std::fill (unit_tangentials[i].begin(),
                             unit_tangentials[i].end(),
                             tang);
                }
              break;
            }
            case 3:
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
                  std::fill (unit_tangentials[i].begin(),
                             unit_tangentials[i].end(), tang1);
                  std::fill (unit_tangentials[GeometryInfo<dim>::faces_per_cell+i].begin(),
                             unit_tangentials[GeometryInfo<dim>::faces_per_cell+i].end(),
                             tang2);
                }
              break;
            }
            default:
              Assert(false,ExcNotImplemented());
            }
        }
    }
}



template<int dim, int spacedim>
MappingManifold<dim,spacedim>::MappingManifold ()
{}



template<int dim, int spacedim>
MappingManifold<dim,spacedim>::MappingManifold (const MappingManifold<dim,spacedim> &)
{}



template<int dim, int spacedim>
Mapping<dim,spacedim> *
MappingManifold<dim,spacedim>::clone () const
{
  return new MappingManifold<dim,spacedim>(*this);
}



template<int dim, int spacedim>
Point<dim>
MappingManifold<dim,spacedim>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &,
                             const Point<spacedim> &) const
{
  Assert(false, ExcNotImplemented());
  return Point<dim>();
}



template<int dim, int spacedim>
Point<spacedim>
MappingManifold<dim,spacedim>::
transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<dim> &p) const
{
  std::vector<Point<spacedim> > vertices;
  std::vector<double> weights;
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
    {
      vertices.push_back(cell->vertex(v));
      weights.push_back(GeometryInfo<dim>::d_linear_shape_function(p,v));
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

      // The contravariant transformation used in the Piola
      // transformation, which requires the determinant of the Jacobi
      // matrix of the transformation.  Because we have no way of
      // knowing here whether the finite elements wants to use the
      // contravariant of the Piola transforms, we add the JxW values
      // to the list of flags to be updated for each cell.
      if (out & update_contravariant_transformation)
        out |= update_JxW_values;

      if (out & update_normal_vectors)
        out |= update_JxW_values;
    }

  // Now throw an exception if we stumble upon something that was not
  // implemented yet
  Assert(!(out & (
             update_jacobian_grads |
             update_jacobian_pushed_forward_grads |
             update_jacobian_pushed_forward_2nd_derivatives |
             update_jacobian_pushed_forward_3rd_derivatives
           )), ExcNotImplemented());

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
     * Some specialization for face Manifolds. In one dimension, there
     * are no Manifolds associated to faces. The mapping argument is
     * only used to help the compiler infer dim and spacedim.
     */
    template<int spacedim>
    const dealii::Manifold<1, spacedim> &
    get_face_manifold(const MappingManifold<1,spacedim> &,
                      const typename dealii::Triangulation<1,spacedim>::cell_iterator &cell,
                      const unsigned int &)
    {
      return cell->get_manifold();
    }



    /**
     * Some specialization for face Manifolds. The mapping argument is
     * only used to help the compiler infer dim and spacedim.
     */
    template<int dim, int spacedim>
    const dealii::Manifold<dim,spacedim> &
    get_face_manifold(const MappingManifold<dim,spacedim> &,
                      const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                      const unsigned int face_no)
    {
      return cell->face(face_no)->get_manifold();
    }



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

      AssertDimension(data.vertices.size(), GeometryInfo<dim>::vertices_per_cell);

      if (update_flags & update_quadrature_points)
        {
          for (unsigned int point=0; point<quadrature_points.size(); ++point)
            {
              quadrature_points[point] = data.manifold->
                                         get_new_point(Quadrature<spacedim>(data.vertices,
                                                       data.cell_manifold_quadrature_weights[point+data_set]));
            }
        }
    }



    /**
     * Update the co- and contravariant matrices as well as their determinant, for the cell
     * described stored in the data object, but only if the update_flags of the @p data
     * argument indicate so.
     */
    template <int dim, int spacedim>
    void
    maybe_update_Jacobians (const typename dealii::QProjector<dim>::DataSetDescriptor          data_set,
                            const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_contravariant_transformation)
        {
          const unsigned int n_q_points = data.contravariant.size();

          std::fill(data.contravariant.begin(), data.contravariant.end(),
                    DerivativeForm<1,dim,spacedim>());

          AssertDimension(GeometryInfo<dim>::vertices_per_cell,
                          data.vertices.size());
          for (unsigned int point=0; point<n_q_points; ++point)
            {
              // Start by figuring out how to compute the direction in
              // the reference space:
              const Point<dim> &p = data.quad.point(point+data_set);

              // And get its image on the manifold:
              const Point<spacedim> P = data.manifold->
                                        get_new_point(Quadrature<spacedim>(data.vertices,
                                                      data.cell_manifold_quadrature_weights[point+data_set]));

              // To compute the Jacobian, we choose dim points aligned
              // with with the dim reference axes, which are still in the
              // given cell, and ask for the tangent vector in these
              // directions. Choosing the points is somewhat arbitrary,
              // so we try to be smart and we pick points which are
              // on the opposite quadrant w.r.t. the evaluation
              // point.
              for (unsigned int i=0; i<dim; ++i)
                {
                  const Point<dim> ei = Point<dim>::unit_vector(i);
                  const double pi = p[i];
                  Assert(pi >=0 && pi <= 1.0,
                         ExcInternalError("Was expecting a quadrature point "
                                          "inside the unit reference element."));
                  const Point<dim> np(pi > .5 ? p-pi *ei : p+(1-pi)*ei);

                  // In the length L, we store also the direction sign,
                  // which is positive, if the coordinate is < .5,
                  double L = pi > .5 ? -pi: 1-pi;

                  // Get the weights to compute the np point in real space
                  for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_cell; ++j)
                    data.vertex_weights[j] = GeometryInfo<dim>::d_linear_shape_function(np, j);

                  const Point<spacedim> NP=
                    data.manifold->get_new_point(Quadrature<spacedim>(data.vertices,
                                                                      data.vertex_weights));

                  Tensor<1,spacedim> T = data.manifold->get_tangent_vector(P, NP);

                  for (unsigned int d=0; d<spacedim; ++d)
                    data.contravariant[point][d][i] = T[d]/L;
                }
            }

          if (update_flags & update_covariant_transformation)
            {
              const unsigned int n_q_points = data.contravariant.size();
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  data.covariant[point] = (data.contravariant[point]).covariant_form();
                }
            }

          if (update_flags & update_volume_elements)
            {
              const unsigned int n_q_points = data.contravariant.size();
              for (unsigned int point=0; point<n_q_points; ++point)
                data.volume_elements[point] = data.contravariant[point].determinant();
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

  data.store_vertices(cell);
  data.manifold = &(cell->get_manifold());

  internal::maybe_compute_q_points<dim,spacedim> (QProjector<dim>::DataSetDescriptor::cell (),
                                                  data,
                                                  output_data.quadrature_points);

  internal::maybe_update_Jacobians<dim,spacedim> (QProjector<dim>::DataSetDescriptor::cell (),
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
                = std::sqrt(determinant(G)) * weights[point];

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
                      Assert(spacedim == dim+1,
                             ExcMessage("There is no (unique) cell normal for "
                                        + Utilities::int_to_string(dim) +
                                        "-dimensional cells in "
                                        + Utilities::int_to_string(spacedim) +
                                        "-dimensional space. This only works if the "
                                        "space dimension is one greater than the "
                                        "dimensionality of the mesh cells."));

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
                        const double area_ratio=GeometryInfo<dim>::subface_ratio(cell->subface_case(face_no),
                                                                                 subface_no);
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
    do_fill_fe_face_values (const dealii::MappingManifold<dim,spacedim>                       &mapping,
                            const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                            const unsigned int                                                 face_no,
                            const unsigned int                                                 subface_no,
                            const typename QProjector<dim>::DataSetDescriptor                  data_set,
                            const Quadrature<dim-1>                                           &quadrature,
                            const typename dealii::MappingManifold<dim,spacedim>::InternalData      &data,
                            internal::FEValues::MappingRelatedData<dim,spacedim>              &output_data)
    {
      data.store_vertices(cell);

      data.manifold = &get_face_manifold(mapping, cell, face_no);

      maybe_compute_q_points<dim,spacedim> (data_set,
                                            data,
                                            output_data.quadrature_points);
      maybe_update_Jacobians<dim,spacedim> (data_set,
                                            data);

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

//--------------------------- Explicit instantiations -----------------------
#include "mapping_manifold.inst"


DEAL_II_NAMESPACE_CLOSE

