// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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


#include <deal.II/base/tensor.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/fe_values.h>

#include <cmath>
#include <algorithm>


DEAL_II_NAMESPACE_OPEN


template<int dim, int spacedim>
const unsigned int MappingCartesian<dim,spacedim>::invalid_face_number;



template<int dim, int spacedim>
MappingCartesian<dim, spacedim>::InternalData::InternalData (const Quadrature<dim> &q)
  :
  quadrature_points (q.get_points ())
{}



template<int dim, int spacedim>
std::size_t
MappingCartesian<dim, spacedim>::InternalData::memory_consumption () const
{
  return (Mapping<dim, spacedim>::InternalDataBase::memory_consumption() +
          MemoryConsumption::memory_consumption (length) +
          MemoryConsumption::memory_consumption (quadrature_points));
}




template<int dim, int spacedim>
UpdateFlags
MappingCartesian<dim, spacedim>::update_once (const UpdateFlags) const
{
  return update_default;
}



template<int dim, int spacedim>
UpdateFlags
MappingCartesian<dim, spacedim>::update_each (const UpdateFlags in) const
{
  UpdateFlags out = in;
  if (out & update_boundary_forms)
    out |= update_normal_vectors;

  return out;
}



template<int dim, int spacedim>
typename Mapping<dim, spacedim>::InternalDataBase *
MappingCartesian<dim, spacedim>::get_data (const UpdateFlags      update_flags,
                                           const Quadrature<dim> &q) const
{
  InternalData *data = new InternalData (q);

  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  return data;
}



template<int dim, int spacedim>
typename Mapping<dim, spacedim>::InternalDataBase *
MappingCartesian<dim, spacedim>::get_face_data (const UpdateFlags update_flags,
                                                const Quadrature<dim-1>& quadrature) const
{
  InternalData *data
    = new InternalData (QProjector<dim>::project_to_all_faces(quadrature));

  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  return data;
}



template<int dim, int spacedim>
typename Mapping<dim, spacedim>::InternalDataBase *
MappingCartesian<dim, spacedim>::get_subface_data (const UpdateFlags update_flags,
                                                   const Quadrature<dim-1> &quadrature) const
{
  InternalData *data
    = new InternalData (QProjector<dim>::project_to_all_subfaces(quadrature));

  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  return data;
}




template<int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::compute_fill (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                               const unsigned int        face_no,
                                               const unsigned int        sub_no,
                                               const CellSimilarity::Similarity cell_similarity,
                                               InternalData             &data,
                                               std::vector<Point<dim> > &quadrature_points,
                                               std::vector<Point<dim> > &normal_vectors) const
{
  const UpdateFlags update_flags(data.current_update_flags());

  // some more sanity checks
  if (face_no != invalid_face_number)
    {
      // Add 1 on both sides of
      // assertion to avoid compiler
      // warning about testing
      // unsigned int < 0 in 1d.
      Assert (face_no+1 < GeometryInfo<dim>::faces_per_cell+1,
              ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));

      // We would like to check for
      // sub_no < cell->face(face_no)->n_children(),
      // but unfortunately the current
      // function is also called for
      // faces without children (see
      // tests/fe/mapping.cc). Therefore,
      // we must use following workaround
      // of two separate assertions
      Assert ((sub_no == invalid_face_number) ||
              cell->face(face_no)->has_children() ||
              (sub_no+1 < GeometryInfo<dim>::max_children_per_face+1),
              ExcIndexRange (sub_no, 0,
                             GeometryInfo<dim>::max_children_per_face));
      Assert ((sub_no == invalid_face_number) ||
              !cell->face(face_no)->has_children() ||
              (sub_no < cell->face(face_no)->n_children()),
              ExcIndexRange (sub_no, 0, cell->face(face_no)->n_children()));
    }
  else
    // invalid face number, so
    // subface should be invalid as
    // well
    Assert (sub_no == invalid_face_number, ExcInternalError());

  // let @p{start} be the origin of a
  // local coordinate system. it is
  // chosen as the (lower) left
  // vertex
  const Point<dim> start = cell->vertex(0);

  // Compute start point and sizes
  // along axes.  Strange vertex
  // numbering makes this complicated
  // again.
  if (cell_similarity != CellSimilarity::translation)
    {
      switch (dim)
        {
        case 1:
          data.length[0] = cell->vertex(1)(0) - start(0);
          break;
        case 2:
          data.length[0] = cell->vertex(1)(0) - start(0);
          data.length[1] = cell->vertex(2)(1) - start(1);
          break;
        case 3:
          data.length[0] = cell->vertex(1)(0) - start(0);
          data.length[1] = cell->vertex(2)(1) - start(1);
          data.length[2] = cell->vertex(4)(2) - start(2);
          break;
        default:
          Assert(false, ExcNotImplemented());
        }
    }


  // transform quadrature point. this
  // is obtained simply by scaling
  // unit coordinates with lengths in
  // each direction
  if (update_flags & update_quadrature_points)
    {
      const typename QProjector<dim>::DataSetDescriptor offset
        = (face_no == invalid_face_number
           ?
           QProjector<dim>::DataSetDescriptor::cell()
           :
           (sub_no == invalid_face_number
            ?
            // called from FEFaceValues
            QProjector<dim>::DataSetDescriptor::face (face_no,
                                                      cell->face_orientation(face_no),
                                                      cell->face_flip(face_no),
                                                      cell->face_rotation(face_no),
                                                      quadrature_points.size())
            :
            // called from FESubfaceValues
            QProjector<dim>::DataSetDescriptor::subface (face_no, sub_no,
                                                         cell->face_orientation(face_no),
                                                         cell->face_flip(face_no),
                                                         cell->face_rotation(face_no),
                                                         quadrature_points.size(),
                                                         cell->subface_case(face_no))
           ));

      for (unsigned int i=0; i<quadrature_points.size(); ++i)
        {
          quadrature_points[i] = start;
          for (unsigned int d=0; d<dim; ++d)
            quadrature_points[i](d) += data.length[d] *
                                       data.quadrature_points[i+offset](d);
        }
    }


  // compute normal vectors. since
  // cells are aligned to coordinate
  // axes, they are simply vectors
  // with exactly one entry equal to
  // 1 or -1. Furthermore, all
  // normals on a face have the same
  // value
  if (update_flags & update_normal_vectors)
    {
      Assert (face_no < GeometryInfo<dim>::faces_per_cell,
              ExcInternalError());

      switch (dim)
        {
        case 1:
        {
          static const Point<dim>
          normals[GeometryInfo<1>::faces_per_cell]
            = { Point<dim>(-1.),
                Point<dim>( 1.)
              };
          std::fill (normal_vectors.begin(),
                     normal_vectors.end(),
                     normals[face_no]);
          break;
        }

        case 2:
        {
          static const Point<dim>
          normals[GeometryInfo<2>::faces_per_cell]
            = { Point<dim>(-1, 0),
                Point<dim>( 1, 0),
                Point<dim>( 0,-1),
                Point<dim>( 0, 1)
              };
          std::fill (normal_vectors.begin(),
                     normal_vectors.end(),
                     normals[face_no]);
          break;
        }

        case 3:
        {
          static const Point<dim>
          normals[GeometryInfo<3>::faces_per_cell]
            = { Point<dim>(-1, 0, 0),
                Point<dim>( 1, 0, 0),
                Point<dim>( 0,-1, 0),
                Point<dim>( 0, 1, 0),
                Point<dim>( 0, 0,-1),
                Point<dim>( 0, 0, 1)
              };
          std::fill (normal_vectors.begin(),
                     normal_vectors.end(),
                     normals[face_no]);
          break;
        }

        default:
          Assert (false, ExcNotImplemented());
        }
    }
}



template<int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::
fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                const Quadrature<dim> &q,
                typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                std::vector<Point<spacedim> > &quadrature_points,
                std::vector<double> &JxW_values,
                std::vector< DerivativeForm<1,dim,spacedim> > &jacobians,
                std::vector<DerivativeForm<2,dim,spacedim> >      &jacobian_grads,
                std::vector<DerivativeForm<1,spacedim,dim> >      &inverse_jacobians,
                std::vector<Point<spacedim> > &,
                CellSimilarity::Similarity &cell_similarity) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<InternalData *> (&mapping_data) != 0, ExcInternalError());
  InternalData &data = static_cast<InternalData &> (mapping_data);

  std::vector<Point<dim> > dummy;

  compute_fill (cell, invalid_face_number, invalid_face_number, cell_similarity,
                data,
                quadrature_points,
                dummy);

  // compute Jacobian
  // determinant. all values are
  // equal and are the product of the
  // local lengths in each coordinate
  // direction
  if (data.current_update_flags() & (update_JxW_values | update_volume_elements))
    if (cell_similarity != CellSimilarity::translation)
      {
        double J = data.length[0];
        for (unsigned int d=1; d<dim; ++d)
          J *= data.length[d];
        data.volume_element = J;
        if (data.current_update_flags() & update_JxW_values)
          for (unsigned int i=0; i<JxW_values.size(); ++i)
            JxW_values[i] = J * q.weight(i);
      }
  // "compute" Jacobian at the quadrature
  // points, which are all the same
  if (data.current_update_flags() & update_jacobians)
    if (cell_similarity != CellSimilarity::translation)
      for (unsigned int i=0; i<jacobians.size(); ++i)
        {
          jacobians[i] =  DerivativeForm<1,dim,spacedim>();
          for (unsigned int j=0; j<dim; ++j)
            jacobians[i][j][j]=data.length[j];
        }
  // "compute" the derivative of the Jacobian
  // at the quadrature points, which are all
  // zero of course
  if (data.current_update_flags() & update_jacobian_grads)
    if (cell_similarity != CellSimilarity::translation)
      for (unsigned int i=0; i<jacobian_grads.size(); ++i)
        jacobian_grads[i]=DerivativeForm<2,dim,spacedim>();
  // "compute" inverse Jacobian at the
  // quadrature points, which are all
  // the same
  if (data.current_update_flags() & update_inverse_jacobians)
    if (cell_similarity != CellSimilarity::translation)
      for (unsigned int i=0; i<inverse_jacobians.size(); ++i)
        {
          inverse_jacobians[i]=Tensor<2,dim>();
          for (unsigned int j=0; j<dim; ++j)
            inverse_jacobians[j][j]=1./data.length[j];
        }
}



// template <>
// void
// MappingCartesian<1,1>::fill_fe_face_values (
//   const typename Triangulation<1,1>::cell_iterator &,
//   const unsigned int,
//   const Quadrature<0>&,
//   typename Mapping<1,1>::InternalDataBase&,
//   std::vector<Point<1> >&,
//   std::vector<double>&,
//   std::vector<Tensor<1,1> >&,
//   std::vector<Point<1> >&) const
// {
//   Assert(false, ExcNotImplemented());
// }


// template <>
// void
// MappingCartesian<1,2>::fill_fe_face_values (
//   const typename Triangulation<1,2>::cell_iterator &,
//   const unsigned int,
//   const Quadrature<0>&,
//   typename Mapping<1,2>::InternalDataBase&,
//   std::vector<Point<1> >&,
//   std::vector<double>&,
//   std::vector<Tensor<1,1> >&,
//   std::vector<Point<2> >&) const
// {
//   Assert(false, ExcNotImplemented());
// }



// template <>
// void
// MappingCartesian<1,1>::fill_fe_subface_values (
//   const typename Triangulation<1,1>::cell_iterator &,
//   const unsigned int,
//   const unsigned int,
//   const Quadrature<0>&,
//   typename Mapping<1,1>::InternalDataBase&,
//   std::vector<Point<1> >&,
//   std::vector<double>&,
//   std::vector<Tensor<1,1> >&,
//   std::vector<Point<1> >&) const
// {
//   Assert(false, ExcNotImplemented());
// }



// template <>
// void
// MappingCartesian<1,2>::fill_fe_subface_values (
//   const typename Triangulation<1,2>::cell_iterator &,
//   const unsigned int,
//   const unsigned int,
//   const Quadrature<0>&,
//   typename Mapping<1,2>::InternalDataBase&,
//   std::vector<Point<1> >&,
//   std::vector<double>&,
//   std::vector<Tensor<1,1> >&,
//   std::vector<Point<2> >&) const
// {
//   Assert(false, ExcNotImplemented());
// }



// Implementation for dim != 1

template<int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::fill_fe_face_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const unsigned int            face_no,
  const Quadrature<dim-1>      &q,
  typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  std::vector<Point<dim> >     &quadrature_points,
  std::vector<double>          &JxW_values,
  std::vector<Tensor<1,dim> >  &boundary_forms,
  std::vector<Point<spacedim> >     &normal_vectors) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<InternalData *> (&mapping_data) != 0,
          ExcInternalError());
  InternalData &data = static_cast<InternalData &> (mapping_data);

  compute_fill (cell, face_no, invalid_face_number,
                CellSimilarity::none,
                data,
                quadrature_points,
                normal_vectors);

  // first compute Jacobian
  // determinant, which is simply the
  // product of the local lengths
  // since the jacobian is diagonal
  double J = 1.;
  for (unsigned int d=0; d<dim; ++d)
    if (d != GeometryInfo<dim>::unit_normal_direction[face_no])
      J *= data.length[d];

  if (data.current_update_flags() & update_JxW_values)
    for (unsigned int i=0; i<JxW_values.size(); ++i)
      JxW_values[i] = J * q.weight(i);

  if (data.current_update_flags() & update_boundary_forms)
    for (unsigned int i=0; i<boundary_forms.size(); ++i)
      boundary_forms[i] = J * normal_vectors[i];

  if (data.current_update_flags() & update_volume_elements)
    {
      J = data.length[0];
      for (unsigned int d=1; d<dim; ++d)
        J *= data.length[d];
      data.volume_element = J;
    }
}



template<int dim, int spacedim>
void
MappingCartesian<dim, spacedim>::fill_fe_subface_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const unsigned int       face_no,
  const unsigned int       sub_no,
  const Quadrature<dim-1> &q,
  typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
  std::vector<Point<dim> >     &quadrature_points,
  std::vector<double>          &JxW_values,
  std::vector<Tensor<1,dim> >  &boundary_forms,
  std::vector<Point<spacedim> >     &normal_vectors) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert (dynamic_cast<InternalData *> (&mapping_data) != 0, ExcInternalError());
  InternalData &data = static_cast<InternalData &> (mapping_data);

  compute_fill (cell, face_no, sub_no, CellSimilarity::none,
                data,
                quadrature_points,
                normal_vectors);

  // first compute Jacobian
  // determinant, which is simply the
  // product of the local lengths
  // since the jacobian is diagonal
  double J = 1.;
  for (unsigned int d=0; d<dim; ++d)
    if (d != GeometryInfo<dim>::unit_normal_direction[face_no])
      J *= data.length[d];

  if (data.current_update_flags() & update_JxW_values)
    {
      // Here,
      // cell->face(face_no)->n_children()
      // would be the right choice,
      // but unfortunately the
      // current function is also
      // called for faces without
      // children (see
      // tests/fe/mapping.cc). Add
      // following switch to avoid
      // diffs in tests/fe/mapping.OK
      const unsigned int n_subfaces=
        cell->face(face_no)->has_children() ?
        cell->face(face_no)->n_children() :
        GeometryInfo<dim>::max_children_per_face;
      for (unsigned int i=0; i<JxW_values.size(); ++i)
        JxW_values[i] = J * q.weight(i) / n_subfaces;
    }

  if (data.current_update_flags() & update_boundary_forms)
    for (unsigned int i=0; i<boundary_forms.size(); ++i)
      boundary_forms[i] = J * normal_vectors[i];

  if (data.current_update_flags() & update_volume_elements)
    {
      J = data.length[0];
      for (unsigned int d=1; d<dim; ++d)
        J *= data.length[d];
      data.volume_element = J;
    }
}




template<int dim, int spacedim>
void
MappingCartesian<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  VectorSlice<std::vector<Tensor<1,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (mapping_data);

  switch (mapping_type)
    {
    case mapping_covariant:
    {
      Assert (data.update_flags & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d=0; d<dim; ++d)
          output[i][d] = input[i][d]/data.length[d];
      return;
    }

    case mapping_contravariant:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d=0; d<dim; ++d)
          output[i][d] = input[i][d]*data.length[d];
      return;
    }
    case mapping_piola:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
      Assert (data.update_flags & update_volume_elements,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d=0; d<dim; ++d)
          output[i][d] = input[i][d] * data.length[d] / data.volume_element;
      return;
    }
    default:
      Assert(false, ExcNotImplemented());
    }
}



template<int dim, int spacedim>
void
MappingCartesian<dim,spacedim>::transform (
  const VectorSlice<const std::vector<DerivativeForm<1, dim,spacedim> > > input,
  VectorSlice<std::vector<Tensor<2,spacedim> > > output,
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
      Assert (data.update_flags & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int d2=0; d2<dim; ++d2)
            output[i][d1][d2] = input[i][d1][d2] / data.length[d2];
      return;
    }

    case mapping_contravariant:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int d2=0; d2<dim; ++d2)
            output[i][d1][d2] = input[i][d1][d2] * data.length[d2];
      return;
    }

    case mapping_covariant_gradient:
    {
      Assert (data.update_flags & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int d2=0; d2<dim; ++d2)
            output[i][d1][d2] = input[i][d1][d2] / data.length[d2] / data.length[d1];
      return;
    }

    case mapping_contravariant_gradient:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int d2=0; d2<dim; ++d2)
            output[i][d1][d2] = input[i][d1][d2] * data.length[d2] / data.length[d1];
      return;
    }

    case mapping_piola:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
      Assert (data.update_flags & update_volume_elements,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int d2=0; d2<dim; ++d2)
            output[i][d1][d2] = input[i][d1][d2] * data.length[d2]
                                / data.length[d1] / data.volume_element;
      return;
    }

    default:
      Assert(false, ExcNotImplemented());
    }
}




template<int dim, int spacedim>
void
MappingCartesian<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<2, dim> > >    input,
  VectorSlice<std::vector<Tensor<2, spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{

  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_piola_gradient:
    {
      Assert (data.update_flags & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
      Assert (data.update_flags & update_volume_elements,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));

      for (unsigned int i=0; i<output.size(); ++i)
        for (unsigned int d1=0; d1<dim; ++d1)
          for (unsigned int d2=0; d2<dim; ++d2)
            output[i][d1][d2] = input[i][d1][d2] * data.length[d2]
                                / data.length[d1] / data.volume_element;
      return;
    }

    default:
      Assert(false, ExcNotImplemented());
    }

}


template <int dim, int spacedim>
Point<spacedim>
MappingCartesian<dim, spacedim>::transform_unit_to_real_cell (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const Point<dim>                                 &p) const
{
  Tensor<1,dim> length;
  const Point<dim> start = cell->vertex(0);
  switch (dim)
    {
    case 1:
      length[0] = cell->vertex(1)(0) - start(0);
      break;
    case 2:
      length[0] = cell->vertex(1)(0) - start(0);
      length[1] = cell->vertex(2)(1) - start(1);
      break;
    case 3:
      length[0] = cell->vertex(1)(0) - start(0);
      length[1] = cell->vertex(2)(1) - start(1);
      length[2] = cell->vertex(4)(2) - start(2);
      break;
    default:
      Assert(false, ExcNotImplemented());
    }

  Point<dim> p_real = cell->vertex(0);
  for (unsigned int d=0; d<dim; ++d)
    p_real(d) += length[d]*p(d);

  return p_real;
}



template<int dim, int spacedim>
Point<dim>
MappingCartesian<dim, spacedim>::transform_real_to_unit_cell (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const Point<spacedim>                            &p) const
{

  if (dim != spacedim)
    Assert(false, ExcNotImplemented());
  const Point<dim> &start = cell->vertex(0);
  Point<dim> real = p;
  real -= start;

  switch (dim)
    {
    case 1:
      real(0) /= cell->vertex(1)(0) - start(0);
      break;
    case 2:
      real(0) /= cell->vertex(1)(0) - start(0);
      real(1) /= cell->vertex(2)(1) - start(1);
      break;
    case 3:
      real(0) /= cell->vertex(1)(0) - start(0);
      real(1) /= cell->vertex(2)(1) - start(1);
      real(2) /= cell->vertex(4)(2) - start(2);
      break;
    default:
      Assert(false, ExcNotImplemented());
    }
  return real;
}


template<int dim, int spacedim>
Mapping<dim, spacedim> *
MappingCartesian<dim, spacedim>::clone () const
{
  return new MappingCartesian<dim, spacedim>(*this);
}


//---------------------------------------------------------------------------
// explicit instantiations
#include "mapping_cartesian.inst"


DEAL_II_NAMESPACE_CLOSE
