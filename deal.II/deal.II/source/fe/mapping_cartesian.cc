//-----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------

#include <base/tensor.h>
#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/mapping_cartesian.h>
#include <fe/fe_values.h>

#include <cmath>
#include <algorithm>




template <int dim>
const unsigned int MappingCartesian<dim>::invalid_face_number;



template<int dim>
MappingCartesian<dim>::InternalData::InternalData (const Quadrature<dim>& q)
		:
		quadrature_points (q.get_points ())
{}



template <int dim>
UpdateFlags
MappingCartesian<dim>::update_once (const UpdateFlags) const
{
  return update_default;
}



template <int dim>
UpdateFlags
MappingCartesian<dim>::update_each (const UpdateFlags in) const
{
  UpdateFlags out = in;
  if (out & update_boundary_forms)
    out |= update_normal_vectors;
  
  return out;
}



template <int dim>
typename Mapping<dim>::InternalDataBase *
MappingCartesian<dim>::get_data (const UpdateFlags      update_flags,
				 const Quadrature<dim> &q) const
{
  //  Assert (flags & update_normal_vectors == 0, ExcNotImplemented());

  InternalData* data = new InternalData (q);

  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  return data;
}



template <int dim>
typename Mapping<dim>::InternalDataBase *
MappingCartesian<dim>::get_face_data (const UpdateFlags update_flags,
				      const Quadrature<dim-1>& quadrature) const
{
  QProjector<dim> q (quadrature, false);
  InternalData* data = new InternalData (q);

  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  return data;
}



template <int dim>
typename Mapping<dim>::InternalDataBase *
MappingCartesian<dim>::get_subface_data (const UpdateFlags update_flags,
					 const Quadrature<dim-1> &quadrature) const
{
  QProjector<dim> q (quadrature, true);
  InternalData* data = new InternalData (q);

  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  return data;
}




template <int dim>
void
MappingCartesian<dim>::compute_fill (const typename DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int        face_no,
				     const unsigned int        sub_no,
				     InternalData             &data,
				     std::vector<Point<dim> > &quadrature_points,
				     std::vector<Point<dim> > &normal_vectors) const
{
  UpdateFlags update_flags(data.current_update_flags());

  const unsigned int npts = quadrature_points.size ();
  unsigned int offset = 0;
  
  if (face_no != invalid_face_number)
    {
				       // Add 1 on both sides of
				       // assertion to avoid compiler
				       // warning about testing
				       // unsigned int < 0 in 1d.
      Assert (face_no+1 < GeometryInfo<dim>::faces_per_cell+1,
	      ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));

      if (sub_no == invalid_face_number)
					 // called from FEFaceValues
	offset = face_no * quadrature_points.size();
      else
	{
					   // called from FESubfaceValue
	  Assert (sub_no+1 < GeometryInfo<dim>::subfaces_per_face+1,
		  ExcIndexRange (sub_no, 0, GeometryInfo<dim>::subfaces_per_face));
	  offset = (face_no * GeometryInfo<dim>::subfaces_per_face + sub_no)
		   * quadrature_points.size();
	}
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
  switch (dim)
    {
      case 1:
	data.length[0] = cell->vertex(1)(0) - start(0);
	break;
      case 2:
	data.length[0] = cell->vertex(1)(0) - start(0);
	data.length[1] = cell->vertex(3)(1) - start(1);
	break;
      case 3:
	data.length[0] = cell->vertex(1)(0) - start(0);
	data.length[1] = cell->vertex(4)(1) - start(1);
	data.length[2] = cell->vertex(3)(2) - start(2);
	break;
      default:
	Assert(false, ExcNotImplemented());
    }
  

				   // transform quadrature point. this
				   // is obtained simply by scaling
				   // unit coordinates with lengths in
				   // each direction
  if (update_flags & update_q_points)
    {
      Assert (quadrature_points.size() == npts,
	      ExcDimensionMismatch(quadrature_points.size(), npts));
      for (unsigned int i=0; i<npts; ++i)
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
				   // 1 or -1
  if (update_flags & update_normal_vectors)
    {
      Assert (normal_vectors.size() == npts,
	      ExcDimensionMismatch(normal_vectors.size(), npts));
      Point<dim> n;
      switch (100*dim+face_no)
	{
					   // 2D
	  case 200:
		n (1) = -1.;
		break;
	  case 201:
		n (0) = 1.;
		break;
	  case 202:
		n (1) = 1.;
		break;
	  case 203:
		n (0) = -1.;
		break;
					   // 3D
	  case 300:
		n (1) = -1.;
		break;
	  case 301:
		n (1) = 1.;
		break;
	  case 302:
		n (2) = -1.;
		break;
	  case 303:
		n (0) = 1.;
		break;
	  case 304:
		n (2) = 1.;
		break;
	  case 305:
		n (0) = -1.;
		break;
	  default:
		Assert (false, ExcInternalError());
	}
				       // furthermore, all normal
				       // vectors on a face are equal
      std::fill (normal_vectors.begin(), normal_vectors.end(), n);
    }
  
  data.first_cell = false;
}


template <int dim>
void
MappingCartesian<dim>::fill_fe_values (const typename DoFHandler<dim>::cell_iterator& cell,
				       const Quadrature<dim>& q, 
				       typename Mapping<dim>::InternalDataBase& mapping_data,
				       std::vector<Point<dim> >& quadrature_points,
				       std::vector<double>& JxW_values) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &data = dynamic_cast<InternalData&> (mapping_data);

  typename std::vector<Point<dim> > dummy;
  
  compute_fill (cell, invalid_face_number, invalid_face_number,
		data,
		quadrature_points,
		dummy);

				   // compute Jacobian
				   // determinant. all values are
				   // equal and are the product of the
				   // local lengths in each coordinate
				   // direction
//TODO:[GK] does it make sense to query the update flags in data if we did not set any in the get_*_data functions?
  if (data.current_update_flags() & update_JxW_values)
    {
      double J = data.length[0];
      for (unsigned int d=1;d<dim;++d)
	J *= data.length[d];
      for (unsigned int i=0; i<JxW_values.size();++i)
	JxW_values[i] = J * q.weight(i);
    }
}



template <int dim>
void
MappingCartesian<dim>::fill_fe_face_values (const typename DoFHandler<dim>::cell_iterator &cell,
					    const unsigned int            face_no,
					    const Quadrature<dim-1>      &q,
					    typename Mapping<dim>::InternalDataBase &mapping_data,
					    std::vector<Point<dim> >     &quadrature_points,
					    std::vector<double>          &JxW_values,
					    std::vector<Tensor<1,dim> >  &boundary_forms,
					    std::vector<Point<dim> >     &normal_vectors) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &data = dynamic_cast<InternalData&> (mapping_data);

  compute_fill (cell, face_no, invalid_face_number,
		data,
		quadrature_points,
		normal_vectors);

				   // first compute Jacobian
				   // determinant, which is simply the
				   // product of the local lengths
				   // since the jacobian is diagonal
  double J = 1.;
  for (unsigned int d=0;d<dim;++d)
    if (d != (normal_directions[face_no]/2))
      J *= data.length[d];
  
//TODO:[GK] does it make sense to query the update flags in data if we did not set any in the get_*_data functions?
  if (data.current_update_flags() & update_JxW_values)
    {
      for (unsigned int i=0; i<JxW_values.size();++i)
	JxW_values[i] = J * q.weight(i);
    }

//TODO:[GK] does it make sense to query the update flags in data if we did not set any in the get_*_data functions?
  if (data.current_update_flags() & update_boundary_forms)
    {
      for (unsigned int i=0; i<boundary_forms.size();++i)
	boundary_forms[i] = J * normal_vectors[i];
    }
}


template <int dim>
void
MappingCartesian<dim>::fill_fe_subface_values (const typename DoFHandler<dim>::cell_iterator &cell,
					const unsigned int       face_no,
					const unsigned int       sub_no,
					const Quadrature<dim-1> &q,
					typename Mapping<dim>::InternalDataBase &mapping_data,
					std::vector<Point<dim> >     &quadrature_points,
					std::vector<double>          &JxW_values,
					std::vector<Tensor<1,dim> >  &boundary_forms,
					std::vector<Point<dim> >     &normal_vectors) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &data = dynamic_cast<InternalData&> (mapping_data);

  compute_fill (cell, face_no, sub_no,
		data,
		quadrature_points,
		normal_vectors);

				   // first compute Jacobian
				   // determinant, which is simply the
				   // product of the local lengths
				   // since the jacobian is diagonal
  double J = 1.;
  for (unsigned int d=0;d<dim;++d)
    if (d != (normal_directions[face_no]/2))
      J *= data.length[d];
  
//TODO:[GK] does it make sense to query the update flags in data if we did not set any in the get_*_data functions?
  if (data.current_update_flags() & update_JxW_values)
    {
      for (unsigned int i=0; i<JxW_values.size();++i)
	JxW_values[i] = J * q.weight(i) / GeometryInfo<dim>::subfaces_per_face;
    }

//TODO:[GK] does it make sense to query the update flags in data if we did not set any in the get_*_data functions?
  if (data.current_update_flags() & update_boundary_forms)
    {
      for (unsigned int i=0; i<boundary_forms.size();++i)
	boundary_forms[i] = J * normal_vectors[i];
    }
}


#if (deal_II_dimension == 1)

template <>
void
MappingCartesian<1>::fill_fe_face_values (const DoFHandler<1>::cell_iterator &,
					  const unsigned,
					  const Quadrature<0>&,
					  Mapping<1>::InternalDataBase&,
					  std::vector<Point<1> >&,
					  std::vector<double>&,
					  std::vector<Tensor<1,1> >&,
					  std::vector<Point<1> >&) const
{
  Assert(false, ExcNotImplemented());
}


template <>
void
MappingCartesian<1>::fill_fe_subface_values (const DoFHandler<1>::cell_iterator &,
					     const unsigned,
					     const unsigned,
					     const Quadrature<0>&,
					     Mapping<1>::InternalDataBase&,
					     std::vector<Point<1> >&,
					     std::vector<double>&,
					     std::vector<Tensor<1,1> >&,
					     std::vector<Point<1> >&) const
{
  Assert(false, ExcNotImplemented());
}
#endif


template <int dim>
void
MappingCartesian<dim>::transform_covariant (std::vector<Tensor<1,dim> >       &dst,
					    const std::vector<Tensor<1,dim> > &src,
					    const typename Mapping<dim>::InternalDataBase &mapping_data,
					    const unsigned int src_offset) const
{
  covariant_transformation(dst, src, mapping_data, src_offset);
}


template <int dim>
void
MappingCartesian<dim>::transform_covariant (std::vector<Point<dim> >       &dst,
					    const std::vector<Point<dim> > &src,
					    const typename Mapping<dim>::InternalDataBase &mapping_data,
					    const unsigned int src_offset) const
{
  covariant_transformation(dst, src, mapping_data, src_offset);
}

template <int dim>
void
MappingCartesian<dim>::transform_contravariant (std::vector<Tensor<1,dim> >       &dst,
						const std::vector<Tensor<1,dim> > &src,
						const typename Mapping<dim>::InternalDataBase &mapping_data,
						const unsigned int src_offset) const
{
  contravariant_transformation(dst, src, mapping_data, src_offset);
}



template <int dim>
void
MappingCartesian<dim>::transform_contravariant (std::vector<Point<dim> >       &dst,
						const std::vector<Point<dim> > &src,
						const typename Mapping<dim>::InternalDataBase &mapping_data,
						const unsigned int src_offset) const
{
  contravariant_transformation(dst, src, mapping_data, src_offset);
}



template <int dim>
Point<dim> MappingCartesian<dim>::transform_unit_to_real_cell (
  const typename Triangulation<dim>::cell_iterator cell,
  const Point<dim> &p) const
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
	length[1] = cell->vertex(3)(1) - start(1);
	break;
      case 3:
	length[0] = cell->vertex(1)(0) - start(0);
	length[1] = cell->vertex(4)(1) - start(1);
	length[2] = cell->vertex(3)(2) - start(2);
	break;
      default:
	Assert(false, ExcNotImplemented());
    }

  Point<dim> p_real = cell->vertex(0);
  for (unsigned int d=0; d<dim; ++d)
    p_real(d) +=length[d]*p(d);

  return p_real;
}



template <int dim>
Point<dim> MappingCartesian<dim>::transform_real_to_unit_cell (
  const typename Triangulation<dim>::cell_iterator cell,
  const Point<dim> &p) const
{
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
	real(1) /= cell->vertex(3)(1) - start(1);
	break;
      case 3:
	real(0) /= cell->vertex(1)(0) - start(0);
	real(1) /= cell->vertex(4)(1) - start(1);
	real(2) /= cell->vertex(3)(2) - start(2);
	break;
      default:
	Assert(false, ExcNotImplemented());
    }
  return real;
}



template <int dim>
template <typename tensor_>
inline
void
MappingCartesian<dim>::contravariant_transformation (std::vector<tensor_>       &dst,
						     const std::vector<tensor_> &src,
						     const typename Mapping<dim>::InternalDataBase &mapping_data,
						     const unsigned int src_offset) const
{
  Assert(tensor_::dimension==dim && tensor_::rank==1, ExcInvalidData());

				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  const InternalData &data = dynamic_cast<const InternalData&> (mapping_data);

  Assert (data.update_flags & update_contravariant_transformation,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());

  typename std::vector<tensor_>::const_iterator vec    = src.begin() + src_offset;
  typename std::vector<tensor_>::iterator       result = dst.begin();
  typename std::vector<tensor_>::const_iterator end    = dst.end();
  
				   // simply scale by Jacobian
				   // (which is diagonal here)
  while (result!=end)
    {
      for (unsigned int d=0; d<dim; ++d)
	(*result)[d] = (*vec)[d]*data.length[d];
      ++vec;
      ++result;
    }
}



template <int dim>
template <typename tensor_>
inline
void
MappingCartesian<dim>::covariant_transformation (std::vector<tensor_>       &dst,
						 const std::vector<tensor_> &src,
						 const typename Mapping<dim>::InternalDataBase &mapping_data,
						 const unsigned int src_offset) const
{
  Assert(tensor_::dimension==dim && tensor_::rank==1, ExcInvalidData());

				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  const InternalData &data = dynamic_cast<const InternalData&> (mapping_data);

  Assert (data.update_flags & update_covariant_transformation,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());
  
  typename std::vector<tensor_>::const_iterator vec    = src.begin() + src_offset;
  typename std::vector<tensor_>::iterator       result = dst.begin();
  typename std::vector<tensor_>::const_iterator end    = dst.end();

				   // simply scale by inverse Jacobian
				   // (which is diagonal here)
  while (result!=end)
    {
      for (unsigned int d=0;d<dim;++d)
	(*result)[d] = (*vec)[d]/data.length[d];
      vec++;
      result++;
    }
}



//----------------------------------------------------------------------//
// explicit instantiations

template class MappingCartesian<deal_II_dimension>;

template void MappingCartesian<deal_II_dimension>::contravariant_transformation (
  std::vector<Tensor<1,deal_II_dimension> >       &dst,
  const std::vector<Tensor<1,deal_II_dimension> > &src,
  const Mapping<deal_II_dimension>::InternalDataBase& internal,
  const unsigned int src_offset) const;

template void MappingCartesian<deal_II_dimension>::contravariant_transformation (
  std::vector<Point<deal_II_dimension> >       &dst,
  const std::vector<Point<deal_II_dimension> > &src,
  const Mapping<deal_II_dimension>::InternalDataBase& internal,
  const unsigned int src_offset) const;

template void MappingCartesian<deal_II_dimension>::covariant_transformation (
  std::vector<Tensor<1,deal_II_dimension> >       &dst,
  const std::vector<Tensor<1,deal_II_dimension> > &src,
  const Mapping<deal_II_dimension>::InternalDataBase& internal,
  const unsigned int src_offset) const;

template void MappingCartesian<deal_II_dimension>::covariant_transformation (
  std::vector<Point<deal_II_dimension> >       &dst,
  const std::vector<Point<deal_II_dimension> > &src,
  const Mapping<deal_II_dimension>::InternalDataBase& internal,
  const unsigned int src_offset) const;

