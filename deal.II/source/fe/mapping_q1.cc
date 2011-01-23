//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/tensor.h>
#include <base/quadrature.h>
#include <base/qprojector.h>
#include <base/memory_consumption.h>
#include <lac/full_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <cmath>
#include <algorithm>
#include <memory>


DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
const unsigned int MappingQ1<dim,spacedim>::n_shape_functions;



template<int dim, int spacedim>
MappingQ1<dim,spacedim>::InternalData::InternalData (const unsigned int n_shape_functions)
		:
		is_mapping_q1_data(true),
		n_shape_functions (n_shape_functions)
{}



template<int dim, int spacedim>
std::size_t
MappingQ1<dim,spacedim>::InternalData::memory_consumption () const
{
  return (Mapping<dim,spacedim>::InternalDataBase::memory_consumption() +
	  MemoryConsumption::memory_consumption (shape_values) +
	  MemoryConsumption::memory_consumption (shape_derivatives) +
	  MemoryConsumption::memory_consumption (covariant) +
	  MemoryConsumption::memory_consumption (contravariant) +
	  MemoryConsumption::memory_consumption (unit_tangentials) +
	  MemoryConsumption::memory_consumption (aux) +
	  MemoryConsumption::memory_consumption (mapping_support_points) +
	  MemoryConsumption::memory_consumption (cell_of_current_support_points) +
	  MemoryConsumption::memory_consumption (is_mapping_q1_data) +
	  MemoryConsumption::memory_consumption (n_shape_functions));
}



template<int dim, int spacedim>
MappingQ1<dim,spacedim>::MappingQ1 ()
{}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_shapes (const std::vector<Point<dim> > &unit_points,
					 InternalData &data) const
{
				   // choose either the function implemented
				   // in this class, or whatever a virtual
				   // function call resolves to
  if (data.is_mapping_q1_data)
    MappingQ1<dim,spacedim>::compute_shapes_virtual(unit_points, data);
  else
    compute_shapes_virtual(unit_points, data);
}


namespace internal
{
  namespace MappingQ1
  {
    template <int spacedim>
    void
    compute_shapes_virtual (const unsigned int            n_shape_functions,
			    const std::vector<Point<1> > &unit_points,
			    typename dealii::MappingQ1<1,spacedim>::InternalData& data)
    {
      const unsigned int n_points=unit_points.size();
      for (unsigned int k = 0 ; k < n_points ; ++k)
	{
	  double x = unit_points[k](0);

	  if (data.shape_values.size()!=0)
	    {
	      Assert(data.shape_values.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.shape(k,0) = 1.-x;
	      data.shape(k,1) = x;
	    }
	  if (data.shape_derivatives.size()!=0)
	    {
	      Assert(data.shape_derivatives.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.derivative(k,0)[0] = -1.;
	      data.derivative(k,1)[0] = 1.;
	    }
	  if (data.shape_second_derivatives.size()!=0)
	    {
					       // the following may or may not
					       // work if dim != spacedim
	      Assert (spacedim == 1, ExcNotImplemented());

	      Assert(data.shape_second_derivatives.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.second_derivative(k,0)[0][0] = 0;
	      data.second_derivative(k,1)[0][0] = 0;
	    }
	}
    }


    template <int spacedim>
    void
    compute_shapes_virtual (const unsigned int            n_shape_functions,
			    const std::vector<Point<2> > &unit_points,
			    typename dealii::MappingQ1<2,spacedim>::InternalData& data)
    {
      const unsigned int n_points=unit_points.size();
      for (unsigned int k = 0 ; k < n_points ; ++k)
	{
	  double x = unit_points[k](0);
	  double y = unit_points[k](1);

	  if (data.shape_values.size()!=0)
	    {
	      Assert(data.shape_values.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.shape(k,0) = (1.-x)*(1.-y);
	      data.shape(k,1) = x*(1.-y);
	      data.shape(k,2) = (1.-x)*y;
	      data.shape(k,3) = x*y;
	    }
	  if (data.shape_derivatives.size()!=0)
	    {
	      Assert(data.shape_derivatives.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.derivative(k,0)[0] = (y-1.);
	      data.derivative(k,1)[0] = (1.-y);
	      data.derivative(k,2)[0] = -y;
	      data.derivative(k,3)[0] = y;
	      data.derivative(k,0)[1] = (x-1.);
	      data.derivative(k,1)[1] = -x;
	      data.derivative(k,2)[1] = (1.-x);
	      data.derivative(k,3)[1] = x;
	    }
	  if (data.shape_second_derivatives.size()!=0)
	    {
					       // the following may or may not
					       // work if dim != spacedim
	      Assert (spacedim == 2, ExcNotImplemented());

	      Assert(data.shape_second_derivatives.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.second_derivative(k,0)[0][0] = 0;
	      data.second_derivative(k,1)[0][0] = 0;
	      data.second_derivative(k,2)[0][0] = 0;
	      data.second_derivative(k,3)[0][0] = 0;
	      data.second_derivative(k,0)[0][1] = 1.;
	      data.second_derivative(k,1)[0][1] = -1.;
	      data.second_derivative(k,2)[0][1] = -1.;
	      data.second_derivative(k,3)[0][1] = 1.;
	      data.second_derivative(k,0)[1][0] = 1.;
	      data.second_derivative(k,1)[1][0] = -1.;
	      data.second_derivative(k,2)[1][0] = -1.;
	      data.second_derivative(k,3)[1][0] = 1.;
	      data.second_derivative(k,0)[1][1] = 0;
	      data.second_derivative(k,1)[1][1] = 0;
	      data.second_derivative(k,2)[1][1] = 0;
	      data.second_derivative(k,3)[1][1] = 0;
	    }
	}
    }



    template <int spacedim>
    void
    compute_shapes_virtual (const unsigned int            n_shape_functions,
			    const std::vector<Point<3> > &unit_points,
			    typename dealii::MappingQ1<3,spacedim>::InternalData& data)
    {
      const unsigned int n_points=unit_points.size();
      for (unsigned int k = 0 ; k < n_points ; ++k)
	{
	  double x = unit_points[k](0);
	  double y = unit_points[k](1);
	  double z = unit_points[k](2);

	  if (data.shape_values.size()!=0)
	    {
	      Assert(data.shape_values.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.shape(k,0) = (1.-x)*(1.-y)*(1.-z);
	      data.shape(k,1) = x*(1.-y)*(1.-z);
	      data.shape(k,2) = (1.-x)*y*(1.-z);
	      data.shape(k,3) = x*y*(1.-z);
	      data.shape(k,4) = (1.-x)*(1.-y)*z;
	      data.shape(k,5) = x*(1.-y)*z;
	      data.shape(k,6) = (1.-x)*y*z;
	      data.shape(k,7) = x*y*z;
	    }
	  if (data.shape_derivatives.size()!=0)
	    {
	      Assert(data.shape_derivatives.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.derivative(k,0)[0] = (y-1.)*(1.-z);
	      data.derivative(k,1)[0] = (1.-y)*(1.-z);
	      data.derivative(k,2)[0] = -y*(1.-z);
	      data.derivative(k,3)[0] = y*(1.-z);
	      data.derivative(k,4)[0] = (y-1.)*z;
	      data.derivative(k,5)[0] = (1.-y)*z;
	      data.derivative(k,6)[0] = -y*z;
	      data.derivative(k,7)[0] = y*z;
	      data.derivative(k,0)[1] = (x-1.)*(1.-z);
	      data.derivative(k,1)[1] = -x*(1.-z);
	      data.derivative(k,2)[1] = (1.-x)*(1.-z);
	      data.derivative(k,3)[1] = x*(1.-z);
	      data.derivative(k,4)[1] = (x-1.)*z;
	      data.derivative(k,5)[1] = -x*z;
	      data.derivative(k,6)[1] = (1.-x)*z;
	      data.derivative(k,7)[1] = x*z;
	      data.derivative(k,0)[2] = (x-1)*(1.-y);
	      data.derivative(k,1)[2] = x*(y-1.);
	      data.derivative(k,2)[2] = (x-1.)*y;
	      data.derivative(k,3)[2] = -x*y;
	      data.derivative(k,4)[2] = (1.-x)*(1.-y);
	      data.derivative(k,5)[2] = x*(1.-y);
	      data.derivative(k,6)[2] = (1.-x)*y;
	      data.derivative(k,7)[2] = x*y;
	    }
	  if (data.shape_second_derivatives.size()!=0)
	    {
					       // the following may or may not
					       // work if dim != spacedim
	      Assert (spacedim == 3, ExcNotImplemented());

	      Assert(data.shape_second_derivatives.size()==n_shape_functions*n_points,
		     ExcInternalError());
	      data.second_derivative(k,0)[0][0] = 0;
	      data.second_derivative(k,1)[0][0] = 0;
	      data.second_derivative(k,2)[0][0] = 0;
	      data.second_derivative(k,3)[0][0] = 0;
	      data.second_derivative(k,4)[0][0] = 0;
	      data.second_derivative(k,5)[0][0] = 0;
	      data.second_derivative(k,6)[0][0] = 0;
	      data.second_derivative(k,7)[0][0] = 0;
	      data.second_derivative(k,0)[1][1] = 0;
	      data.second_derivative(k,1)[1][1] = 0;
	      data.second_derivative(k,2)[1][1] = 0;
	      data.second_derivative(k,3)[1][1] = 0;
	      data.second_derivative(k,4)[1][1] = 0;
	      data.second_derivative(k,5)[1][1] = 0;
	      data.second_derivative(k,6)[1][1] = 0;
	      data.second_derivative(k,7)[1][1] = 0;
	      data.second_derivative(k,0)[2][2] = 0;
	      data.second_derivative(k,1)[2][2] = 0;
	      data.second_derivative(k,2)[2][2] = 0;
	      data.second_derivative(k,3)[2][2] = 0;
	      data.second_derivative(k,4)[2][2] = 0;
	      data.second_derivative(k,5)[2][2] = 0;
	      data.second_derivative(k,6)[2][2] = 0;
	      data.second_derivative(k,7)[2][2] = 0;

	      data.second_derivative(k,0)[0][1] = (1.-z);
	      data.second_derivative(k,1)[0][1] = -(1.-z);
	      data.second_derivative(k,2)[0][1] = -(1.-z);
	      data.second_derivative(k,3)[0][1] = (1.-z);
	      data.second_derivative(k,4)[0][1] = z;
	      data.second_derivative(k,5)[0][1] = -z;
	      data.second_derivative(k,6)[0][1] = -z;
	      data.second_derivative(k,7)[0][1] = z;
	      data.second_derivative(k,0)[1][0] = (1.-z);
	      data.second_derivative(k,1)[1][0] = -(1.-z);
	      data.second_derivative(k,2)[1][0] = -(1.-z);
	      data.second_derivative(k,3)[1][0] = (1.-z);
	      data.second_derivative(k,4)[1][0] = z;
	      data.second_derivative(k,5)[1][0] = -z;
	      data.second_derivative(k,6)[1][0] = -z;
	      data.second_derivative(k,7)[1][0] = z;

	      data.second_derivative(k,0)[0][2] = (1.-y);
	      data.second_derivative(k,1)[0][2] = -(1.-y);
	      data.second_derivative(k,2)[0][2] = y;
	      data.second_derivative(k,3)[0][2] = -y;
	      data.second_derivative(k,4)[0][2] = -(1.-y);
	      data.second_derivative(k,5)[0][2] = (1.-y);
	      data.second_derivative(k,6)[0][2] = -y;
	      data.second_derivative(k,7)[0][2] = y;
	      data.second_derivative(k,0)[2][0] = (1.-y);
	      data.second_derivative(k,1)[2][0] = -(1.-y);
	      data.second_derivative(k,2)[2][0] = y;
	      data.second_derivative(k,3)[2][0] = -y;
	      data.second_derivative(k,4)[2][0] = -(1.-y);
	      data.second_derivative(k,5)[2][0] = (1.-y);
	      data.second_derivative(k,6)[2][0] = -y;
	      data.second_derivative(k,7)[2][0] = y;

	      data.second_derivative(k,0)[1][2] = (1.-x);
	      data.second_derivative(k,1)[1][2] = x;
	      data.second_derivative(k,2)[1][2] = -(1.-x);
	      data.second_derivative(k,3)[1][2] = -x;
	      data.second_derivative(k,4)[1][2] = -(1.-x);
	      data.second_derivative(k,5)[1][2] = -x;
	      data.second_derivative(k,6)[1][2] = (1.-x);
	      data.second_derivative(k,7)[1][2] = x;
	      data.second_derivative(k,0)[2][1] = (1.-x);
	      data.second_derivative(k,1)[2][1] = x;
	      data.second_derivative(k,2)[2][1] = -(1.-x);
	      data.second_derivative(k,3)[2][1] = -x;
	      data.second_derivative(k,4)[2][1] = -(1.-x);
	      data.second_derivative(k,5)[2][1] = -x;
	      data.second_derivative(k,6)[2][1] = (1.-x);
	      data.second_derivative(k,7)[2][1] = x;
	    }
	}
    }
  }
}


template<int dim, int spacedim>
void
MappingQ1<dim, spacedim>::
compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
			InternalData& data) const
{
  internal::MappingQ1::
    compute_shapes_virtual<spacedim> (n_shape_functions,
				      unit_points, data);
}



template<int dim, int spacedim>
UpdateFlags
MappingQ1<dim,spacedim>::update_once (const UpdateFlags in) const
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



template<int dim, int spacedim>
UpdateFlags
MappingQ1<dim,spacedim>::update_each (const UpdateFlags in) const
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
				   // iterate. since there are 4
				   // if-clauses in the loop, it will
				   // take at most 3 iterations to
				   // converge. do them:
  for (unsigned int i=0; i<4; ++i)
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
				       // Therefore these values have to
				       // updated for each cell.
      if (out & update_contravariant_transformation)
	out |= update_JxW_values;

      if (out & update_normal_vectors)
	out |= update_JxW_values;
    }

  return out;
}


template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_data (const UpdateFlags      update_flags,
				       const Quadrature<dim> &q,
				       const unsigned int     n_original_q_points,
				       InternalData          &data) const
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

  compute_shapes (q.get_points(), data);
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ1<dim,spacedim>::get_data (const UpdateFlags update_flags,
				   const Quadrature<dim>& q) const
{
  InternalData* data = new InternalData(n_shape_functions);
  compute_data (update_flags, q, q.size(), *data);
  return data;
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_face_data (const UpdateFlags update_flags,
					    const Quadrature<dim>& q,
					    const unsigned int n_original_q_points,
					    InternalData& data) const
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
	      static const int tangential_orientation[4]={-1,1,1,-1};
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



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ1<dim,spacedim>::get_face_data (const UpdateFlags        update_flags,
					const Quadrature<dim-1> &quadrature) const
{
  InternalData* data = new InternalData(n_shape_functions);
  compute_face_data (update_flags,
		     QProjector<dim>::project_to_all_faces(quadrature),
		     quadrature.size(),
		     *data);

  return data;
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ1<dim,spacedim>::get_subface_data (const UpdateFlags update_flags,
					   const Quadrature<dim-1>& quadrature) const
{
  InternalData* data = new InternalData(n_shape_functions);
  compute_face_data (update_flags,
		     QProjector<dim>::project_to_all_subfaces(quadrature),
		     quadrature.size(),
		     *data);

  return data;
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_fill (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
			               const unsigned int  n_q_points,
       		                       const DataSetDescriptor  data_set,
				       const CellSimilarity::Similarity cell_similarity,
	      	                       InternalData  &data,
		                       std::vector<Point<spacedim> > &quadrature_points) const
{
  const UpdateFlags update_flags(data.current_update_flags());

				   // if necessary, recompute the
				   // support points of the
				   // transformation of this cell
				   // (note that we need to first
				   // check the triangulation pointer,
				   // since otherwise the second test
				   // might trigger an exception if
				   // the triangulations are not the
				   // same)
  if ((data.mapping_support_points.size() == 0)
      ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation())
      ||
      (cell != data.cell_of_current_support_points))
    {
      compute_mapping_support_points(cell, data.mapping_support_points);
      data.cell_of_current_support_points = cell;
    }

                                   // first compute quadrature points
  if (update_flags & update_quadrature_points)
    {
      Assert (quadrature_points.size() == n_q_points,
	      ExcDimensionMismatch(quadrature_points.size(), n_q_points));
      std::fill(quadrature_points.begin(), quadrature_points.end(),
		Point<spacedim>());

      for (unsigned int point=0; point<n_q_points; ++point)
	for (unsigned int k=0; k<data.n_shape_functions; ++k)
	  quadrature_points[point]
	    += data.shape(point+data_set,k) * data.mapping_support_points[k];
    }

                                   // then Jacobians
  if (update_flags & update_contravariant_transformation)
    {
      Assert (data.contravariant.size() == n_q_points,
	      ExcDimensionMismatch(data.contravariant.size(), n_q_points));

				   // if the current cell is just a
				   // translation of the previous one, no
				   // need to recompute jacobians...
      if (cell_similarity != CellSimilarity::translation)
	{
	  std::fill(data.contravariant.begin(), data.contravariant.end(),
		    Tensor<2,spacedim>());
	  for (unsigned int point=0; point<n_q_points; ++point)
	    for (unsigned int k=0; k<data.n_shape_functions; ++k)
	      {
				   // some compilers seem to have problems
				   // to detected that the two innermost
				   // loops just use the data of the same
				   // tensors, so get a reference to them by
				   // hand
		const Tensor<1,dim> &data_derv = data.derivative(point+data_set, k);
		const Tensor<1,spacedim> &supp_pts = data.mapping_support_points[k];

		for (unsigned int i=0; i<spacedim; ++i)
		  for (unsigned int j=0; j<dim; ++j)
		    data.contravariant[point][i][j] += data_derv[j] * supp_pts[i];
	      }
	}
    }

  if (update_flags & update_covariant_transformation)
    {
      Assert (data.covariant.size() == n_q_points,
	      ExcDimensionMismatch(data.covariant.size(), n_q_points));
      if (cell_similarity != CellSimilarity::translation)
	{
	  if (dim == spacedim)
				    // invert contravariant for
			            // covariant transformation
				    // matrices
	    for (unsigned int point=0; point<n_q_points; ++point)
	      data.covariant[point] = invert(data.contravariant[point]);

	  else if (dim == spacedim - 1)
	    {
				    // CODIMENSION 1
				    // calculate left-inversion of the
				    // contravariant matrices to obtain
				    // covariant ones (auxiliary
				    // rectangular fullmatrices are used)
	      FullMatrix<double> contravariant_matrix(spacedim,dim);
	      FullMatrix<double> covariant_matrix(dim,spacedim);

	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  contravariant_matrix.
		    copy_from(data.contravariant[point],0,spacedim-1,0,dim-1);
		  covariant_matrix.
		    left_invert(contravariant_matrix);
		  covariant_matrix.
		    copy_to(data.covariant[point],0,dim-1,0,spacedim-1);
		}
	    }
	}
    }

  if (update_flags & update_volume_elements)
    if (cell_similarity != CellSimilarity::translation)
      for (unsigned int point=0; point<n_q_points; ++point)
	data.volume_elements[point] = determinant(data.contravariant[point]);
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  std::vector<Point<spacedim> > &a) const
{
  a.resize(GeometryInfo<dim>::vertices_per_cell);

  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a[i] = cell->vertex(i);
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::fill_fe_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const Quadrature<dim>                                     &q,
  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
  std::vector<Point<spacedim> >                             &quadrature_points,
  std::vector<double>                                       &JxW_values,
  std::vector<Tensor<2,spacedim> >                          &jacobians,
  std::vector<Tensor<3,spacedim> >                          &jacobian_grads,
  std::vector<Tensor<2,spacedim> >                          &inverse_jacobians,
  std::vector<Point<spacedim> >                             &normal_vectors,
  CellSimilarity::Similarity                           &cell_similarity) const
{
				   // ensure that the following cast
				   // is really correct:
  Assert (dynamic_cast<InternalData *>(&mapping_data) != 0,
	  ExcInternalError());
  InternalData &data = static_cast<InternalData&>(mapping_data);

  const unsigned int n_q_points=q.size();

  compute_fill (cell, n_q_points, DataSetDescriptor::cell (), cell_similarity,
                data, quadrature_points);


  const UpdateFlags update_flags(data.current_update_flags());
  const std::vector<double> &weights=q.get_weights();

				   // Multiply quadrature weights by
				   // Jacobian determinants or by norm
				   // of the cross product of the
				   // columns of the Jacobian and
				   // store the cell normal vectors in
				   // the case <2,3>
  if (update_flags & (update_normal_vectors
		      | update_JxW_values))
    {
      Assert (JxW_values.size() == n_q_points,
	      ExcDimensionMismatch(JxW_values.size(), n_q_points));

      Assert( !(update_flags & update_normal_vectors ) ||
	      (normal_vectors.size() == n_q_points),
	      ExcDimensionMismatch(normal_vectors.size(), n_q_points));

      if (cell_similarity != CellSimilarity::translation)
	for (unsigned int point=0; point<n_q_points; ++point)
	  {
	    if (dim==spacedim)
					       // if dim==spacedim,
					       // then there is no
					       // cell normal to
					       // compute. since this
					       // is for FEValues (and
					       // not FEFaceValues),
					       // there are also no
					       // face normals to
					       // compute
	      JxW_values[point]
		= determinant(data.contravariant[point])*weights[point];
	    else
	      {
		if (cell_similarity == CellSimilarity::inverted_translation)
		  {
						     // we only need to flip the normal
		    if (update_flags & update_normal_vectors)
		      normal_vectors[point] *= -1.;
		  }
		else
		  {
						     // temporarily
						     // transpose the
						     // matrix so that
						     // we can refer
						     // to its columns
						     // using a single
						     // index
		    data.contravariant[point] = transpose(data.contravariant[point]);

						     // compute the normal
						     // vector to this cell
						     // and put it into the
						     // last row of
						     // data.contravariant
		    if ( (dim==1) && (spacedim==2) )
		      cross_product(data.contravariant[point][1],
				    -data.contravariant[point][0]);
		    else if ( (dim==2) && (spacedim==3) )
		      cross_product(data.contravariant[point][2],
				    data.contravariant[point][0],
				    data.contravariant[point][1]);
		    else
		      Assert (false, ExcNotImplemented());

						     // det(J) is now the norm
						     // of the cross product
						     // of all the mapped unit
						     // tangential vectors
		    JxW_values[point]
		      = data.contravariant[point][spacedim-1].norm()*weights[point];

						     // in order to compute
						     // the normal vector,
						     // normalize the cross
						     // product of mapped unit
						     // vectors
		    data.contravariant[point][spacedim-1]
		      /= data.contravariant[point][spacedim-1].norm();

		    if (update_flags & update_normal_vectors)
		      {
			normal_vectors[point]
			  = data.contravariant[point][spacedim-1];

			if (cell->direction_flag() == false)
			  normal_vectors[point] *= -1.;
		      }

						     // un-transpose
						     // the matrix
						     // again
		    data.contravariant[point] = transpose(data.contravariant[point]);
		  }
	      }
	  }
    }
				   // copy values from InternalData to
				   // vector given by reference
  if (update_flags & update_jacobians)
    {
      Assert (jacobians.size() == n_q_points,
	      ExcDimensionMismatch(jacobians.size(), n_q_points));
      if (cell_similarity != CellSimilarity::translation)
	for (unsigned int point=0; point<n_q_points; ++point)
	  jacobians[point] = data.contravariant[point];
    }
				   // calculate values of the
				   // derivatives of the Jacobians. do
				   // it here, since we only do it for
				   // cells, not faces.
  if (update_flags & update_jacobian_grads)
    {
      Assert (jacobian_grads.size() == n_q_points,
	      ExcDimensionMismatch(jacobian_grads.size(), n_q_points));

      if (cell_similarity != CellSimilarity::translation)
	{
	  std::fill(jacobian_grads.begin(),
		    jacobian_grads.end(),
		    Tensor<3,spacedim>());

	  for (unsigned int point=0; point<n_q_points; ++point)
	    for (unsigned int k=0; k<data.n_shape_functions; ++k)
	      for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=0; j<dim; ++j)
		  for (unsigned int l=0; l<dim; ++l)
		    jacobian_grads[point][i][j][l]
		      += (data.second_derivative(point+DataSetDescriptor::cell (), k)[j][l]
			  *
			  data.mapping_support_points[k][i]);
	}
    }
				   // copy values from InternalData to vector
				   // given by reference
  if (update_flags & update_inverse_jacobians)
    {
      Assert (inverse_jacobians.size() == n_q_points,
	      ExcDimensionMismatch(inverse_jacobians.size(), n_q_points));
      if (cell_similarity != CellSimilarity::translation)
	for (unsigned int point=0; point<n_q_points; ++point)
	  inverse_jacobians[point] = transpose(data.covariant[point]);
    }
}



template <>
void
MappingQ1<1,1>::fill_fe_face_values (
  const Triangulation<1,1>::cell_iterator &,
  const unsigned,
  const Quadrature<0>&,
  Mapping<1,1>::InternalDataBase&,
  std::vector<Point<1> >&,
  std::vector<double>&,
  std::vector<Tensor<1,1> >&,
  std::vector<Point<1> >&) const
{
  Assert(false, ExcNotImplemented());
}



template <>
void
MappingQ1<1,2>::fill_fe_face_values (
  const Triangulation<1,2>::cell_iterator &,
  const unsigned,
  const Quadrature<0>&,
  Mapping<1,2>::InternalDataBase&,
  std::vector<Point<2> >&,
  std::vector<double>&,
  std::vector<Tensor<1,2> >&,
  std::vector<Point<2> >&) const
{
  Assert(false, ExcNotImplemented());
}



template <>
void
MappingQ1<1,1>::fill_fe_subface_values (
  const Triangulation<1,1>::cell_iterator &,
  const unsigned,
  const unsigned,
  const Quadrature<0>&,
  Mapping<1,1>::InternalDataBase&,
  std::vector<Point<1> >&,
  std::vector<double>&,
  std::vector<Tensor<1,1> >&,
  std::vector<Point<1> >&) const
{
  Assert(false, ExcNotImplemented());
}



template <>
void
MappingQ1<1,2>::fill_fe_subface_values (
  const Triangulation<1,2>::cell_iterator &,
  const unsigned,
  const unsigned,
  const Quadrature<0>&,
  Mapping<1,2>::InternalDataBase&,
  std::vector<Point<2> >&,
  std::vector<double>&,
  std::vector<Tensor<1,2> >&,
  std::vector<Point<2> >&) const
{
  Assert(false, ExcNotImplemented());
}



namespace internal
{
  namespace
  {
    template <int spacedim>
    void
    compute_fill_face (const dealii::MappingQ1<1,spacedim> &,
		       const typename dealii::Triangulation<1,spacedim>::cell_iterator &,
		       const unsigned int,
		       const unsigned int,
		       const unsigned int,
		       const std::vector<double> &,
		       typename dealii::MappingQ1<1,spacedim>::InternalData &,
		       std::vector<double> &,
		       std::vector<Tensor<1,spacedim> > &,
		       std::vector<Point<spacedim> > &)
    {
      Assert(false, ExcNotImplemented());
    }



    template <int dim, int spacedim>
    void
    compute_fill_face (const dealii::MappingQ1<dim,spacedim> &mapping,
		       const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
		       const unsigned int               face_no,
		       const unsigned int               subface_no,
		       const unsigned int               n_q_points,
		       const std::vector<double>        &weights,
		       typename dealii::MappingQ1<dim,spacedim>::InternalData &data,
		       std::vector<double>              &JxW_values,
		       std::vector<Tensor<1,spacedim> > &boundary_forms,
		       std::vector<Point<spacedim> >    &normal_vectors)
    {
      const UpdateFlags update_flags(data.current_update_flags());

      if (update_flags & update_boundary_forms)
      	{
      	  Assert (boundary_forms.size()==n_q_points,
      		  ExcDimensionMismatch(boundary_forms.size(), n_q_points));
      	  if (update_flags & update_normal_vectors)
      	    Assert (normal_vectors.size()==n_q_points,
      		    ExcDimensionMismatch(normal_vectors.size(), n_q_points));
      	  if (update_flags & update_JxW_values)
      	    Assert (JxW_values.size() == n_q_points,
      		    ExcDimensionMismatch(JxW_values.size(), n_q_points));

					   // map the unit tangentials to the
					   // real cell
	  for (unsigned int d=0; d<dim-1; ++d)
	    {
	      Assert (face_no+GeometryInfo<dim>::faces_per_cell*d <
		      data.unit_tangentials.size(),
		      ExcInternalError());
	      Assert (data.aux[d].size() <=
		      data.unit_tangentials[face_no+GeometryInfo<dim>::faces_per_cell*d].size(),
		      ExcInternalError());

	      mapping.transform (data.unit_tangentials[face_no+GeometryInfo<dim>::faces_per_cell*d],
				 data.aux[d],
				 data,
				 mapping_contravariant);
	    }

					   // if dim==spacedim, we can use the
					   // unit tangentials to compute the
					   // boundary form by simply taking
					   // the cross product
	  if (dim == spacedim)
	    {
	      for (unsigned int i=0; i<n_q_points; ++i)
		if (dim == 2)
      	  	  cross_product (boundary_forms[i], data.aux[0][i]);
		else if (dim == 3)
      	  	  cross_product (boundary_forms[i], data.aux[0][i], data.aux[1][i]);
		else
		  Assert(false, ExcNotImplemented());
      	    }
	  else
	    {
					       // in the codim-one case, the
					       // boundary form results from
					       // the cross product of all the
					       // face tangential vectors and
					       // the cell normal vector
					       //
					       // to compute the cell normal,
					       // use the same method used in
					       // fill_fe_values for cells
					       // above
	      Assert (data.contravariant.size() == n_q_points,
		      ExcInternalError());
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  Tensor<1,spacedim> cell_normal;

						   // temporarily
						   // transpose the
						   // matrix so that
						   // we can refer
						   // to its columns
						   // using a single
						   // index
		  data.contravariant[point] = transpose(data.contravariant[point]);

						   // compute the normal
						   // vector to this cell
						   // and put it into the
						   // last row of
						   // data.contravariant
		  if ( (dim==1) && (spacedim==2) )
		    cross_product(cell_normal,
				  -data.contravariant[point][0]);
		  else if ( (dim==2) && (spacedim==3) )
		    cross_product(cell_normal,
				  data.contravariant[point][0],
				  data.contravariant[point][1]);
		  else
		    Assert (false, ExcNotImplemented());

						   // then compute the face
						   // normal from the face
						   // tangent and the cell
						   // normal:
		  if ( (dim==1) && (spacedim==2) )
		    {
						       // need to think how to
						       // do this. the code
						       // should be something
						       // like the below, but
						       // I can't test it
						       // right now since we
						       // can't use
						       // FEFaceValues in 1d
		      Assert (false, ExcNotImplemented());
		      cross_product (boundary_forms[point],
				     (face_no == 0 ? 1 : -1) * cell_normal);
		    }
		  else if ( (dim==2) && (spacedim==3) )
		    cross_product (boundary_forms[point],
				   data.aux[0][point], cell_normal);
		  else
		    Assert (false, ExcNotImplemented());

						   // un-transpose
						   // the matrix
						   // again
		  data.contravariant[point] = transpose(data.contravariant[point]);
		}
	    }


      	  if (update_flags & (update_normal_vectors
      	  		      | update_JxW_values))
      	    for (unsigned int i=0;i<boundary_forms.size();++i)
      	      {
      	  	if (update_flags & update_JxW_values)
      	  	  {
      	  	    JxW_values[i] = boundary_forms[i].norm() * weights[i];
      	  	    if (subface_no!=deal_II_numbers::invalid_unsigned_int)
      	  	      {
      	  		const double area_ratio=GeometryInfo<dim>::subface_ratio(
      	  		  cell->subface_case(face_no), subface_no);
      	  		JxW_values[i] *= area_ratio;
      	  	      }
      	  	  }

      	  	if (update_flags & update_normal_vectors)
      	  	  normal_vectors[i] = boundary_forms[i] / boundary_forms[i].norm();
      	      }
	}
    }
  }
}


template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_fill_face (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const unsigned int               face_no,
  const unsigned int               subface_no,
  const unsigned int               n_q_points,
  const DataSetDescriptor          data_set,
  const std::vector<double>        &weights,
  InternalData                     &data,
  std::vector<Point<spacedim> >         &quadrature_points,
  std::vector<double>              &JxW_values,
  std::vector<Tensor<1,spacedim> > &boundary_forms,
  std::vector<Point<spacedim> >    &normal_vectors) const
{
  compute_fill (cell, n_q_points, data_set, CellSimilarity::none,
		data, quadrature_points);
  internal::compute_fill_face (*this,
			       cell, face_no, subface_no, n_q_points,
			       weights, data,
			       JxW_values, boundary_forms, normal_vectors);
}


template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::
fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
		     const unsigned int                               face_no,
		     const Quadrature<dim-1>                          &q,
		     typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
		     std::vector<Point<spacedim> >                         &quadrature_points,
		     std::vector<double>                              &JxW_values,
		     std::vector<Tensor<1,spacedim> >                      &boundary_forms,
		     std::vector<Point<spacedim> >                    &normal_vectors) const
{
				   // ensure that the following cast
				   // is really correct:
  Assert (dynamic_cast<InternalData *>(&mapping_data) != 0,
	  ExcInternalError());
  InternalData &data = static_cast<InternalData&>(mapping_data);

  const unsigned int n_q_points = q.size();

  compute_fill_face (cell, face_no, deal_II_numbers::invalid_unsigned_int,
		     n_q_points,
		     DataSetDescriptor::face (face_no,
                                              cell->face_orientation(face_no),
                                              cell->face_flip(face_no),
                                              cell->face_rotation(face_no),
                                              n_q_points),
		     q.get_weights(),
		     data,
		     quadrature_points,
		     JxW_values,
		     boundary_forms,
		     normal_vectors);
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::
fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
			const unsigned int       face_no,
			const unsigned int       sub_no,
			const Quadrature<dim-1> &q,
			typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
			std::vector<Point<spacedim> >     &quadrature_points,
			std::vector<double>          &JxW_values,
			std::vector<Tensor<1,spacedim> >  &boundary_forms,
			std::vector<Point<spacedim> >     &normal_vectors) const
{
				   // ensure that the following cast
				   // is really correct:
  Assert (dynamic_cast<InternalData *>(&mapping_data) != 0,
	  ExcInternalError());
  InternalData &data = static_cast<InternalData&>(mapping_data);

  const unsigned int n_q_points = q.size();

  compute_fill_face (cell, face_no, sub_no,
		     n_q_points,
		     DataSetDescriptor::subface (face_no, sub_no,
						 cell->face_orientation(face_no),
						 cell->face_flip(face_no),
						 cell->face_rotation(face_no),
						 n_q_points,
						 cell->subface_case(face_no)),
		     q.get_weights(),
		     data,
		     quadrature_points,
		     JxW_values,
		     boundary_forms,
		     normal_vectors);
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<1, dim> > > input,
  VectorSlice<std::vector<Tensor<1, spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
	  ExcInternalError());
  const InternalData &data = static_cast<const InternalData&>(mapping_data);

  Tensor<1, spacedim> auxiliary;

  switch (mapping_type)
    {
      case mapping_covariant:
      {
	Assert (data.update_flags & update_covariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());

	for (unsigned int i=0; i<output.size(); ++i)
	  {
	    for (unsigned int d=0;d<dim;++d)
	      auxiliary[d] = input[i][d];
	    contract (output[i], auxiliary, data.covariant[i]);
	  }
	return;
      }

      case mapping_contravariant:
      {
	Assert (data.update_flags & update_contravariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());

	for (unsigned int i=0; i<output.size(); ++i)
	  {
	    for (unsigned int d=0;d<dim;++d)
	      auxiliary[d] = input[i][d];
	    contract (output[i], data.contravariant[i], auxiliary);
	  }
	return;
      }

      case mapping_piola:
      {
	Assert (data.update_flags & update_contravariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());
	Assert (data.update_flags & update_volume_elements,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());

	for (unsigned int i=0; i<output.size(); ++i)
	  {
	    for (unsigned int d=0;d<dim;++d)
	      auxiliary[d] = input[i][d] / data.volume_elements[i];
	    contract (output[i], data.contravariant[i], auxiliary);
	  }
	return;
      }

      default:
	    Assert(false, ExcNotImplemented());
    }
}


template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<2, dim> > > input,
  VectorSlice<std::vector<Tensor<2, spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
	  ExcInternalError());
  const InternalData &data = static_cast<const InternalData&>(mapping_data);

  Tensor<2, spacedim> aux1;
  Tensor<2, spacedim> aux2;

  switch (mapping_type)
    {
      case mapping_covariant:
      {
	Assert (data.update_flags & update_covariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());

	for (unsigned int i=0; i<output.size(); ++i)
	  {
	    for (unsigned int d1=0;d1<dim;++d1)
	      for (unsigned int d2=0;d2<dim;++d2)
		aux1[d1][d2] = input[i][d1][d2];
	    contract (output[i], aux1, data.covariant[i]);
	  }
	return;
      }

      case mapping_contravariant:
      {
	Assert (data.update_flags & update_contravariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());

	for (unsigned int i=0; i<output.size(); ++i)
	  {
	    for (unsigned int d1=0;d1<dim;++d1)
	      for (unsigned int d2=0;d2<dim;++d2)
		aux1[d1][d2] = input[i][d1][d2];
	    contract (output[i], data.contravariant[i], aux1);
	  }
	return;
      }

      case mapping_covariant_gradient:
      {
	Assert (data.update_flags & update_contravariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());

	for (unsigned int i=0; i<output.size(); ++i)
	  {
	    for (unsigned int d1=0;d1<dim;++d1)
	      for (unsigned int d2=0;d2<dim;++d2)
		aux1[d1][d2] = input[i][d1][d2];

	    contract(aux2, aux1, data.covariant[i]);
	    contract(output[i], data.covariant[i], aux2);
	  }

	return;
      }

      case mapping_contravariant_gradient:
      {
	Assert (data.update_flags & update_covariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());
	Assert (data.update_flags & update_contravariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());

	for (unsigned int i=0; i<output.size(); ++i)
	  {
	    for (unsigned int d1=0;d1<dim;++d1)
	      for (unsigned int d2=0;d2<dim;++d2)
		aux1[d1][d2] = input[i][d1][d2];

	    contract(aux2, aux1, data.covariant[i]);
	    contract(output[i], data.contravariant[i], aux2);
	  }

	return;
      }

      case mapping_piola:
      {
	Assert (data.update_flags & update_covariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());
	Assert (data.update_flags & update_contravariant_transformation,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());
	Assert (data.update_flags & update_volume_elements,
		typename FEValuesBase<dim>::ExcAccessToUninitializedField());

	for (unsigned int i=0; i<output.size(); ++i)
	  {
	    for (unsigned int d1=0;d1<dim;++d1)
	      for (unsigned int d2=0;d2<dim;++d2)
		aux1[d1][d2] = input[i][d1][d2] / data.volume_elements[i];

	    contract (aux2, aux1, data.covariant[i]);
	    contract (output[i], data.contravariant[i], aux2);
	  }
	return;
      }
      default:
	    Assert(false, ExcNotImplemented());
    }
}



template<int dim, int spacedim>
Point<spacedim>
MappingQ1<dim,spacedim>::transform_unit_to_real_cell (
  const typename Triangulation<dim,spacedim>::cell_iterator& cell,
  const Point<dim>& p) const
{
				   // Use the get_data function to
				   // create an InternalData with data
				   // vectors of the right size and
				   // transformation shape values
				   // already computed at point p.
  const Quadrature<dim> point_quadrature(p);

  std::auto_ptr<InternalData>
    mdata (dynamic_cast<InternalData *> (
             get_data(update_transformation_values, point_quadrature)));

				   // compute the mapping support
				   // points
  compute_mapping_support_points(cell, mdata->mapping_support_points);

  return transform_unit_to_real_cell_internal(*mdata);
}



template<int dim, int spacedim>
Point<spacedim>
MappingQ1<dim,spacedim>::
transform_unit_to_real_cell_internal (const InternalData &data) const
{
  const unsigned int n_mapping_points=data.mapping_support_points.size();
  Assert(data.shape_values.size()==n_mapping_points, ExcInternalError());

				   // use now the InternalData to
				   // compute the point in real space.
  Point<spacedim> p_real;
  for (unsigned int i=0; i<data.mapping_support_points.size(); ++i)
    p_real += data.mapping_support_points[i] * data.shape(0,i);

  return p_real;
}



template<int dim, int spacedim>
Point<dim>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<spacedim>                            &p) const
{
				   // Let the start value of the
				   // newton iteration be the center
				   // of the unit cell
  Point<dim> p_unit;
  for (unsigned int i=0; i<dim; ++i)
    p_unit(i)=0.5;

				   // Use the get_data function to
				   // create an InternalData with data
				   // vectors of the right size and
				   // transformation shape values and
				   // derivatives already computed at
				   // point p_unit.
  const Quadrature<dim> point_quadrature(p_unit);
  std::auto_ptr<InternalData>
    mdata(dynamic_cast<InternalData *> (
            MappingQ1<dim,spacedim>::get_data(update_transformation_values
                                     | update_transformation_gradients,
                                     point_quadrature)));

  MappingQ1<dim,spacedim>::compute_mapping_support_points (cell,
							   mdata->mapping_support_points);
  Assert(mdata->mapping_support_points.size() ==
         GeometryInfo<dim>::vertices_per_cell,
	 ExcInternalError());

				   // perform the newton iteration.
  transform_real_to_unit_cell_internal(cell, p, *mdata, p_unit);

  return p_unit;
}



template<>
void
MappingQ1<1,2>::
transform_real_to_unit_cell_internal
(const Triangulation<1,2>::cell_iterator &,
 const Point<2>                            &,
 InternalData                                     &,
 Point<1>                                       &) const
{
	Assert(false, ExcNotImplemented());
}



template<>
void
MappingQ1<2,3>::
transform_real_to_unit_cell_internal
(const Triangulation<2,3>::cell_iterator &,
 const Point<3>                            &,
 InternalData                                     &,
 Point<2>                                       &) const
{
	Assert(false, ExcNotImplemented());
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell_internal
(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
 const Point<spacedim>                            &p,
 InternalData                                     &mdata,
 Point<dim>                                       &p_unit) const
{
  const unsigned int n_shapes=mdata.shape_values.size();
  Assert(n_shapes!=0, ExcInternalError());
  Assert(mdata.shape_derivatives.size()==n_shapes, ExcInternalError());

  std::vector<Point<spacedim> > &points=mdata.mapping_support_points;
  Assert(points.size()==n_shapes, ExcInternalError());

				   // Newton iteration to solve
				   // f(x)=p(x)-p=0
				   // x_{n+1}=x_n-[f'(x)]^{-1}f(x)

				   // The start value is set to be the
				   // center of the unit cell.

				   // The shape values and derivatives
				   // of the mapping at this point are
				   // previously computed.


				   // f(x)
  Point<spacedim> p_real(transform_unit_to_real_cell_internal(mdata));
  Point<spacedim> f = p_real-p;

  const double eps=1e-15*cell->diameter();
  unsigned int loop=0;
  while (f.square()>eps*eps && loop++<10)
    {
				       // f'(x)
      Tensor<2,dim> df;
      for (unsigned int k=0; k<mdata.n_shape_functions; ++k)
	{
	  const Tensor<1,dim> &grad_transform=mdata.derivative(0,k);
	  const Point<spacedim> &point=points[k];

	  for (unsigned int i=0; i<dim; ++i)
	    for (unsigned int j=0; j<dim; ++j)
	      df[i][j]+=point[i]*grad_transform[j];
	}

				       // Solve  [f'(x)]d=f(x)
      Tensor<1,dim> d;
      Tensor<2,dim> df_1;

      df_1 = invert(df);
      contract (d, df_1, static_cast<const Tensor<1,dim>&>(f));

				       // update of p_unit
      p_unit -= d;
				       // shape values and derivatives
				       // at new p_unit point
      compute_shapes(std::vector<Point<dim> > (1, p_unit), mdata);

				       // f(x)
      p_real = transform_unit_to_real_cell_internal(mdata);
      f = p_real-p;
    }
}



template<int dim, int spacedim>
Mapping<dim,spacedim> *
MappingQ1<dim,spacedim>::clone () const
{
  return new MappingQ1<dim,spacedim>(*this);
}

//---------------------------------------------------------------------------


template<int dim, int spacedim>
MappingQ1<dim,spacedim> StaticMappingQ1<dim,spacedim>::mapping;



//--------------------------- Explicit instantiations -----------------------
#include "mapping_q1.inst"


DEAL_II_NAMESPACE_CLOSE
