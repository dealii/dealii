//----------------------------  mapping_q1.cc  ---------------------------
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
//----------------------------  mapping_q1.cc  ---------------------------

#include <base/tensor.h>
#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/mapping_q1.h>
#include <fe/fe_values.h>

#include <cmath>
#include <algorithm>




template <int dim>
const unsigned int MappingQ1<dim>::n_shape_functions;



template<int dim>
MappingQ1<dim>::InternalData::InternalData (const unsigned int n_shape_functions):
		is_mapping_q1_data(true),
		n_shape_functions (n_shape_functions)
{}


template<int dim> inline
double
MappingQ1<dim>::InternalData::shape (const unsigned int qpoint,
				     const unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
	 ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0, shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}



template<int dim> inline
double&
MappingQ1<dim>::InternalData::shape (unsigned int qpoint,
				     unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_values.size(),
	 ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0, shape_values.size()));
  return shape_values [qpoint*n_shape_functions + shape_nr];
}


template<int dim> inline
Tensor<1,dim>
MappingQ1<dim>::InternalData::derivative (unsigned int qpoint,
					  unsigned int shape_nr) const
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
	 ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0, shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}



template<int dim> inline
Tensor<1,dim>&
MappingQ1<dim>::InternalData::derivative (unsigned int qpoint,
					  unsigned int shape_nr)
{
  Assert(qpoint*n_shape_functions + shape_nr < shape_derivatives.size(),
	 ExcIndexRange(qpoint*n_shape_functions + shape_nr, 0, shape_derivatives.size()));
  return shape_derivatives [qpoint*n_shape_functions + shape_nr];
}


template<int dim>
void
MappingQ1<dim>::compute_shapes (const std::vector<Point<dim> > &unit_points,
				InternalData &data) const
{
  if (data.is_mapping_q1_data)
    MappingQ1<dim>::compute_shapes_virtual(unit_points, data);
  else
    compute_shapes_virtual(unit_points, data);
}


#if (deal_II_dimension == 1)

template<>
const unsigned int MappingQ1<deal_II_dimension>::vertex_mapping[2] =
{ 0, 1
};


template<>
void
MappingQ1<1>::compute_shapes_virtual (const std::vector<Point<1> > &unit_points,
				      InternalData& data) const
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
    }
}

#endif

#if (deal_II_dimension == 2)

template<> const unsigned int
MappingQ1<2>::vertex_mapping[4] =
{
  0, 1, 3, 2
};


template<>
void
MappingQ1<2>::compute_shapes_virtual (const std::vector<Point<2> > &unit_points,
				      InternalData &data) const
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
    }
}

#endif

#if (deal_II_dimension == 3)

template<>
const unsigned int MappingQ1<deal_II_dimension>::vertex_mapping[8] =
{
  0, 1, 4, 5, 3, 2, 7, 6
};


template<>
void
MappingQ1<3>::compute_shapes_virtual (const std::vector<Point<3> > &unit_points,
				      InternalData &data) const
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
    }
}

#endif

template <int dim>
UpdateFlags
MappingQ1<dim>::update_once (const UpdateFlags in) const
{
  UpdateFlags out = UpdateFlags(in & (update_transformation_values
				      | update_transformation_gradients));

				   // Shape function values
  if (in & update_q_points)
    out |= update_transformation_values;

				   // Shape function gradients
  if (in & (update_covariant_transformation
	    | update_contravariant_transformation
	    | update_JxW_values
	    | update_gradients
	    | update_boundary_forms
	    | update_normal_vectors))
    out |= update_transformation_gradients;

  //  cerr << "Once: " << hex << out << dec << endl;

  return out;
}

template <int dim>
UpdateFlags
MappingQ1<dim>::update_each (const UpdateFlags in) const
{
				   // Select flags of concern for the
				   // transformation.
  UpdateFlags out = UpdateFlags(in & (update_q_points
				      | update_covariant_transformation
				      | update_contravariant_transformation
				      | update_JxW_values
				      | update_boundary_forms
				      | update_normal_vectors));

  //  cerr << "Mapping-each " << hex << in << ' ' << out;
  
  if (in & update_normal_vectors)
    out |= update_boundary_forms;
  
  if (in & (update_covariant_transformation
	    | update_JxW_values
	    | update_boundary_forms
	    | update_normal_vectors))
    out |= update_contravariant_transformation;

  //  cerr << "  " << hex << out << dec << endl;
  
  return out;
}


template <int dim>
void
MappingQ1<dim>::compute_data (const UpdateFlags      update_flags,
			      const Quadrature<dim> &q,
			      const unsigned int     n_original_q_points,
			      InternalData          &data) const
{
  const unsigned int npts = q.n_quadrature_points;

  data.update_once = update_once(update_flags);
  data.update_each = update_each(update_flags);
  data.update_flags = data.update_once | data.update_each;

  const UpdateFlags flags(data.update_flags);
  
  if (flags & update_transformation_values)
    data.shape_values.resize(data.n_shape_functions * npts);

  if (flags & update_transformation_gradients)
    data.shape_derivatives.resize(data.n_shape_functions * npts);

  if (flags & update_covariant_transformation)
    data.covariant.resize(n_original_q_points);

  if (flags & update_contravariant_transformation)
    data.contravariant.resize(n_original_q_points);

  compute_shapes (q.get_points(), data);
}


template <int dim>
Mapping<dim>::InternalDataBase*
MappingQ1<dim>::get_data (const UpdateFlags update_flags,
			  const Quadrature<dim>& q) const
{
  InternalData* data = new InternalData(n_shape_functions);
  compute_data (update_flags, q, q.n_quadrature_points, *data);
  return data;
}


template <int dim>
void
MappingQ1<dim>::compute_face_data (UpdateFlags update_flags,
				   const Quadrature<dim>& q,
				   const unsigned int n_original_q_points,
				   InternalData& data) const
{
  if (update_flags & update_JxW_values)
    update_flags |= update_boundary_forms;
  
   compute_data (update_flags, q, n_original_q_points, data);

#if (deal_II_dimension>1)
  if (data.update_flags & update_boundary_forms)
    {
      data.aux.resize(dim-1);
      for (unsigned int i=0;i<dim-1;++i)
	data.aux[i].resize(n_original_q_points);
      
				       // Compute tangentials to the
				       // unit cell.
      const unsigned int nfaces = GeometryInfo<dim>::faces_per_cell;
      data.unit_tangentials.resize(nfaces*(dim-1));
      for (unsigned int i=0;i<nfaces;++i)
	{
					   // Base index of the
					   // non-zero entry of the
					   // tangential vector.  Add
					   // dim so we can subtract 1
					   // without getting negative
					   // values.
	  unsigned int nindex = normal_directions[i]/2 + dim;

					   // First tangential has a
					   // non-zero in component
					   // (i+1)%dim, if normal is
					   // non-zero in i.
	  Tensor<1,dim> tangential;
	  tangential[(nindex+1)%dim] = (normal_directions[i]%2) ? -1 : 1;
	  data.unit_tangentials[i].resize(n_original_q_points);
	  fill (data.unit_tangentials[i].begin(),
		data.unit_tangentials[i].end(),
		tangential);
	  
	  if (dim>2)
	    {
					       // Second tangential
					       // has a non-zero in
					       // component (i-1)%dim,
					       // if normal is
					       // non-zero in
					       // i. Creates a
					       // right-handed system.
	      Tensor<1,dim> tangential;
	      tangential[(nindex-1)%dim] = 1.;
	      data.unit_tangentials[i+nfaces].resize(n_original_q_points);
	      fill (data.unit_tangentials[i+nfaces].begin(),
		    data.unit_tangentials[i+nfaces].end(),
		    tangential);
	    }
	}
    }
#endif
}

  

template <int dim>
Mapping<dim>::InternalDataBase*
MappingQ1<dim>::get_face_data (const UpdateFlags        update_flags,
			       const Quadrature<dim-1> &quadrature) const
{
  InternalData* data = new InternalData(n_shape_functions);
  QProjector<dim> q (quadrature, false);
  compute_face_data (update_flags, q, quadrature.n_quadrature_points, *data);

  return data;
}



template <int dim>
Mapping<dim>::InternalDataBase*
MappingQ1<dim>::get_subface_data (const UpdateFlags update_flags,
				  const Quadrature<dim-1>& quadrature) const
{
  InternalData* data = new InternalData(n_shape_functions);
  QProjector<dim> q (quadrature, true);
  compute_face_data (update_flags, q, quadrature.n_quadrature_points, *data);

  return data;
}




template <int dim>
void
MappingQ1<dim>::compute_fill (const typename DoFHandler<dim>::cell_iterator &cell,
			      const unsigned int   npts,
			      const unsigned int   offset,
			      InternalData        &data,
			      std::vector<Point<dim> > &quadrature_points) const
{
  const UpdateFlags update_flags(data.current_update_flags());

  if (update_flags & update_q_points)
    {
      Assert (quadrature_points.size() == npts,
	      ExcDimensionMismatch(quadrature_points.size(), npts));
      fill(quadrature_points.begin(),
	   quadrature_points.end(),
	   Point<dim>());
    }

  if (update_flags & update_covariant_transformation)
    {
      Assert (data.covariant.size() == npts,
	      ExcDimensionMismatch(data.covariant.size(), npts));
    }
  
  if (update_flags & update_contravariant_transformation)
    {
      Assert (data.contravariant.size() == npts,
	      ExcDimensionMismatch(data.contravariant.size(), npts));
      fill(data.contravariant.begin(),
	   data.contravariant.end(),
	   Tensor<2,dim>());
    }
  
  if (update_flags & update_jacobian_grads)
    {
      Assert(false, ExcNotImplemented());
//      Assert (covariant_grads.size () == npts,
//	      ExcDimensionMismatch(covariant_grads.size(), npts));
    }
  
  std::vector<Point<dim> > &a=data.mapping_support_points;
  
				   // store all Lagrangian
				   // support points in a
  if (a.size()==0
      || (&cell->get_triangulation() !=
	  &data.cell_of_current_support_points->get_triangulation())
      || (cell!=data.cell_of_current_support_points))
    {
      compute_mapping_support_points(cell, a);
      data.cell_of_current_support_points=cell;
    }
  
  for (unsigned int point=0; point<npts; ++point)
    {
				       // First, compute function
				       // values and derivatives
      for (unsigned int k=0; k<data.n_shape_functions; ++k)
	{
	  if (update_flags & update_q_points)
	    {
	      quadrature_points[point]
		+= data.shape(point+offset,k) * a[k];
	    }
	  if (update_flags & update_contravariant_transformation)
	    {
	      for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=0; j<dim; ++j)
		  data.contravariant[point][i][j]
		    += data.derivative(point+offset, k)[j]
		    * a[k][i];
	    }
	}

				       // Invert contravariant for
				       // covariant transformation
				       // matrices
      if (update_flags & update_covariant_transformation)
	data.covariant[point]
	  = invert(data.contravariant[point]);
    }
  data.first_cell = false;
}


template <int dim>
void
MappingQ1<dim>::compute_mapping_support_points(
  const typename Triangulation<dim>::cell_iterator &cell,
  std::vector<Point<dim> > &a) const
{
  a.resize(GeometryInfo<dim>::vertices_per_cell);

  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a[i] = cell->vertex(vertex_mapping[i]);
}



template <int dim>
void
MappingQ1<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				const Quadrature<dim>                &q,
				Mapping<dim>::InternalDataBase      &mapping_data,
				std::vector<Point<dim> >                  &quadrature_points,
				std::vector<double>                       &JxW_values) const
{
  InternalData *data_ptr = dynamic_cast<InternalData *> (&mapping_data);
  Assert(data_ptr!=0, ExcInternalError());
  InternalData &data=*data_ptr;

  const unsigned int npts=q.n_quadrature_points;
  
  compute_fill (cell,
		npts,
		0,
		data,
		quadrature_points);

  
  const UpdateFlags update_flags(data.current_update_flags());
  const std::vector<double> &weights=q.get_weights();

				   // Multiply quadrature weights
                                  // by Jaconian determinants
 //TODO: compute Jacobi determinants directly, if co/contravariant is not needed
  if (update_flags & update_JxW_values)
    {      
      Assert (JxW_values.size() == npts,
	      ExcDimensionMismatch(JxW_values.size(), npts));
      for (unsigned int point=0; point<npts; ++point)
       JxW_values[point]
	 = determinant(data.contravariant[point])*weights[point];
    }
}



template <int dim>
void
MappingQ1<dim>::compute_fill_face (const typename DoFHandler<dim>::cell_iterator &cell,
				   const unsigned int      face_no,
				   const bool              is_subface,
				   const unsigned int      npts,
				   const unsigned int      offset,
				   const std::vector<double>   &weights,
				   InternalData           &data,
				   std::vector<Point<dim> >    &quadrature_points,
				   std::vector<double>         &JxW_values,
				   std::vector<Tensor<1,dim> > &boundary_forms,
				   std::vector<Point<dim> >    &normal_vectors) const
{
  compute_fill (cell,
		npts,
		offset,
		data,
		quadrature_points);

  const UpdateFlags update_flags(data.current_update_flags());
  
  if (update_flags & update_boundary_forms)
    {
      Assert (boundary_forms.size()==npts,
	      ExcDimensionMismatch(boundary_forms.size(), npts));
      if (update_flags & update_normal_vectors)
	Assert (normal_vectors.size()==npts,
		ExcDimensionMismatch(normal_vectors.size(), npts));
      if (update_flags & update_JxW_values)
	Assert (JxW_values.size() == npts,
		ExcDimensionMismatch(JxW_values.size(), npts));
      
      
      transform_contravariant(data.aux[0],
			      data.unit_tangentials[face_no],
			      data, 0);

      typename std::vector<Tensor<1,dim> >::iterator
	result = boundary_forms.begin();
      typename std::vector<Tensor<1,dim> >::const_iterator
	end = boundary_forms.end();
      typename std::vector<Tensor<1,dim> >::const_iterator
	tang1 = data.aux[0].begin();
      
      switch (dim)
	{
	  case 2:
	  {
	    for (; result != end; ++result, ++tang1)
	      cross_product (*result, *tang1);
	    break;
	  };

	  case 3:
	  {
	    transform_contravariant(data.aux[1],
				    data.unit_tangentials[
				      face_no+GeometryInfo<dim>::faces_per_cell],
				    data, 0);
	    typename std::vector<Tensor<1,dim> >::const_iterator
	      tang2 = data.aux[1].begin();
	    for (;result != end; ++result, ++tang1, ++tang2)
	      cross_product (*result, *tang1, *tang2);
	    break;
	  };

	  default:
		Assert(false, ExcNotImplemented());
	}
      
      if (update_flags & (update_normal_vectors
			  | update_JxW_values))
	for (unsigned int i=0;i<boundary_forms.size();++i)
	  {
	    double f = sqrt(contract(boundary_forms[i],
				     boundary_forms[i]));
	    if (update_flags & update_JxW_values)
	      {
		JxW_values[i] = f * weights[i];
		if (is_subface)
		  JxW_values[i]/=GeometryInfo<dim>::subfaces_per_face;
	      }
	    if (update_flags & update_normal_vectors)
	      {
		normal_vectors[i] = boundary_forms[i];
		normal_vectors[i] /= f;
	      }
	  }
    }
}


template <int dim>
void
MappingQ1<dim>::fill_fe_face_values (const typename DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int       face_no,
				     const Quadrature<dim-1> &q,
				     typename Mapping<dim>::InternalDataBase &mapping_data,
				     std::vector<Point<dim> >     &quadrature_points,
				     std::vector<double>          &JxW_values,
				     std::vector<Tensor<1,dim> >  &boundary_forms,
				     std::vector<Point<dim> >     &normal_vectors) const
{
  InternalData *data_ptr = dynamic_cast<InternalData *> (&mapping_data);
  Assert(data_ptr!=0, ExcInternalError());
  InternalData &data=*data_ptr;

  const unsigned int npts=q.n_quadrature_points;
  const unsigned int offset=face_no*npts;
  
  compute_fill_face (cell, face_no, false,
		     npts,
		     offset,
		     q.get_weights(),
		     data,
		     quadrature_points,
		     JxW_values,
		     boundary_forms,
		     normal_vectors);
}


template <int dim>
void
MappingQ1<dim>::fill_fe_subface_values (const typename DoFHandler<dim>::cell_iterator &cell,
					const unsigned int       face_no,
					const unsigned int       sub_no,
					const Quadrature<dim-1> &q,
					typename Mapping<dim>::InternalDataBase &mapping_data,
					std::vector<Point<dim> >     &quadrature_points,
					std::vector<double>          &JxW_values,
					std::vector<Tensor<1,dim> >  &boundary_forms,
					std::vector<Point<dim> >     &normal_vectors) const
{
  InternalData *data_ptr = dynamic_cast<InternalData *> (&mapping_data);
  Assert(data_ptr!=0, ExcInternalError());
  InternalData &data=*data_ptr;

  const unsigned int npts=q.n_quadrature_points;
  const unsigned int offset=
    (face_no*GeometryInfo<dim>::subfaces_per_face + sub_no)*npts;
  
  compute_fill_face (cell, face_no, true,
		     npts,
		     offset,
		     q.get_weights(),
		     data,
		     quadrature_points,
		     JxW_values,
		     boundary_forms,
		     normal_vectors);
}


#if (deal_II_dimension == 1)

template <>
void
MappingQ1<1>::compute_fill_face (const DoFHandler<1>::cell_iterator &,
				 const unsigned int,
				 const bool,
				 const unsigned int,
				 const unsigned int,
				 const std::vector<double> &,
				 InternalData &,
				 std::vector<Point<1> > &,
				 std::vector<double> &,
				 std::vector<Tensor<1,1> > &,
				 std::vector<Point<1> > &) const
{
  Assert(false, ExcNotImplemented());
}


template <>
void
MappingQ1<1>::fill_fe_face_values (const DoFHandler<1>::cell_iterator &,
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
MappingQ1<1>::fill_fe_subface_values (const DoFHandler<1>::cell_iterator &,
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
MappingQ1<dim>::transform_covariant (std::vector<Tensor<1,dim> >       &dst,
				     const std::vector<Tensor<1,dim> > &src,
				     const Mapping<dim>::InternalDataBase &mapping_data,
				     const unsigned int src_offset) const
{
  covariant_transformation(dst, src, mapping_data, src_offset);
}


template <int dim>
void
MappingQ1<dim>::transform_covariant (std::vector<Point<dim> >       &dst,
				     const std::vector<Point<dim> > &src,
				     const Mapping<dim>::InternalDataBase &mapping_data,
				     const unsigned int src_offset) const
{
  covariant_transformation(dst, src, mapping_data, src_offset);
}

template <int dim>
void
MappingQ1<dim>::transform_contravariant (std::vector<Tensor<1,dim> >       &dst,
					 const std::vector<Tensor<1,dim> > &src,
					 const Mapping<dim>::InternalDataBase &mapping_data,
					 const unsigned int src_offset) const
{
  contravariant_transformation(dst, src, mapping_data, src_offset);
}


template <int dim>
void
MappingQ1<dim>::transform_contravariant (std::vector<Point<dim> >       &dst,
					 const std::vector<Point<dim> > &src,
					 const Mapping<dim>::InternalDataBase &mapping_data,
					 const unsigned int src_offset) const
{
  contravariant_transformation(dst, src, mapping_data, src_offset);
}



template <int dim>
Point<dim> MappingQ1<dim>::transform_unit_to_real_cell (
  const typename Triangulation<dim>::cell_iterator cell,
  const Point<dim> &p) const
{
  				   // Use the get_cell_data function
				   // to create an InternalData with
				   // data vectors already of the
				   // right size and with mapping
				   // support points set.
  InternalData *mdata=get_cell_data(cell, update_transformation_values);
  Assert(mdata!=0, ExcInternalError());

  return transform_unit_to_real_cell_internal(cell, p, *mdata);
}


template <int dim>
Point<dim> MappingQ1<dim>::transform_unit_to_real_cell_internal (
  const typename Triangulation<dim>::cell_iterator cell,
  const Point<dim> &p,
  const InternalData &data) const
{
  const unsigned int n_mapping_points=data.mapping_support_points.size();
  Assert(data.shape_values.size()==n_mapping_points, ExcInternalError());
  
				   // use now the InternalData to
				   // compute the point in real space.
  Point<dim> p_real;
  for (unsigned int i=0; i<data.mapping_support_points.size(); ++i)
    p_real+=data.mapping_support_points[i] * data.shape(0,i);

  return p_real;
}



template <int dim>
Point<dim> MappingQ1<dim>::transform_real_to_unit_cell (
  const typename Triangulation<dim>::cell_iterator cell,
  const Point<dim> &p) const
{
				   // Use the get_cell_data function
				   // to create an InternalData with
				   // data vectors already of the
				   // right size and with mapping
				   // support points set.
  InternalData *mdata=get_cell_data(cell,
				    update_transformation_values
				    | update_transformation_gradients);
  Assert(mdata!=0, ExcInternalError());
  std::vector<Point<dim> > &points=mdata->mapping_support_points;
  
				   // Newton iteration to solve
				   // f(x)=p(x)-p=0
				   // x_{n+1}=x_n-[f'(x)]^{-1}f(x)
  
				   // Let the start value be the
				   // center of the unit cell
				   // (p_unit stands for x)
  Point<dim> p_unit;
  for (unsigned int i=0; i<dim; ++i)
    p_unit(i)=0.5;
  
				   // compute shape values and
				   // derivatives of the mapping
  compute_shapes(std::vector<Point<dim> > (1, p_unit), *mdata);
  
				   // f(x)
  Point<dim> p_real(transform_unit_to_real_cell_internal(cell, p_unit, *mdata));
  Point<dim> f = p_real-p;

  const double eps=1e-15*cell->diameter();
  unsigned int loop=0;
  while (f.square()>eps*eps && loop++<10)
    {      
				       // f'(x)
      Tensor<2,dim> df;
      for (unsigned int k=0; k<mdata->n_shape_functions; ++k)
	{
	  const Tensor<1,dim> &grad_transform=mdata->derivative(0,k);
	  const Point<dim> &point=points[k];
	  
	  for (unsigned int i=0; i<dim; ++i)
	    for (unsigned int j=0; j<dim; ++j)
	      df[i][j]+=point[i]*grad_transform[j];
	}
      
				       // Solve  [f'(x)]d=f(x)
      Point<dim> d;
      Tensor<2,dim> df_1;

      df_1 = invert(df);
      contract (d, df_1, f);

				       // update of p_unit
      p_unit -= d;
				       // shape values and derivatives
				       // at new p_unit point
      compute_shapes(std::vector<Point<dim> > (1, p_unit), *mdata);
      
				       // f(x)
      p_real=transform_unit_to_real_cell_internal(cell, p_unit, *mdata);
      f = p_real-p;
    }

  delete mdata;
  return p_unit;
}


template <int dim>
MappingQ1<dim>::InternalData*
MappingQ1<dim>::get_cell_data(const typename Triangulation<dim>::cell_iterator cell,
			      const UpdateFlags update_flags) const
{
  static Point<dim> dummy_p;
  static Quadrature<dim> dummy_quadrature(dummy_p);

  InternalData *mdata=dynamic_cast<InternalData *> (
    get_data(update_flags, dummy_quadrature));
  Assert(mdata!=0, ExcInternalError());

				   // compute the mapping support
				   // points
  std::vector<Point<dim> > &points=mdata->mapping_support_points;
  compute_mapping_support_points(cell, points);

  return mdata;
}

  



template <int dim>
template <typename tensor_>
inline
void
MappingQ1<dim>::contravariant_transformation (std::vector<tensor_>       &dst,
					      const std::vector<tensor_> &src,
					      const Mapping<dim>::InternalDataBase &mapping_data,
					      const unsigned int src_offset) const
{
  Assert(tensor_::dimension==dim && tensor_::rank==1, Mapping<dim>::ExcInvalidData());
  const InternalData* data_ptr = dynamic_cast<const InternalData *> (&mapping_data);
  Assert(data_ptr!=0, ExcInternalError());
  const InternalData &data=*data_ptr;

  Assert (data.update_flags & update_contravariant_transformation,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());
  
  Assert (src.size() + src_offset >= data.contravariant.size(),
	  ExcDimensionMismatch(src.size(), data.contravariant.size()));
  Assert (dst.size() == data.contravariant.size(),
	  ExcDimensionMismatch(dst.size(), data.contravariant.size()));

  typename std::vector<tensor_>::const_iterator vec = src.begin()+src_offset;
  typename std::vector<Tensor<2,dim> >::const_iterator tensor = data.contravariant.begin();
  typename std::vector<tensor_>::iterator result = dst.begin();
  typename std::vector<tensor_>::const_iterator end = dst.end();
  
  while (result!=end)
    {
      contract (*(result++), *(tensor++), *(vec++));
    }
}


template <int dim>
template <typename tensor_>
inline
void
MappingQ1<dim>::covariant_transformation (std::vector<tensor_>       &dst,
					  const std::vector<tensor_> &src,
					  const Mapping<dim>::InternalDataBase &mapping_data,
					  const unsigned int src_offset) const
{
  Assert(tensor_::dimension==dim && tensor_::rank==1, Mapping<dim>::ExcInvalidData());
  const InternalData *data_ptr = dynamic_cast<const InternalData *> (&mapping_data);
  Assert(data_ptr!=0, ExcInternalError());
  const InternalData &data=*data_ptr;

  Assert (data.update_flags & update_covariant_transformation,
	  typename FEValuesBase<dim>::ExcAccessToUninitializedField());
  
  Assert (src.size() + src_offset >= data.contravariant.size(),
	  ExcDimensionMismatch(src.size() + src_offset, data.contravariant.size()));
  Assert (dst.size() == data.contravariant.size(),
	  ExcDimensionMismatch(dst.size() + src_offset, data.contravariant.size()));

  typename std::vector<tensor_>::const_iterator vec = src.begin() + src_offset;
  typename std::vector<Tensor<2,dim> >::const_iterator tensor = data.covariant.begin();
  typename std::vector<tensor_>::iterator result = dst.begin();
  const typename std::vector<tensor_>::const_iterator end = dst.end();
  
  while (result!=end)
    {
      contract (*(result++), *(vec++), *(tensor++));
    }
}



//----------------------------------------------------------------------//

template class MappingQ1<deal_II_dimension>;

template void MappingQ1<deal_II_dimension>::contravariant_transformation (
  std::vector<Tensor<1,deal_II_dimension> >       &dst,
  const std::vector<Tensor<1,deal_II_dimension> > &src,
  const Mapping<deal_II_dimension>::InternalDataBase& internal,
  const unsigned int src_offset) const;

template void MappingQ1<deal_II_dimension>::contravariant_transformation (
  std::vector<Point<deal_II_dimension> >       &dst,
  const std::vector<Point<deal_II_dimension> > &src,
  const Mapping<deal_II_dimension>::InternalDataBase& internal,
  const unsigned int src_offset) const;

template void MappingQ1<deal_II_dimension>::covariant_transformation (
  std::vector<Tensor<1,deal_II_dimension> >       &dst,
  const std::vector<Tensor<1,deal_II_dimension> > &src,
  const Mapping<deal_II_dimension>::InternalDataBase& internal,
  const unsigned int src_offset) const;

template void MappingQ1<deal_II_dimension>::covariant_transformation (
  std::vector<Point<deal_II_dimension> >       &dst,
  const std::vector<Point<deal_II_dimension> > &src,
  const Mapping<deal_II_dimension>::InternalDataBase& internal,
  const unsigned int src_offset) const;

