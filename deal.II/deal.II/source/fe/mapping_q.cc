//----------------------------  mapping_q.cc  ---------------------------
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
//----------------------------  mapping_q.cc  ---------------------------

#include <fe/mapping_q.h>
#include <fe/fe_q.h>
#include <base/quadrature.h>
#include <base/quadrature_lib.h>
#include <base/tensor_product_polynomials.h>
#include <lac/full_matrix.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <dofs/dof_accessor.h>
#include <fe/mapping_q1.h>

#include <numeric>



template<int dim>
MappingQ<dim>::InternalData::InternalData (const unsigned int n_shape_functions):
		MappingQ1<dim>::InternalData(n_shape_functions),
//TODO: in 1d, use_mapping_q1_on_current_cell is always false, but should be true.
		use_mapping_q1_on_current_cell(false),
		mapping_q1_data(1 << dim)
{
  is_mapping_q1_data=false;
}


#if deal_II_dimension == 1

// in 1d, it is irrelevant which polynomial degree to use, since all
// cells are scaled linearly
template<>
MappingQ<1>::MappingQ (const unsigned int):
		laplace_on_quad_vector(0),
		laplace_on_hex_vector(0),
		degree(1),
		n_inner(0),
		n_outer(0),
		tensor_pols(0),
		n_shape_functions(2),
		renumber(0),
//TODO: why have two ways to compute? if they both work, choose one and remove the other    
		alternative_normals_computation(false),
//TODO: remove use_mapping_q_on_all_cells as it is set to false in the constructor and never set again    
		use_mapping_q_on_all_cells(false)
{}


template<>
MappingQ<1>::~MappingQ ()
{}

#endif



template<typename number>
static number power(const number x, const unsigned int y)
{
  number value=1;
  for (unsigned int i=0; i<y; ++i)
    value*=x;
  return value;
}


template<int dim>
MappingQ<dim>::MappingQ (const unsigned int p):
		laplace_on_quad_vector(0),
		laplace_on_hex_vector(0),
		degree(p),
		n_inner(power(degree-1, dim)),
		n_outer((dim==2) ? 4+4*(degree-1)
			:8+12*(degree-1)+6*(degree-1)*(degree-1)),
		tensor_pols(0),
		n_shape_functions(0),
		renumber(0),
//TODO: why have two ways to compute? if they both work, choose one and remove the other    
		alternative_normals_computation(false),
//TODO: remove use_mapping_q_on_all_cells as it is set to false in the constructor and never set again    
		use_mapping_q_on_all_cells(false)
{
				   // Construct the tensor product
				   // polynomials used as shape
				   // functions for the Qp mapping of
				   // cells at the boundary.
  std::vector<LagrangeEquidistant> v;
  for (unsigned int i=0; i<=degree; ++i)
    v.push_back(LagrangeEquidistant(degree,i));

  tensor_pols = new TensorProductPolynomials<dim> (v);
  n_shape_functions=tensor_pols->n_tensor_product_polynomials();
  Assert(n_inner+n_outer==n_shape_functions, ExcInternalError());
  
				   // build the renumbering of the
				   // shape functions of the Qp
				   // mapping.
  renumber.resize(n_shape_functions,0);
  std::vector<unsigned int> dpo(dim+1, 1);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i] = dpo[i-1]*(degree-1);
  FiniteElementData<dim> fe_data(dpo, 1);
  FE_Q<dim>::build_renumbering (fe_data, p, renumber);

				   // build laplace_on_quad_vector
  if (degree>1)
    {
      if (dim >= 2)
	set_laplace_on_quad_vector(laplace_on_quad_vector);
      if (dim >= 3)
	set_laplace_on_hex_vector(laplace_on_hex_vector);
    }
}


template<int dim>
MappingQ<dim>::~MappingQ ()
{
  delete tensor_pols;
}



#if deal_II_dimension == 1

template<>
void
MappingQ<1>::compute_shapes_virtual (const std::vector<Point<1> > &unit_points,
				     MappingQ1<1>::InternalData   &data) const
{
  MappingQ1<1>::compute_shapes_virtual(unit_points, data);
}

#endif



template<int dim>
void
MappingQ<dim>::compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
				       MappingQ1<dim>::InternalData &data) const
{
  
  const unsigned int n_points=unit_points.size();
  std::vector<double> values;
  std::vector<Tensor<1,dim> > grads;
  if (data.shape_values.size()!=0)
    {
      Assert(data.shape_values.size()==n_shape_functions*n_points,
	     ExcInternalError());
      values.resize(n_shape_functions);
    }
  if (data.shape_derivatives.size()!=0)
    {
      Assert(data.shape_derivatives.size()==n_shape_functions*n_points,
	     ExcInternalError());
      grads.resize(n_shape_functions);
    }
  
				   // dummy variable of size 0
  std::vector<Tensor<2,dim> > grad2;

  
  if (data.shape_values.size()!=0 || data.shape_derivatives.size()!=0)
    for (unsigned int point=0; point<n_points; ++point)
      {
	tensor_pols->compute(unit_points[point], values, grads, grad2);
	
	if (data.shape_values.size()!=0)
	  for (unsigned int i=0; i<n_shape_functions; ++i)
	    data.shape(point,renumber[i]) = values[i];
	
	if (data.shape_derivatives.size()!=0)
	  for (unsigned int i=0; i<n_shape_functions; ++i)
	    data.derivative(point,renumber[i]) = grads[i];
      }
}


template <int dim>
UpdateFlags
MappingQ<dim>::update_each (const UpdateFlags in) const
{
  UpdateFlags out=MappingQ1<dim>::update_each(in);

  if (in & update_normal_vectors)
    if (alternative_normals_computation)
      out |= update_covariant_transformation;

  return out;
}



template <int dim>
Mapping<dim>::InternalDataBase*
MappingQ<dim>::get_data (const UpdateFlags update_flags,
			 const Quadrature<dim> &quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  compute_data (update_flags, quadrature,
		quadrature.n_quadrature_points, *data);
  if (!use_mapping_q_on_all_cells)
    compute_data (update_flags, quadrature,
		  quadrature.n_quadrature_points, data->mapping_q1_data);
  return data;
}



template <int dim>
Mapping<dim>::InternalDataBase*
MappingQ<dim>::get_face_data (const UpdateFlags update_flags,
			      const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  QProjector<dim> q (quadrature, false);
  compute_face_data (update_flags, q,
		     quadrature.n_quadrature_points, *data);
  if (!use_mapping_q_on_all_cells)
    MappingQ1<dim>::compute_face_data (update_flags, q,
				       quadrature.n_quadrature_points,
				       data->mapping_q1_data);

  return data;
}



template <int dim>
Mapping<dim>::InternalDataBase*
MappingQ<dim>::get_subface_data (const UpdateFlags update_flags,
				 const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  QProjector<dim> q (quadrature, true);
  compute_face_data (update_flags, q,
		     quadrature.n_quadrature_points, *data);
  if (!use_mapping_q_on_all_cells)
    MappingQ1<dim>::compute_face_data (update_flags, q,
				       quadrature.n_quadrature_points,
				       data->mapping_q1_data);

  return data;
}



template <int dim>
void
MappingQ<dim>::compute_face_data (UpdateFlags update_flags,
				  const Quadrature<dim>& q,
				  const unsigned int n_original_q_points,
				  MappingQ1<dim>::InternalData& mapping_q1_data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &data = dynamic_cast<InternalData&> (mapping_q1_data);
  
  MappingQ1<dim>::compute_face_data(update_flags, q,
				    n_original_q_points, data);

//TODO: externalize this to a proper template function. or: scrap the whole function for 1d since it doesn't make much sense anyway?
#if (deal_II_dimension>1)
  if ((data.update_flags & update_normal_vectors)
      && alternative_normals_computation)
    {
      const unsigned int nfaces = GeometryInfo<dim>::faces_per_cell;
      data.unit_normals.resize(nfaces);
      std::vector<Tensor<1,dim> > n(nfaces);
      if (dim==2)
	{
	  n[0][1]=-1;
	  n[1][0]=1;
	  n[2][1]=1;
	  n[3][0]=-1;
	}
      else if (dim==3)
	{
	  n[0][1]=-1;
	  n[1][1]=1;
	  n[2][2]=-1;
	  n[3][0]=1;
	  n[4][2]=1;
	  n[5][0]=-1;		  
	}
      else
	Assert(false, ExcNotImplemented());
      
      for (unsigned int i=0; i<nfaces; ++i)
	{
	  data.unit_normals[i].resize(n_original_q_points);
	  fill (data.unit_normals[i].begin(),
		data.unit_normals[i].end(),
		n[i]);
	}
    }
#endif
}



template <int dim>
void
MappingQ<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
			       const Quadrature<dim>                &q,
			       Mapping<dim>::InternalDataBase       &mapping_data,
			       std::vector<Point<dim> >             &quadrature_points,
			       std::vector<double>                  &JxW_values) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &data = dynamic_cast<InternalData&> (mapping_data);

  data.use_mapping_q1_on_current_cell=!(use_mapping_q_on_all_cells
					|| cell->has_boundary_lines());
  
  if (data.use_mapping_q1_on_current_cell)
    MappingQ1<dim>::fill_fe_values(cell, q, data.mapping_q1_data,
				   quadrature_points, JxW_values);
  else
    MappingQ1<dim>::fill_fe_values(cell, q, data,
				   quadrature_points, JxW_values);
}


template <int dim>
void
MappingQ<dim>::fill_fe_face_values (const typename DoFHandler<dim>::cell_iterator &cell,
				    const unsigned int       face_no,
				    const Quadrature<dim-1> &q,
				    typename Mapping<dim>::InternalDataBase &mapping_data,
				    std::vector<Point<dim> >     &quadrature_points,
				    std::vector<double>          &JxW_values,
				    std::vector<Tensor<1,dim> >  &exterior_forms,
				    std::vector<Point<dim> >     &normal_vectors) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &data = dynamic_cast<InternalData&> (mapping_data);
  
//TODO: shouldn't we ask whether the face is at the boundary, rather than the cell?  
  data.use_mapping_q1_on_current_cell=!(use_mapping_q_on_all_cells
					|| cell->has_boundary_lines());

  const unsigned int npts=q.n_quadrature_points;
  const unsigned int offset=face_no*npts;

  if (data.use_mapping_q1_on_current_cell)
    MappingQ1<dim>::compute_fill_face (cell, face_no, false,
				       npts, offset, q.get_weights(),
				       data.mapping_q1_data,
				       quadrature_points, JxW_values,
				       exterior_forms, normal_vectors);
  else
    compute_fill_face (cell, face_no, false,
		       npts, offset, q.get_weights(),
		       data,
		       quadrature_points, JxW_values,
		       exterior_forms, normal_vectors);
}


template <int dim>
void
MappingQ<dim>::fill_fe_subface_values (const typename DoFHandler<dim>::cell_iterator &cell,
				       const unsigned int       face_no,
				       const unsigned int       sub_no,
				       const Quadrature<dim-1> &q,
				       typename Mapping<dim>::InternalDataBase &mapping_data,
				       std::vector<Point<dim> >     &quadrature_points,
				       std::vector<double>          &JxW_values,
				       std::vector<Tensor<1,dim> >  &exterior_forms,
				       std::vector<Point<dim> >     &normal_vectors) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &data = dynamic_cast<InternalData&> (mapping_data);

//TODO: shouldn't we ask whether the face is at the boundary, rather than the cell?  
  data.use_mapping_q1_on_current_cell=!(use_mapping_q_on_all_cells
					|| cell->has_boundary_lines());

  const unsigned int npts=q.n_quadrature_points;
  const unsigned int offset=
    (face_no*GeometryInfo<dim>::subfaces_per_face + sub_no)*npts;

  if (data.use_mapping_q1_on_current_cell)
    MappingQ1<dim>::compute_fill_face (cell, face_no, true,
				       npts, offset, q.get_weights(),
				       data.mapping_q1_data,
				       quadrature_points, JxW_values,
				       exterior_forms, normal_vectors);
  else
    compute_fill_face (cell, face_no, true,
		       npts, offset, q.get_weights(),
		       data,
		       quadrature_points, JxW_values,
		       exterior_forms, normal_vectors);
}


#if deal_II_dimension==1

template <>
void
MappingQ<1>::set_laplace_on_quad_vector(std::vector<std::vector<double> > &) const
{
  Assert(false, ExcInternalError());
}

#else

template <int dim>
void
MappingQ<dim>::set_laplace_on_quad_vector(std::vector<std::vector<double> > &loqvs) const
{
  Assert(degree>1, ExcInternalError());
  const unsigned int n_inner_2d=(degree-1)*(degree-1);
  const unsigned int n_outer_2d=4+4*(degree-1);

				   // first check whether we have
				   // precomputed the values for some
				   // polynomial degree
  double const *loqv_ptr=0;
  if (degree==2)
    {
      static const double loqv2[1*8]
	={1/16., 1/16., 1/16., 1/16., 3/16., 3/16., 3/16., 3/16.};
      loqv_ptr=&loqv2[0];
    }
  else if (degree==3)
    {
      static const double loqv3[4*12]
	={80/1053., 1/81., 11/1053., 1/81., 25/117., 44/351.,
	  7/117., 16/351., 7/117., 16/351., 25/117., 44/351.,
	  1/81., 80/1053., 1/81., 11/1053., 44/351., 25/117.,
	  25/117., 44/351., 16/351., 7/117., 7/117., 16/351.,
	  1/81., 11/1053., 1/81., 80/1053., 7/117., 16/351.,
	  16/351., 7/117., 25/117., 44/351., 44/351., 25/117.,
	  11/1053., 1/81., 80/1053., 1/81., 16/351., 7/117.,
	  44/351., 25/117., 44/351., 25/117., 16/351., 7/117.};
      
      loqv_ptr=&loqv3[0];
    }

  if (loqv_ptr!=0)
    {
				       // precomputed. copy values to
				       // the loqvs array
      loqvs.resize(n_inner_2d);
      for (unsigned int unit_point=0; unit_point<n_inner_2d; ++unit_point)
	{
	  loqvs[unit_point].resize(n_outer_2d, 0);
	  for (unsigned int k=0; k<n_outer_2d; ++k)
	    loqvs[unit_point][k]=loqv_ptr[unit_point*n_outer_2d+k];
	}
    }
  else
    {
				   // not precomputed, then do so now
      if (dim==2)
	compute_laplace_vector(loqvs);
//TODO: what if dim==3?
    }

				   // the sum of weights of the points
				   // at the outer rim should be
				   // one. check this
  for (unsigned int unit_point=0; unit_point<loqvs.size(); ++unit_point)
    Assert(fabs(std::accumulate(loqvs[unit_point].begin(),
				loqvs[unit_point].end(),0.)-1)<1e-13,
	   ExcInternalError());
  
				   // TEST output
  if (false)
    {
      std::cout << "degree=" << degree << std::endl;
      for (unsigned int unit_point=0; unit_point<loqvs.size(); ++unit_point)
	for (unsigned int k=0; k<n_outer_2d; ++k)
	  std::cout << loqvs[unit_point][k] << std::endl;
    }
}

#endif


#if deal_II_dimension==3

template <>
void
MappingQ<3>::set_laplace_on_hex_vector(std::vector<std::vector<double> > &lohvs) const
{
  Assert(degree>1, ExcInternalError());

				   // first check whether we have
				   // precomputed the values for some
				   // polynomial degree
  double const *lohv_ptr=0;
  if (degree==2)
    {
      static const double loqv2[26]
  	={1/128., 1/128., 1/128., 1/128., 1/128., 1/128., 1/128., 1/128.,
	  7/192., 7/192., 7/192., 7/192., 7/192., 7/192., 7/192., 7/192.,
	  7/192., 7/192., 7/192., 7/192.,
	  1/12., 1/12., 1/12., 1/12., 1/12., 1/12.};
      
      lohv_ptr=&loqv2[0];
    }
  
  if (lohv_ptr!=0)
    {
				       // precomputed. copy values to
				       // the lohvs array
      lohvs.resize(n_inner);
      for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
	{
	  lohvs[unit_point].resize(n_outer, 0);
	  for (unsigned int k=0; k<n_outer; ++k)
	    lohvs[unit_point][k]=lohv_ptr[unit_point*n_outer+k];
	}
    }
  else
				   // not precomputed, then do so now
    compute_laplace_vector(lohvs);
    
				   // the sum of weights of the points
				   // at the outer rim should be
				   // one. check this
  for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
    Assert(fabs(std::accumulate(lohvs[unit_point].begin(),
				lohvs[unit_point].end(),0.) - 1)<1e-13,
	   ExcInternalError());
  
				   // TEST output
  if (false)
    {
      std::cout << "degree=" << degree << std::endl;
      for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
	for (unsigned int k=0; k<n_outer; ++k)
	  std::cout << lohvs[unit_point][k] << std::endl;
    }
}

#endif


template <int dim>
void
MappingQ<dim>::set_laplace_on_hex_vector(std::vector<std::vector<double> > &) const
{
  Assert(false, ExcInternalError());
}




#if deal_II_dimension==1

template <>
void
MappingQ<1>::compute_laplace_vector(std::vector<std::vector<double> > &) const
{
  Assert(false, ExcInternalError());
}

#else


template <int dim>
void
MappingQ<dim>::compute_laplace_vector(std::vector<std::vector<double> > &lvs) const
{
  Assert(lvs.size()==0, ExcInternalError());
  Assert(dim==2 || dim==3, ExcNotImplemented());
  Assert(degree>1, ExcInternalError());

				   // compute the shape
				   // gradients at the quadrature
				   // points on the unit cell
  const QGauss4<dim> quadrature;
  const unsigned int n_q_points=quadrature.n_quadrature_points;
  
  InternalData quadrature_data(n_shape_functions);
  quadrature_data.shape_derivatives.resize(n_shape_functions * n_q_points);
  compute_shapes(quadrature.get_points(), quadrature_data);
  
				   // Compute the stiffness matrix of
				   // the inner dofs
  FullMatrix<double> S(n_inner);
  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<n_inner; ++i)
      for (unsigned int j=0; j<n_inner; ++j)
	S(i,j)+=contract(quadrature_data.derivative(point, n_outer+i),
			 quadrature_data.derivative(point, n_outer+j))
		*quadrature.weight(point);
  
				   // Compute the components of T to be the
				   // product of gradients of inner and
				   // outer shape functions.
  FullMatrix<double> T(n_inner, n_outer);
  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<n_inner; ++i)
      for (unsigned int k=0; k<n_outer; ++k)
	T(i,k)+=contract(quadrature_data.derivative(point, n_outer+i),
			 quadrature_data.derivative(point, k))
		*quadrature.weight(point);
  
  FullMatrix<double> S_1(n_inner);
  S_1.invert(S);
  
  FullMatrix<double> S_1_T(n_inner, n_outer);
  
				   // S:=S_1*T
  S_1.mmult(S_1_T,T);
  
				   // compute the inner
				   // unit_support_points
  std::vector<Point<dim> > inner_unit_support_points(n_inner);
  const double step = 1./degree;
  const unsigned int z_end=(dim==3) ? degree : 2;
  unsigned int iall=0;
  for (unsigned int iz=1; iz<z_end; ++iz)
    for (unsigned int iy=1; iy<degree; ++iy)
      for (unsigned int ix=1; ix<degree; ++ix, ++iall)
	{
	  Point<dim> &p=inner_unit_support_points[iall];
	  p(0)=ix*step;
	  p(1)=iy*step;
	  if (dim==3)
	    p(2)=iz*step;
	}
  Assert(iall==n_inner, ExcInternalError());
  
				   // Compute the shape values at
				   // the inner
				   // unit_support_points
  InternalData support_data(n_shape_functions);
  support_data.shape_values.resize(n_shape_functions * n_inner);
  
  compute_shapes(inner_unit_support_points, support_data);
  
				   // Resize and initialize the
				   // lvs
  lvs.resize(n_inner);
  for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
    lvs[unit_point].resize(n_outer, 0);
  
				   // fill this vector
  for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
    {
      std::vector<double> &lv=lvs[unit_point];
      for (unsigned int k=0; k<n_outer; ++k)
	{
	  double sum=0;
	  for (unsigned int i=0; i<n_inner; ++i)
	    sum+=support_data.shape(unit_point, n_outer+i)
		 * S_1_T(i,k);
	      
	  lv[k]=-sum+support_data.shape(unit_point, k);
	}
    }
}

#endif



template <int dim>
void
MappingQ<dim>::apply_laplace_vector(const std::vector<std::vector<double> > &lvs,
				    std::vector<Point<dim> > &a) const
{
  Assert(lvs.size()!=0, ExcLaplaceVectorNotSet(degree));
  const unsigned int n_inner_apply=lvs.size();
  Assert(n_inner_apply==n_inner || n_inner_apply==(degree-1)*(degree-1),
	 ExcInternalError());
  const unsigned int n_outer_apply=lvs[0].size();
  Assert(a.size()==n_outer_apply, ExcDimensionMismatch(a.size(), n_outer_apply));

  for (unsigned int unit_point=0; unit_point<n_inner_apply; ++unit_point)
    {
      const std::vector<double> &lv=lvs[unit_point];
      Assert(lv.size()==n_outer_apply, ExcInternalError());
      Point<dim> p;
      for (unsigned int k=0; k<n_outer_apply; ++k)
	p+=lv[k]*a[k];

      a.push_back(p);
    }
}


template <int dim>
void
MappingQ<dim>::compute_mapping_support_points(
  const typename Triangulation<dim>::cell_iterator &cell,
  std::vector<Point<dim> > &a) const
{
  if (use_mapping_q_on_all_cells || cell->has_boundary_lines())
    compute_support_points_laplace(cell, a);
//  compute_support_points_simple(cell, a);
  else
    {
      a.resize(GeometryInfo<dim>::vertices_per_cell);
      
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	a[i] = cell->vertex(vertex_mapping[i]);
    }
}

  
template <int dim>
void
MappingQ<dim>::compute_support_points_laplace(const typename Triangulation<dim>::cell_iterator &cell,
					      std::vector<Point<dim> > &a) const
{
  a.resize(GeometryInfo<dim>::vertices_per_cell);
				   // the vertices first
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a[i] = cell->vertex(i);
  
  if (degree>1)
    {
      if (dim==1)
	{
	  Assert(false, ExcNotImplemented());
	}
      else
	{
					   // then the points on lines
					   // (for dim=2,3)
	  add_line_support_points (cell, a);
	  
	  if (dim==2)
	    apply_laplace_vector(laplace_on_quad_vector,a);
	  else if (dim==3)
	    {
	      add_face_support_points(cell, a);
	      
              apply_laplace_vector(laplace_on_hex_vector, a);
	    }
	}
    }    
}


template <int dim>
void
MappingQ<dim>::compute_support_points_simple(const typename Triangulation<dim>::cell_iterator &cell,
					     std::vector<Point<dim> > &a) const
{
				   // the vertices first
// TODO: check for size of 'a' first?
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a.push_back(cell->vertex(i));
  
  if (degree>1)
    {
				       // then the points on lines
				       // (for dim=2,3)
      add_line_support_points (cell, a);
      
				       // then the points on quads
				       // (for dim=3)
      fill_quad_support_points_simple (cell, a);
      
				       // then the points in cell. for
				       // this we need the midpoint of
				       // the points already in @p{a}
      const Point<dim> middle = std::accumulate(a.begin(), a.end(), Point<dim>())
				/ a.size();

      switch (degree)
	{
	  case 2:
	  {
	    a.push_back(middle);
	    break;
	  };

	  case 3:
	  {
					     // The four points in the
					     // cell are located at
					     // the midpoint between
					     // the middle point and
					     // the 4 vertices
	    
					     // TODO: better position of
					     // points: transform them by
					     // a Q2 transformation.
	    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	      a.push_back(middle*2./3.+cell->vertex(vertex_mapping[i])/3.);
	    break;
	  };

	  case 4:
	  {
	    Assert(a.size()==16, ExcInternalError());
	    a.insert(a.end(), 9, Point<dim>());
	    
	    const unsigned int inner_map[8]=
	    { 0, 1, 2, 5, 8, 7, 6, 3 };
	    
	    
					     // The nine points in the
					     // cell are located at the
					     // midpoint between the
					     // middle point and (the 4
					     // vertices and the face
					     // midpoints)
	    
	    a[16+4]=middle;
	    for (unsigned int i=0, j=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	      {
		a[16+inner_map[j++]]=(middle+cell->vertex(i))/2.;
		a[16+inner_map[j++]]=(middle+(cell->vertex(i)+cell->vertex((i+1)%4))/2.)/2.;
	      }
	    break;
	  };

	  default:
		Assert(false, ExcNotImplemented());
	};
    };
  
  Assert(a.size()==n_shape_functions, ExcInternalError());
}



#if deal_II_dimension==1

template <>
void
MappingQ<1>::add_line_support_points (const Triangulation<1>::cell_iterator &,
				      std::vector<Point<1> > &) const
{
				   // there are no points on bounding
				   // lines which are to be added
}

#endif


template <int dim>
void
MappingQ<dim>::add_line_support_points (const Triangulation<dim>::cell_iterator &cell,
					std::vector<Point<dim> > &a) const
{
  const Boundary<dim> *boundary = 0;

  std::vector<Point<dim> > line_points;
  if (degree>2)
    line_points.resize(degree-1);

				   // loop over each of the lines, and
				   // if it is at the boundary, then
				   // first get the boundary
				   // description and second compute
				   // the points on it
  for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
    {
      const typename Triangulation<dim>::line_iterator line = cell->line(line_no);
      if (line->at_boundary())
	boundary=&line->get_triangulation().get_boundary(line->boundary_indicator());
      else
	boundary=&straight_boundary;

				       // if we only need the
				       // midpoint, then ask for
				       // it. otherwise call the more
				       // complicated functions
      if (degree==2)
	a.push_back(boundary->get_new_point_on_line(line));
      else
	{
	  boundary->get_intermediate_points_on_line (line, line_points);
	  a.insert (a.end(), line_points.begin(), line_points.end());
	} 
    }
}




#if deal_II_dimension==3

//TODO: rename function to add_quad_support_points, to unify notation
template<>
void
MappingQ<3>::add_face_support_points(const Triangulation<3>::cell_iterator &cell,
				     std::vector<Point<3> > &a) const
{
  const unsigned int faces_per_cell    = GeometryInfo<3>::faces_per_cell,
		     vertices_per_face = GeometryInfo<3>::vertices_per_face,
		     lines_per_face    = GeometryInfo<3>::lines_per_face,
		     vertices_per_cell = GeometryInfo<3>::vertices_per_cell;
  
  static const unsigned int face_vertex_to_cell_vertex
    [faces_per_cell][vertices_per_face]={{0,1,2,3},
					 {4,5,6,7},
					 {0,1,5,4},
					 {1,5,6,2},
					 {3,2,6,7},
					 {0,4,7,3}};
  
  static const unsigned int face_line_to_cell_line
    [faces_per_cell][lines_per_face]={{0,1,2,3},
				      {4,5,6,7},
				      {0,9,4,8},
				      {9,5,10,1},
				      {2,10,6,11},
				      {8,7,11,3}};
  
  
				   // loop over all faces and collect points on them
  for (unsigned int face_no=0; face_no<faces_per_cell; ++face_no)
    {
      const Triangulation<3>::face_iterator face=cell->face(face_no);
      
      for (unsigned int i=0; i<vertices_per_face; ++i)
	Assert(face->vertex_index(i)==
	       cell->vertex_index(face_vertex_to_cell_vertex[face_no][i]),
	       ExcInternalError());
      
      for (unsigned int i=0; i<lines_per_face; ++i)
	Assert(face->line(i)==
	       cell->line(face_line_to_cell_line[face_no][i]),
	       ExcInternalError());

				       // if face at boundary, then
				       // ask boundary object to
				       // return intermediate points
				       // on it
      if (face->at_boundary())
	{
	  std::vector<Point<3> > quad_points ((degree-1)*(degree-1));

	  face->get_triangulation().get_boundary(face->boundary_indicator())
	    .get_intermediate_points_on_quad (face, quad_points);
	  a.insert (a.end(), quad_points.begin(), quad_points.end());
	}
      else
	{
					   // face is not at boundary,
					   // but maybe some of its
					   // lines are. count them
	  unsigned int lines_at_boundary=0;
	  for (unsigned int i=0; i<lines_per_face; ++i)
	    if (face->line(i)->at_boundary())
	      ++lines_at_boundary;
	  
	  Assert(lines_at_boundary<lines_per_face, ExcInternalError());

					   // if at least one of the
					   // lines bounding this quad
					   // is at the boundary, then
					   // collect points
					   // separately
	  if (lines_at_boundary>0)
	    {
					       // sort the points into b
//TODO: this is not thread-safe!!! b might be used for objects with
//TODO: different degrees at the same time!
	      static std::vector<Point<3> > b;
	      b.resize(4*degree);
	      Assert(4*degree==vertices_per_face+lines_per_face*(degree-1),
		     ExcDimensionMismatch(4*degree,
					  vertices_per_face+lines_per_face*(degree-1)));
	      for (unsigned int i=0; i<vertices_per_face; ++i)
		b[i]=a[face_vertex_to_cell_vertex[face_no][i]];
		      
	      for (unsigned int i=0; i<lines_per_face; ++i)
		for (unsigned int j=0; j<degree-1; ++j)
		  b[vertices_per_face+i*(degree-1)+j]=
		    a[vertices_per_cell+face_line_to_cell_line[face_no][i]*(degree-1)+j];
		  
	      apply_laplace_vector(laplace_on_quad_vector, b);
	      Assert(b.size()==4*degree+(degree-1)*(degree-1),
		     ExcDimensionMismatch(b.size(), 4*degree+(degree-1)*(degree-1)));
	      
	      for (unsigned int i=0; i<(degree-1)*(degree-1); ++i)
		a.push_back(b[4*degree+i]);
	    }
	  else
	    {
					       // face is entirely in
					       // the interior. get
					       // intermediate points
					       // from a straight
					       // boundary object
	      std::vector<Point<3> > quad_points ((degree-1)*(degree-1));
	      
	      straight_boundary.get_intermediate_points_on_quad (face, quad_points);
	      a.insert (a.end(), quad_points.begin(), quad_points.end());
	    }
	}
    }
}

#endif

template<int dim>
void
MappingQ<dim>::add_face_support_points(const typename Triangulation<dim>::cell_iterator &,
				       std::vector<Point<dim> > &) const
{
  Assert(false, ExcInternalError());
}




#if deal_II_dimension==3

template <>
void
MappingQ<3>::fill_quad_support_points_simple (const Triangulation<3>::cell_iterator &cell,
					      std::vector<Point<3> > &a) const
{
  const Boundary<3> *boundary = 0;

  std::vector<Point<3> > quad_points;
  Assert(degree>1, ExcInternalError());
  quad_points.resize((degree-1)*(degree-1));
  
  for (unsigned int quad_no=0; quad_no<GeometryInfo<3>::quads_per_cell; ++quad_no)
    {
      const Triangulation<3>::quad_iterator quad = cell->face(quad_no);
      if (quad->at_boundary())
	boundary=&quad->get_triangulation().get_boundary(quad->boundary_indicator());
      else
	boundary=&straight_boundary;

      boundary->get_intermediate_points_on_quad (quad, quad_points);
      a.insert (a.end(), quad_points.begin(), quad_points.end());
    }
}

#endif

template <int dim>
void
MappingQ<dim>::fill_quad_support_points_simple (const Triangulation<dim>::cell_iterator &,
						std::vector<Point<dim> > &) const
{}



//TODO: remove call of cross_product for dim==2
#if deal_II_dimension==2

void cross_product (Tensor<1,2> &, const Tensor<1,2> &, const Tensor<1,2> &)
{
  Assert(false, ExcInternalError());
}

#endif

template <int dim>
void
MappingQ<dim>::compute_fill_face (const typename DoFHandler<dim>::cell_iterator &cell,
				  const unsigned int            face_no,
				  const bool                    is_subface,
				  const unsigned int            npts,
				  const unsigned int            offset,
				  const std::vector<double>         &weights,
				  MappingQ1<dim>::InternalData &mapping_q1_data,
				  std::vector<Point<dim> >          &quadrature_points,
				  std::vector<double>               &JxW_values,
				  std::vector<Tensor<1,dim> >       &boundary_forms,
				  std::vector<Point<dim> >          &normal_vectors) const
{
  MappingQ1<dim>::compute_fill_face (cell, face_no, is_subface,
				     npts,
				     offset,
				     weights,
				     mapping_q1_data,
				     quadrature_points,
				     JxW_values,
				     boundary_forms,
				     normal_vectors);
  
  const UpdateFlags update_flags(mapping_q1_data.current_update_flags());

  if ((update_flags & update_normal_vectors)
      && alternative_normals_computation)
    {
      InternalData *data_ptr = dynamic_cast<InternalData *> (&mapping_q1_data);
      Assert(data_ptr!=0, ExcInternalError());
      InternalData &data=*data_ptr;

      transform_covariant(normal_vectors,
			  data.unit_normals[face_no],
			  data, 0);
      
      for (unsigned int i=0; i<normal_vectors.size(); ++i)
	normal_vectors[i] /= sqrt(normal_vectors[i].square());
    }
}



template <int dim>
void
MappingQ<dim>::transform_covariant (std::vector<Tensor<1,dim> >       &dst,
				    const std::vector<Tensor<1,dim> > &src,
				    const Mapping<dim>::InternalDataBase &mapping_data,
				    const unsigned int src_offset) const
{
  const MappingQ1<dim>::InternalData *q1_data_ptr =
    dynamic_cast<const MappingQ1<dim>::InternalData *> (&mapping_data);
  Assert(q1_data_ptr!=0, ExcInternalError());
  const MappingQ1<dim>::InternalData &q1_data=*q1_data_ptr;

  if (q1_data.is_mapping_q1_data)
    covariant_transformation(dst, src, q1_data, src_offset);
  else
    {
      const InternalData *data_ptr = dynamic_cast<const InternalData *> (q1_data_ptr);
      Assert(data_ptr!=0, ExcInternalError());
      const InternalData &data=*data_ptr;

      if (data.use_mapping_q1_on_current_cell)
	covariant_transformation(dst, src, data.mapping_q1_data, src_offset);
      else
	covariant_transformation(dst, src, data, src_offset);    
    }
}


template <int dim>
void
MappingQ<dim>::transform_covariant (std::vector<Point<dim> >       &dst,
				    const std::vector<Point<dim> > &src,
				    const Mapping<dim>::InternalDataBase &mapping_data,
				    const unsigned int src_offset) const
{
  const MappingQ1<dim>::InternalData *q1_data_ptr =
    dynamic_cast<const MappingQ1<dim>::InternalData *> (&mapping_data);
  Assert(q1_data_ptr!=0, ExcInternalError());
  const MappingQ1<dim>::InternalData &q1_data=*q1_data_ptr;

  if (q1_data.is_mapping_q1_data)
    covariant_transformation(dst, src, q1_data, src_offset);
  else
    {
      const InternalData *data_ptr = dynamic_cast<const InternalData *> (q1_data_ptr);
      Assert(data_ptr!=0, ExcInternalError());
      const InternalData &data=*data_ptr;

      if (data.use_mapping_q1_on_current_cell)
	covariant_transformation(dst, src, data.mapping_q1_data, src_offset);
      else
	covariant_transformation(dst, src, data, src_offset);    
    }  
}


template <int dim>
void
MappingQ<dim>::transform_contravariant (std::vector<Tensor<1,dim> >       &dst,
					const std::vector<Tensor<1,dim> > &src,
					const Mapping<dim>::InternalDataBase &mapping_data,
					const unsigned int src_offset) const
{
  const MappingQ1<dim>::InternalData *q1_data_ptr =
    dynamic_cast<const MappingQ1<dim>::InternalData *> (&mapping_data);
  Assert(q1_data_ptr!=0, ExcInternalError());
  const MappingQ1<dim>::InternalData &q1_data=*q1_data_ptr;

  if (q1_data.is_mapping_q1_data)
    contravariant_transformation(dst, src, q1_data, src_offset);
  else
    {
      const InternalData *data_ptr = dynamic_cast<const InternalData *> (q1_data_ptr);
      Assert(data_ptr!=0, ExcInternalError());
      const InternalData &data=*data_ptr;

      if (data.use_mapping_q1_on_current_cell)
	contravariant_transformation(dst, src, data.mapping_q1_data, src_offset);
      else
	contravariant_transformation(dst, src, data, src_offset);    
    }
}


template <int dim>
void
MappingQ<dim>::transform_contravariant (std::vector<Point<dim> >       &dst,
					const std::vector<Point<dim> > &src,
					const Mapping<dim>::InternalDataBase &mapping_data,
					const unsigned int src_offset) const
{
  const MappingQ1<dim>::InternalData *q1_data_ptr =
    dynamic_cast<const MappingQ1<dim>::InternalData *> (&mapping_data);
  Assert(q1_data_ptr!=0, ExcInternalError());
  const MappingQ1<dim>::InternalData &q1_data=*q1_data_ptr;

  if (q1_data.is_mapping_q1_data)
    contravariant_transformation(dst, src, q1_data, src_offset);
  else
    {
      const InternalData *data_ptr = dynamic_cast<const InternalData *> (q1_data_ptr);
      Assert(data_ptr!=0, ExcInternalError());
      const InternalData &data=*data_ptr;

      if (data.use_mapping_q1_on_current_cell)
	contravariant_transformation(dst, src, data.mapping_q1_data, src_offset);
      else
	contravariant_transformation(dst, src, data, src_offset);    
    }  
}



// explicit instantiation
template class MappingQ<deal_II_dimension>;
