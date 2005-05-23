//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <fe/mapping_q.h>
#include <fe/fe_q.h>
#include <base/polynomial.h>
#include <base/quadrature.h>
#include <base/quadrature_lib.h>
#include <base/memory_consumption.h>
#include <base/tensor_product_polynomials.h>
#include <lac/full_matrix.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_tools.h>

#include <numeric>
#include <memory>


template <int dim>
const bool MappingQ<dim>::use_mapping_q_on_all_cells;


template<int dim>
MappingQ<dim>::InternalData::InternalData (const unsigned int n_shape_functions)
		:
		MappingQ1<dim>::InternalData(n_shape_functions),
				use_mapping_q1_on_current_cell(false),
				mapping_q1_data(1 << dim)
{
  this->is_mapping_q1_data=false;
}



template<int dim>
unsigned int
MappingQ<dim>::InternalData::memory_consumption () const 
{
  return (MappingQ1<dim>::InternalData::memory_consumption () +
	  MemoryConsumption::memory_consumption (unit_normals) +
	  MemoryConsumption::memory_consumption (use_mapping_q1_on_current_cell) +
	  MemoryConsumption::memory_consumption (mapping_q1_data));
}



#if deal_II_dimension == 1

// in 1d, it is irrelevant which polynomial degree to use, since all
// cells are scaled linearly
template<>
MappingQ<1>::MappingQ (const unsigned int):
		degree(1),
		n_inner(0),
		n_outer(0),
		tensor_pols(0),
		n_shape_functions(2),
		renumber(0)
{}


template<>
MappingQ<1>::~MappingQ ()
{}

#endif



template<typename number>
number power(const number x, const unsigned int y)
{
				   // since the power to which x is
				   // raised is usually the space
				   // dimension, and since this is
				   // rarely larger than three, the
				   // following code is optimal and
				   // cannot be further optimized by
				   // grouping of operands to reduce
				   // the number of multiplications
				   // from O(x) to O(log x)
  number value=1;
  for (unsigned int i=0; i<y; ++i)
    value *= x;
  return value;
}


template<int dim>
MappingQ<dim>::MappingQ (const unsigned int p)
                :
		degree(p),
		n_inner(power(degree-1, dim)),
		n_outer((dim==2) ? 4+4*(degree-1)
			:8+12*(degree-1)+6*(degree-1)*(degree-1)),
		tensor_pols(0),
		n_shape_functions(power(degree+1,dim)),
		renumber(0)
{
				   // Construct the tensor product
				   // polynomials used as shape
				   // functions for the Qp mapping of
				   // cells at the boundary.
  std::vector<Polynomials::LagrangeEquidistant> v;
  for (unsigned int i=0; i<=degree; ++i)
    v.push_back(Polynomials::LagrangeEquidistant(degree,i));

  tensor_pols = new TensorProductPolynomials<dim> (v);
  Assert (n_shape_functions==tensor_pols->n(),
	  ExcInternalError());
  Assert(n_inner+n_outer==n_shape_functions, ExcInternalError());
  
				   // build the renumbering of the
				   // shape functions of the Qp
				   // mapping.
  renumber.resize(n_shape_functions,0);
  FETools::lexicographic_to_hierarchic_numbering (FE_Q<dim>(degree),
						  renumber);

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
				       typename MappingQ1<dim>::InternalData &data) const
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
typename Mapping<dim>::InternalDataBase *
MappingQ<dim>::get_data (const UpdateFlags update_flags,
			 const Quadrature<dim> &quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  this->compute_data (update_flags, quadrature,
                      quadrature.n_quadrature_points, *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_data (update_flags, quadrature,
                        quadrature.n_quadrature_points, data->mapping_q1_data);
  return data;
}



template <int dim>
typename Mapping<dim>::InternalDataBase *
MappingQ<dim>::get_face_data (const UpdateFlags update_flags,
			      const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_faces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.n_quadrature_points, *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_face_data (update_flags, q,
                             quadrature.n_quadrature_points,
                             data->mapping_q1_data);
  return data;
}



template <int dim>
typename Mapping<dim>::InternalDataBase *
MappingQ<dim>::get_subface_data (const UpdateFlags update_flags,
				 const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_subfaces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.n_quadrature_points, *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_face_data (update_flags, q,
                             quadrature.n_quadrature_points,
                             data->mapping_q1_data);
  return data;
}



template <int dim>
void
MappingQ<dim>::fill_fe_values (const typename Triangulation<dim>::cell_iterator &cell,
			       const Quadrature<dim>                &q,
			       typename Mapping<dim>::InternalDataBase       &mapping_data,
			       std::vector<Point<dim> >             &quadrature_points,
			       std::vector<double>                  &JxW_values) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &data = dynamic_cast<InternalData&> (mapping_data);

				   // check whether this cell needs
				   // the full mapping or can be
				   // treated by a reduced Q1 mapping,
				   // e.g. if the cell is in the
				   // interior of the domain
  data.use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
					  || cell->has_boundary_lines());

				   // depending on this result, use
				   // this or the other data object
				   // for the mapping
  typename MappingQ1<dim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    p_data=&data;
  
  MappingQ1<dim>::fill_fe_values(cell, q, *p_data,
				 quadrature_points, JxW_values);
}



template <int dim>
void
MappingQ<dim>::fill_fe_face_values (const typename Triangulation<dim>::cell_iterator &cell,
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
  
				   // check whether this cell needs
				   // the full mapping or can be
				   // treated by a reduced Q1 mapping,
				   // e.g. if the cell is entirely in
				   // the interior of the domain. note
				   // that it is not sufficient to ask
				   // whether the present _face_ is in
				   // the interior, as the mapping on
				   // the face depends on the mapping
				   // of the cell, which in turn
				   // depends on the fact whether
				   // _any_ of the faces of this cell
				   // is at the boundary, not only the
				   // present face
  data.use_mapping_q1_on_current_cell=!(use_mapping_q_on_all_cells
					|| cell->has_boundary_lines());

				   // depending on this result, use
				   // this or the other data object
				   // for the mapping
  typename MappingQ1<dim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    p_data=&data;

  const unsigned int n_q_points=q.n_quadrature_points;
  this->compute_fill_face (cell, face_no, false,
                           n_q_points,
                           QProjector<dim>::DataSetDescriptor::
                           face (face_no, cell->face_orientation(face_no),
                                 n_q_points),
                           q.get_weights(),
                           *p_data,
                           quadrature_points, JxW_values,
                           exterior_forms, normal_vectors);
}


template <int dim>
void
MappingQ<dim>::fill_fe_subface_values (const typename Triangulation<dim>::cell_iterator &cell,
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

				   // check whether this cell needs
				   // the full mapping or can be
				   // treated by a reduced Q1 mapping,
				   // e.g. if the cell is entirely in
				   // the interior of the domain. note
				   // that it is not sufficient to ask
				   // whether the present _face_ is in
				   // the interior, as the mapping on
				   // the face depends on the mapping
				   // of the cell, which in turn
				   // depends on the fact whether
				   // _any_ of the faces of this cell
				   // is at the boundary, not only the
				   // present face
  data.use_mapping_q1_on_current_cell=!(use_mapping_q_on_all_cells
					|| cell->has_boundary_lines());

				   // depending on this result, use
				   // this or the other data object
				   // for the mapping
  typename MappingQ1<dim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    p_data=&data;

  const unsigned int n_q_points=q.n_quadrature_points;
  this->compute_fill_face (cell, face_no, true,
                           n_q_points,
                           QProjector<dim>::DataSetDescriptor::
                           sub_face (face_no, sub_no,
                                     cell->face_orientation(face_no),
                                     n_q_points),
                           q.get_weights(),
                           *p_data,
                           quadrature_points, JxW_values,
                           exterior_forms, normal_vectors);
}


#if deal_II_dimension==1

template <>
void
MappingQ<1>::set_laplace_on_quad_vector(Table<2,double> &) const
{
  Assert(false, ExcInternalError());
}

#else

template <int dim>
void
MappingQ<dim>::set_laplace_on_quad_vector(Table<2,double> &loqvs) const
{
  Assert(degree>1, ExcInternalError());
  const unsigned int n_inner_2d=(degree-1)*(degree-1);
  const unsigned int n_outer_2d=4+4*(degree-1);

				   // first check whether we have precomputed
				   // the values for some polynomial degree;
				   // the sizes of arrays is
				   // n_inner_2d*n_outer_2d
  double const *loqv_ptr=0;
  switch (degree)
    {
                                       // for degree==1, we shouldn't have to
                                       // compute any support points, since
                                       // all of them are on the vertices
      
      case 2:
      {
                                         // (checked these values against the
                                         // output of compute_laplace_vector
                                         // again, and found they're indeed
                                         // right -- just in case someone
                                         // wonders where they come from --
                                         // WB)
	static const double loqv2[1*8]
	  ={1/16., 1/16., 1/16., 1/16., 3/16., 3/16., 3/16., 3/16.};
	loqv_ptr=&loqv2[0];
        Assert (sizeof(loqv2)/sizeof(loqv2[0]) ==
                n_inner_2d * n_outer_2d,
                ExcInternalError());
	
	break;
      }

      case 3:
      {
                                         // (same as above)
	static const double loqv3[4*12]
	  ={80/1053., 1/81., 11/1053., 1/81., 25/117., 44/351.,
	    7/117., 16/351., 7/117., 16/351., 25/117., 44/351.,
	    1/81., 80/1053., 1/81., 11/1053., 44/351., 25/117.,
	    25/117., 44/351., 16/351., 7/117., 7/117., 16/351.,
	    1/81., 11/1053., 1/81., 80/1053., 7/117., 16/351.,
	    16/351., 7/117., 25/117., 44/351., 44/351., 25/117.,
	    11/1053., 1/81., 80/1053., 1/81., 16/351., 7/117.,
	    44/351., 25/117., 44/351., 25/117., 16/351., 7/117.};
        Assert (sizeof(loqv3)/sizeof(loqv3[0]) ==
                n_inner_2d * n_outer_2d,
                ExcInternalError());
        
	loqv_ptr=&loqv3[0];
	
	break;
      }


      case 4:
      {
	static const double loqv4[9*16]
	  ={0.074059218503115712, -0.0010757446289059922,
            0.0019142922390714627, -0.0010757446289060069,
            0.22312738654318912, 0.13468513060151868,
            0.03812914216116723, 0.029131600026332517,
            0.02200737428129396, 0.016008355644312244,
            0.029131600026332527, 0.022007374281293915,
            0.016008355644312224, 0.22312738654318917,
            0.13468513060151874, 0.038129142161167237,
            
            0.0066480315133420603, 0.0066480315133427186,
            0.0028734528616576258, 0.0028734528616574584,
            0.15277168188202372, 0.23481527607092728,
            0.15277168188202397, 0.079035726825841868,
            0.059692382812499903, 0.036198648174158153,
            0.024962693117977778, 0.040819489554071289,
            0.024962693117977889, 0.079035726825843783,
            0.059692382812500312, 0.036198648174158236,
            
            -0.0010757446289069229, 0.074059218503115892,
            -0.0010757446289058842, 0.0019142922390713395,
            0.038129142161167293, 0.13468513060151846,
            0.22312738654318981, 0.22312738654318776,
            0.13468513060151835, 0.038129142161167202,
            0.016008355644312171, 0.022007374281293943,
            0.029131600026332617, 0.029131600026335094,
            0.02200737428129395, 0.016008355644312293,
            
            0.0066480315133420733, 0.0028734528616574727,
            0.0028734528616576362, 0.0066480315133427255,
            0.079035726825843755, 0.059692382812500257,
            0.036198648174158243, 0.024962693117977792,
            0.04081948955407131, 0.024962693117977903,
            0.079035726825841868, 0.059692382812499799,
            0.036198648174158098, 0.15277168188202384,
            0.23481527607092734, 0.15277168188202397,
            
            0.011067708333333018, 0.011067708333333358,
            0.011067708333333736, 0.01106770833333337,
            0.067708333333334217, 0.1035156250000009,
            0.067708333333334356, 0.067708333333333759,
            0.10351562499999903, 0.067708333333333995,
            0.067708333333333703, 0.10351562499999885,
            0.067708333333333898, 0.067708333333334245,
            0.10351562500000108, 0.067708333333334397,
            
            0.0028734528616571847, 0.0066480315133423621,
            0.0066480315133430378, 0.0028734528616573343,
            0.036198648174158195, 0.059692382812500278,
            0.079035726825844074, 0.15277168188202336,
            0.23481527607092598, 0.15277168188202367,
            0.036198648174158042, 0.059692382812499861,
            0.079035726825842201, 0.024962693117977788,
            0.040819489554074009, 0.024962693117977879,
            
            -0.0010757446289069133, 0.0019142922390713401,
            -0.0010757446289058651, 0.07405921850311592,
            0.029131600026335087, 0.022007374281293911,
            0.016008355644312279, 0.016008355644312175,
            0.022007374281293957, 0.029131600026332641,
            0.22312738654318776, 0.13468513060151827,
            0.038129142161167182, 0.038129142161167286,
            0.13468513060151852, 0.22312738654318992,
            
            0.0028734528616571756, 0.0028734528616573209,
            0.0066480315133430369, 0.0066480315133423742,
            0.024962693117977757, 0.040819489554073919,
            0.024962693117977847, 0.036198648174158049,
            0.059692382812499924, 0.079035726825842215,
            0.15277168188202328, 0.23481527607092575,
            0.15277168188202356, 0.036198648174158167,
            0.059692382812500319, 0.079035726825844088,
            
            0.0019142922390712369, -0.0010757446289068027,
            0.07405921850311617, -0.0010757446289067784,
            0.016008355644312283, 0.022007374281293967,
            0.02913160002633523, 0.038129142161167258,
            0.13468513060151821, 0.22312738654318864,
            0.038129142161167265, 0.13468513060151813,
            0.22312738654318859, 0.016008355644312276,
            0.022007374281294009, 0.029131600026335244
          };
        
        Assert (sizeof(loqv4)/sizeof(loqv4[0]) ==
                n_inner_2d * n_outer_2d,
                ExcInternalError());
        
	loqv_ptr=&loqv4[0];
	
	break;
      }
      
				       // no other cases implemented,
				       // so simply fall through
      default:
            break;
    }
  
  if (loqv_ptr!=0)
    {
				       // precomputed. copy values to
				       // the loqvs array
      loqvs.reinit(n_inner_2d, n_outer_2d);
      for (unsigned int unit_point=0; unit_point<n_inner_2d; ++unit_point)
	for (unsigned int k=0; k<n_outer_2d; ++k)
	  loqvs[unit_point][k]=loqv_ptr[unit_point*n_outer_2d+k];
    }
  else
    {
				       // not precomputed, then do so now
      if (dim==2)
	compute_laplace_vector(loqvs);
      else
					 // computing the Laplace vector for
					 // faces is not supported in 3d at
					 // present. presumably, doing so
					 // would not be so hard: we would
					 // only have to call the function in
					 // 2d, i.e. the quad(=face) values in
					 // 3d are equal to the quad(=cell)
					 // values in 2d. however, that would
					 // require us to link in the 2d
					 // library, which is kind of awkward
					 // (note that compute_laplace_vector
					 // really makes use of a lot of 2d
					 // stuff, such as FEValues etc). an
					 // alternative would be to precompute
					 // the values of this array for a
					 // couple of higher mapping orders,
					 // pin down their values and insert
					 // them into the array above.
	Assert (false, ExcNotImplemented());
  }

				   // the sum of weights of the points
				   // at the outer rim should be
				   // one. check this
  for (unsigned int unit_point=0; unit_point<loqvs.n_rows(); ++unit_point)
    Assert(std::fabs(std::accumulate(loqvs[unit_point].begin(),
				     loqvs[unit_point].end(),0.)-1)<1e-13,
	   ExcInternalError());
}

#endif


#if deal_II_dimension==3

template <>
void
MappingQ<3>::set_laplace_on_hex_vector(Table<2,double> &lohvs) const
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
      lohvs.reinit(n_inner, n_outer);
      for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
	for (unsigned int k=0; k<n_outer; ++k)
	  lohvs[unit_point][k]=lohv_ptr[unit_point*n_outer+k];
    }
  else
				     // not precomputed, then do so now
    compute_laplace_vector(lohvs);
    
				   // the sum of weights of the points
				   // at the outer rim should be
				   // one. check this
  for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
    Assert(std::fabs(std::accumulate(lohvs[unit_point].begin(),
				     lohvs[unit_point].end(),0.) - 1)<1e-13,
	   ExcInternalError());
}

#endif


template <int dim>
void
MappingQ<dim>::set_laplace_on_hex_vector(Table<2,double> &) const
{
  Assert(false, ExcInternalError());
}




#if deal_II_dimension==1

template <>
void
MappingQ<1>::compute_laplace_vector(Table<2,double> &) const
{
  Assert(false, ExcInternalError());
}

#else


template <int dim>
void
MappingQ<dim>::compute_laplace_vector(Table<2,double> &lvs) const
{
  Assert(lvs.n_rows()==0, ExcInternalError());
  Assert(dim==2 || dim==3, ExcNotImplemented());

                                   // for degree==1, we shouldn't have to
                                   // compute any support points, since all of
                                   // them are on the vertices
  Assert(degree>1, ExcInternalError());

				   // compute the shape
				   // gradients at the quadrature
				   // points on the unit cell
  const QGauss<dim> quadrature(degree+1);
  const unsigned int n_q_points=quadrature.n_quadrature_points;
  
  InternalData quadrature_data(n_shape_functions);
  quadrature_data.shape_derivatives.resize(n_shape_functions * n_q_points);
  this->compute_shapes(quadrature.get_points(), quadrature_data);
  
				   // Compute the stiffness matrix of
				   // the inner dofs
  FullMatrix<double> S(n_inner);
  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<n_inner; ++i)
      for (unsigned int j=0; j<n_inner; ++j)
	S(i,j) += contract(quadrature_data.derivative(point, n_outer+i),
			   quadrature_data.derivative(point, n_outer+j))
		  * quadrature.weight(point);
  
				   // Compute the components of T to be the
				   // product of gradients of inner and
				   // outer shape functions.
  FullMatrix<double> T(n_inner, n_outer);
  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<n_inner; ++i)
      for (unsigned int k=0; k<n_outer; ++k)
	T(i,k) += contract(quadrature_data.derivative(point, n_outer+i),
			   quadrature_data.derivative(point, k))
		  *quadrature.weight(point);
  
  FullMatrix<double> S_1(n_inner);
  S_1.invert(S);
  
  FullMatrix<double> S_1_T(n_inner, n_outer);
  
				   // S:=S_1*T
  S_1.mmult(S_1_T,T);
  
				   // Resize and initialize the
				   // lvs
  lvs.reinit (n_inner, n_outer);
  for (unsigned int i=0; i<n_inner; ++i)
    for (unsigned int k=0; k<n_outer; ++k)
      lvs(i,k) = -S_1_T(i,k);
}

#endif



template <int dim>
void
MappingQ<dim>::apply_laplace_vector(const Table<2,double> &lvs,
				    std::vector<Point<dim> > &a) const
{
                                   // check whether the data we need
                                   // is really available. if you fail
                                   // here and if
                                   // lvs==laplace_on_quad_vector in
                                   // the calling function, then we
                                   // didn't compute the quad laplace
                                   // vector. this is mentioned in the
                                   // constructor of this class,
                                   // although I don't understand the
                                   // reason for not aborting there
                                   // any more [WB]
  Assert(lvs.n_rows()!=0, ExcLaplaceVectorNotSet(degree));
  
  const unsigned int n_inner_apply=lvs.n_rows();
  Assert(n_inner_apply==n_inner || n_inner_apply==(degree-1)*(degree-1),
	 ExcInternalError());
  const unsigned int n_outer_apply=lvs.n_cols();
  Assert(a.size()==n_outer_apply,
	 ExcDimensionMismatch(a.size(), n_outer_apply));

				   // compute each inner point as
				   // linear combination of the outer
				   // points. the weights are given by
				   // the lvs entries, the outer
				   // points are the first (existing)
				   // elements of a
  for (unsigned int unit_point=0; unit_point<n_inner_apply; ++unit_point)
    {
      Assert(lvs.n_cols()==n_outer_apply, ExcInternalError());
      Point<dim> p;
      for (unsigned int k=0; k<n_outer_apply; ++k)
	p+=lvs[unit_point][k]*a[k];

      a.push_back(p);
    }
}


template <int dim>
void
MappingQ<dim>::compute_mapping_support_points(
  const typename Triangulation<dim>::cell_iterator &cell,
  std::vector<Point<dim> > &a) const
{
				   // if this is a cell for which we
				   // want to compute the full
				   // mapping, then get them from the
				   // following function
  if (use_mapping_q_on_all_cells || cell->has_boundary_lines())
    compute_support_points_laplace(cell, a);
  else
				     // otherwise: use a Q1 mapping
				     // for which the mapping shape
				     // function support points are
				     // simply the vertices of the
				     // cell
    {
      a.resize(GeometryInfo<dim>::vertices_per_cell);
      
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	a[i] = cell->vertex(this->vertex_mapping[i]);
    }
}

  
template <int dim>
void
MappingQ<dim>::compute_support_points_laplace(const typename Triangulation<dim>::cell_iterator &cell,
					      std::vector<Point<dim> > &a) const
{
				   // in any case, we need the
				   // vertices first
  a.resize(GeometryInfo<dim>::vertices_per_cell);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a[i] = cell->vertex(i);
  
  if (degree>1)
    switch (dim)
      {
	case 2:
					   // in 2d, add the
					   // points on the four
					   // bounding lines to
					   // the exterior (outer)
					   // points
	  add_line_support_points (cell, a);
	  apply_laplace_vector (laplace_on_quad_vector,a);
	  break;

	case 3:
					   // in 3d also add the
					   // points located on
					   // the boundary faces
	  add_line_support_points (cell, a);
	  add_quad_support_points (cell, a);
	  apply_laplace_vector (laplace_on_hex_vector, a);
	  break;
	      
	default:
	  Assert(false, ExcNotImplemented());
	  break;
      };
}





#if deal_II_dimension==1

template <>
void
MappingQ<1>::add_line_support_points (const Triangulation<1>::cell_iterator &,
				      std::vector<Point<1> > &) const
{
				   // there are no points on bounding
				   // lines which are to be added
  const unsigned int dim=1;
  Assert (dim > 1, ExcImpossibleInDim(dim));
}

#endif


template <int dim>
void
MappingQ<dim>::add_line_support_points (const typename Triangulation<dim>::cell_iterator &cell,
					std::vector<Point<dim> > &a) const
{
  static const StraightBoundary<dim> straight_boundary;
				   // if we only need the midpoint,
				   // then ask for it.
  if (degree==2)
    {
      for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
	{
	  const typename Triangulation<dim>::line_iterator line = cell->line(line_no);
	  const Boundary<dim> * const boundary
	    = (line->at_boundary() ?
	       &line->get_triangulation().get_boundary(line->boundary_indicator()) :
	       &straight_boundary);
	  
	  a.push_back(boundary->get_new_point_on_line(line));
	};
    }
  else
				     // otherwise call the more
				     // complicated functions and ask
				     // for inner points from the
				     // boundary description
    {
      std::vector<Point<dim> > line_points (degree-1);
      
				       // loop over each of the lines,
				       // and if it is at the
				       // boundary, then first get the
				       // boundary description and
				       // second compute the points on
				       // it
      for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
	{
	  const typename Triangulation<dim>::line_iterator line = cell->line(line_no);
	  
	  const Boundary<dim> * const boundary
	    = (line->at_boundary() ?
	       &line->get_triangulation().get_boundary(line->boundary_indicator()) :
	       &straight_boundary);
	  
	  boundary->get_intermediate_points_on_line (line, line_points);
	  a.insert (a.end(), line_points.begin(), line_points.end());
	}
    }
}




#if deal_II_dimension==3


template<>
void
MappingQ<3>::
add_quad_support_points(const Triangulation<3>::cell_iterator &cell,
                        std::vector<Point<3> >                &a) const
{
  const unsigned int faces_per_cell    = GeometryInfo<3>::faces_per_cell,
		     vertices_per_face = GeometryInfo<3>::vertices_per_face,
		     lines_per_face    = GeometryInfo<3>::lines_per_face,
		     vertices_per_cell = GeometryInfo<3>::vertices_per_cell;

                                   // mapping of faces to vertex
                                   // indices for properly oriented
                                   // faces is given by
                                   // GeometryInfo<3>::face_to_cell_vertices
				   //
                                   // same if the face is in opposite
                                   // orientation
  static const unsigned int face_to_cell_vertices2
    [faces_per_cell][vertices_per_face]={{0,3,2,1},
					 {4,7,6,5},
					 {0,4,5,1},
					 {1,2,6,5},
					 {3,7,6,2},
					 {0,3,7,4}};

  static const StraightBoundary<3> straight_boundary;
				   // used if face quad at boundary or
				   // entirely in the interior of the
				   // domain
  std::vector<Point<3> > quad_points ((degree-1)*(degree-1));
				   // used if only one line of face
				   // quad is at boundary
  std::vector<Point<3> > b(4*degree);
  
  
				   // loop over all faces and collect
				   // points on them
  for (unsigned int face_no=0; face_no<faces_per_cell; ++face_no)
    {
      const Triangulation<3>::face_iterator face = cell->face(face_no);

                                       // select the correct mappings
                                       // for the present face
      const bool face_orientation = cell->face_orientation(face_no);

                                       // some sanity checks up front
      for (unsigned int i=0; i<vertices_per_face; ++i)
        Assert(face->vertex_index(i)==
               cell->vertex_index(face_orientation ?
				  GeometryInfo<3>::face_to_cell_vertices(face_no, i) :
				  face_to_cell_vertices2[face_no][i]),
               ExcInternalError());
      

				       // indices of the lines that
				       // bound a properly oriented
				       // face are given by
				       // GeometryInfo<3>::
				       // face_to_cell_lines(face_no, line_no)
				       //
				       // same for faces in wrong
				       // orientation is given by
				       // GeometryInfo<3>::
				       // face_to_cell_lines(face_no, 3-line_no)
      for (unsigned int i=0; i<lines_per_face; ++i)
        Assert(face->line(i)==
               cell->line(face_orientation ?
			  GeometryInfo<3>::face_to_cell_lines(face_no, i) :
			  GeometryInfo<3>::face_to_cell_lines(face_no, 3-i)),
               ExcInternalError());
      
				       // if face at boundary, then
				       // ask boundary object to
				       // return intermediate points
				       // on it
      if (face->at_boundary())
	{
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
	  
	  Assert(lines_at_boundary<=lines_per_face, ExcInternalError());

					   // if at least one of the
					   // lines bounding this quad
					   // is at the boundary, then
					   // collect points
					   // separately
	  if (lines_at_boundary>0)
	    {
					       // call of function
					       // apply_laplace_vector
					       // increases size of b
					       // about 1. There
					       // resize b for the
					       // case the mentioned
					       // function was already
					       // called.
	      b.resize(4*degree);
	      
					       // b is of size
					       // 4*degree, make sure
					       // that this is the
					       // right size
	      Assert(b.size()==vertices_per_face+lines_per_face*(degree-1),
		     ExcDimensionMismatch(b.size(),
                                          vertices_per_face+lines_per_face*(degree-1)));
	      
					       // sort the points into b
              for (unsigned int i=0; i<vertices_per_face; ++i)
                b[i]=a[face_orientation ?
		      GeometryInfo<3>::face_to_cell_vertices(face_no, i) :
		      face_to_cell_vertices2[face_no][i]];
		      
              for (unsigned int i=0; i<lines_per_face; ++i)
                for (unsigned int j=0; j<degree-1; ++j)
                  b[vertices_per_face+i*(degree-1)+j]=
                    a[vertices_per_cell+
		     (face_orientation ?
		      GeometryInfo<3>::face_to_cell_lines(face_no, i) :
		      GeometryInfo<3>::face_to_cell_lines(face_no, 3-i))*(degree-1)+j];

					       // Now b includes the
					       // right order of
					       // support points on
					       // the quad to apply
					       // the laplace vector
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
	      straight_boundary.get_intermediate_points_on_quad (face, quad_points);
	      a.insert (a.end(), quad_points.begin(), quad_points.end());
	    }
	}
    }
}

#endif


template<int dim>
void
MappingQ<dim>::
add_quad_support_points(const typename Triangulation<dim>::cell_iterator &,
                        std::vector<Point<dim> > &) const
{
  Assert (dim > 2, ExcImpossibleInDim(dim));
}



template <int dim>
void
MappingQ<dim>::transform_covariant (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<1,dim> > > output,
  const typename Mapping<dim>::InternalDataBase &mapping_data) const
{
  const typename MappingQ1<dim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());

  Assert (output.size() + offset <= input.size(), ExcInternalError());
  
  typename std::vector<Tensor<2,dim> >::const_iterator tensor;

  if (q1_data->is_mapping_q1_data)
    tensor = q1_data->covariant.begin();
  else
    {
      const InternalData *data = dynamic_cast<const InternalData *> (q1_data);
      Assert(data!=0, ExcInternalError());

      if (data->use_mapping_q1_on_current_cell)
	tensor = data->mapping_q1_data.covariant.begin();
      else
	tensor = data->covariant.begin();    
    }

  for (unsigned int i=0; i<output.size(); ++i)
    contract (output[i], input[i+offset], *(tensor++));
}



template <int dim>
void
MappingQ<dim>::transform_covariant (
  const VectorSlice<const std::vector<Tensor<2,dim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<2,dim> > > output,
  const typename Mapping<dim>::InternalDataBase &mapping_data) const
{
  const typename MappingQ1<dim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());

  Assert (output.size() + offset <= input.size(), ExcInternalError());

  typename std::vector<Tensor<2,dim> >::const_iterator tensor;

  if (q1_data->is_mapping_q1_data)
    tensor = q1_data->covariant.begin();
  else
    {
      const InternalData *data = dynamic_cast<const InternalData *> (q1_data);
      Assert(data!=0, ExcInternalError());

      if (data->use_mapping_q1_on_current_cell)
	tensor = data->mapping_q1_data.covariant.begin();
      else
	tensor = data->covariant.begin();
    }

  for (unsigned int i=0; i<output.size(); ++i)
    contract (output[i], input[i+offset], *(tensor++));
}



template <int dim>
void
MappingQ<dim>::transform_contravariant (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<1,dim> > > output,
  const typename Mapping<dim>::InternalDataBase &mapping_data) const
{
  const typename MappingQ1<dim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());
  
  Assert (output.size() + offset <= input.size(), ExcInternalError());

  typename std::vector<Tensor<2,dim> >::const_iterator tensor;

  if (q1_data->is_mapping_q1_data)
    tensor = q1_data->contravariant.begin();
  else
    {
      const InternalData *data = dynamic_cast<const InternalData *> (q1_data);
      Assert(data!=0, ExcInternalError());

      if (data->use_mapping_q1_on_current_cell)
	tensor = data->mapping_q1_data.contravariant.begin();
      else
	tensor = data->contravariant.begin();    
    }
  
  for (unsigned int i=0; i<output.size(); ++i)
    contract (output[i], *(tensor++), input[i+offset]);
}



template <int dim>
void
MappingQ<dim>::transform_contravariant (
  const VectorSlice<const std::vector<Tensor<2,dim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<2,dim> > > output,
  const typename Mapping<dim>::InternalDataBase &mapping_data) const
{
  const typename MappingQ1<dim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());
  
  Assert (output.size() + offset <= input.size(), ExcInternalError());

  typename std::vector<Tensor<2,dim> >::const_iterator tensor;

  if (q1_data->is_mapping_q1_data)
    tensor = q1_data->contravariant.begin();
  else
    {
      const InternalData *data = dynamic_cast<const InternalData *> (q1_data);
      Assert(data!=0, ExcInternalError());

      if (data->use_mapping_q1_on_current_cell)
	tensor = data->mapping_q1_data.contravariant.begin();
      else
	tensor = data->contravariant.begin();    
    }

  for (unsigned int i=0; i<output.size(); ++i)
    contract (output[i], *(tensor++), input[i+offset]);
}



template <int dim>
Point<dim>
MappingQ<dim>::
transform_unit_to_real_cell (const typename Triangulation<dim>::cell_iterator &cell,
                             const Point<dim>                                 &p) const
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
  
  mdata->use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
					    || cell->has_boundary_lines());

  typename MappingQ1<dim>::InternalData
    *p_data = (mdata->use_mapping_q1_on_current_cell ?
               &mdata->mapping_q1_data :
               &*mdata);

  compute_mapping_support_points(cell, p_data->mapping_support_points);
  
  return this->transform_unit_to_real_cell_internal(*p_data);
}



template <int dim>
Point<dim>
MappingQ<dim>::
transform_real_to_unit_cell (const typename Triangulation<dim>::cell_iterator &cell,
                             const Point<dim>                                 &p) const
{
				   // first a Newton iteration based
				   // on a Q1 mapping
  Point<dim> p_unit = MappingQ1<dim>::transform_real_to_unit_cell(cell, p);
  
                                   // then a Newton iteration based on
                                   // the full MappingQ if we need
                                   // this
  if (cell->has_boundary_lines() || use_mapping_q_on_all_cells)
    {
      const Quadrature<dim> point_quadrature(p_unit);
      std::auto_ptr<InternalData>
        mdata (dynamic_cast<InternalData *> (
                 get_data(update_transformation_values |
                          update_transformation_gradients,
                          point_quadrature)));
      
      mdata->use_mapping_q1_on_current_cell = false;

      std::vector<Point<dim> > &points = mdata->mapping_support_points;
      compute_mapping_support_points (cell, points);

      this->transform_real_to_unit_cell_internal(cell, p, *mdata, p_unit);
    }
  
  return p_unit;
}



template <int dim>
unsigned int
MappingQ<dim>::get_degree() const
{
  return degree;
}

  
// explicit instantiation
template class MappingQ<deal_II_dimension>;
