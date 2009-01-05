//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/utilities.h>
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
#include <fe/mapping_q.h>
#include <fe/fe_q.h>

#include <numeric>
#include <memory>

DEAL_II_NAMESPACE_OPEN


template<int dim, int spacedim>
MappingQ<dim,spacedim>::InternalData::InternalData (const unsigned int n_shape_functions)
		:
		MappingQ1<dim,spacedim>::InternalData(n_shape_functions),
                use_mapping_q1_on_current_cell(false),
		mapping_q1_data(1 << dim)
{
  this->is_mapping_q1_data=false;
}



template<int dim, int spacedim>
unsigned int
MappingQ<dim,spacedim>::InternalData::memory_consumption () const 
{
  return (MappingQ1<dim,spacedim>::InternalData::memory_consumption () +
	  MemoryConsumption::memory_consumption (unit_normals) +
	  MemoryConsumption::memory_consumption (use_mapping_q1_on_current_cell) +
	  MemoryConsumption::memory_consumption (mapping_q1_data));
}



#if deal_II_dimension == 1

// in 1d, it is irrelevant which polynomial degree to use, since all
// cells are scaled linearly
template<>
MappingQ<1>::MappingQ (const unsigned int,
		       const bool /*use_mapping_q_on_all_cells*/)
		:
		degree(1),
		n_inner(0),
		n_outer(0),
		tensor_pols(0),
		n_shape_functions(2),
		renumber(0),
		use_mapping_q_on_all_cells (false),
		feq(degree)
{}


template<>
MappingQ<1>::MappingQ (const MappingQ<1> &m):
		MappingQ1<1> (),
		degree(1),
		n_inner(0),
		n_outer(0),
		tensor_pols(0),
		n_shape_functions(2),
		renumber(0),
		use_mapping_q_on_all_cells (m.use_mapping_q_on_all_cells),
		feq(degree)
{}


template<>
MappingQ<1>::~MappingQ ()
{}

#endif



namespace
{
  template <int dim>
  std::vector<unsigned int>
  get_dpo_vector (const unsigned int degree)
  {
    std::vector<unsigned int> dpo(dim+1, 1U);
    for (unsigned int i=1; i<dpo.size(); ++i)
      dpo[i]=dpo[i-1]*(degree-1);
    return dpo;
  }
}




template<int dim, int spacedim>
MappingQ<dim,spacedim>::MappingQ (const unsigned int p,
			 const bool use_mapping_q_on_all_cells)
                :
		degree(p),
		n_inner(Utilities::fixed_power<dim>(degree-1)),
		n_outer((dim==2) ? 4+4*(degree-1)
			:8+12*(degree-1)+6*(degree-1)*(degree-1)),
		tensor_pols(0),
		n_shape_functions(Utilities::fixed_power<dim>(degree+1)),
		renumber(FETools::
			 lexicographic_to_hierarchic_numbering (
			   FiniteElementData<dim> (get_dpo_vector<dim>(degree), 1,
						   degree))),
		use_mapping_q_on_all_cells (use_mapping_q_on_all_cells),
		feq(degree)
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
  
				   // build laplace_on_quad_vector
  if (degree>1)
    {
      if (dim >= 2)
	set_laplace_on_quad_vector(laplace_on_quad_vector);
      if (dim >= 3)
	set_laplace_on_hex_vector(laplace_on_hex_vector);
    }
}


template<int dim, int spacedim>
MappingQ<dim,spacedim>::MappingQ (const MappingQ<dim,spacedim> &mapping)
		:
		MappingQ1<dim,spacedim>(),
		degree(mapping.degree),
		n_inner(mapping.n_inner),
		n_outer(n_outer),
		tensor_pols(0),
		n_shape_functions(mapping.n_shape_functions),
		renumber(mapping.renumber),
		use_mapping_q_on_all_cells (mapping.use_mapping_q_on_all_cells),
		feq(degree)
{
  tensor_pols=new TensorProductPolynomials<dim> (*mapping.tensor_pols);
  laplace_on_quad_vector=mapping.laplace_on_quad_vector;
  laplace_on_hex_vector=mapping.laplace_on_hex_vector;
}


template<int dim, int spacedim>
MappingQ<dim,spacedim>::~MappingQ ()
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



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
				       typename MappingQ1<dim,spacedim>::InternalData &data) const
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
  
//				   // dummy variable of size 0
  std::vector<Tensor<2,dim> > grad2;
  if (data.shape_second_derivatives.size()!=0)
    {
      Assert(data.shape_second_derivatives.size()==n_shape_functions*n_points,
	     ExcInternalError());
      grad2.resize(n_shape_functions);
    }

  
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

      	if (data.shape_second_derivatives.size()!=0)
	  for (unsigned int i=0; i<n_shape_functions; ++i)
	    data.second_derivative(point,renumber[i]) = grad2[i];
      }
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ<dim,spacedim>::get_data (const UpdateFlags update_flags,
			 const Quadrature<dim> &quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  this->compute_data (update_flags, quadrature,
                      quadrature.size(), *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_data (update_flags, quadrature,
                        quadrature.size(), data->mapping_q1_data);
  return data;
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ<dim,spacedim>::get_face_data (const UpdateFlags update_flags,
			      const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_faces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.size(), *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_face_data (update_flags, q,
                             quadrature.size(),
                             data->mapping_q1_data);
  return data;
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ<dim,spacedim>::get_subface_data (const UpdateFlags update_flags,
				 const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_subfaces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.size(), *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_face_data (update_flags, q,
                             quadrature.size(),
                             data->mapping_q1_data);
  return data;
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::fill_fe_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const Quadrature<dim>                                     &q,
  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
  std::vector<Point<spacedim> >                             &quadrature_points,
  std::vector<double>                                       &JxW_values,
  std::vector<Tensor<2,spacedim> >                          &jacobians,
  std::vector<Tensor<3,spacedim> >                          &jacobian_grads,
  std::vector<Tensor<2,spacedim> >                          &inverse_jacobians,
  std::vector<Point<spacedim> >                             &cell_normal_vectors) const
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
  typename MappingQ1<dim,spacedim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    p_data=&data;
  
  MappingQ1<dim,spacedim>::fill_fe_values(cell, q, *p_data,
	        			  quadrature_points, JxW_values,
				          jacobians, jacobian_grads, inverse_jacobians,
					  cell_normal_vectors);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::fill_fe_face_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const unsigned int       face_no,
  const Quadrature<dim-1> &q,
  typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  std::vector<Point<dim> >     &quadrature_points,
  std::vector<double>          &JxW_values,
  std::vector<Tensor<1,dim> >  &exterior_forms,
  std::vector<Point<dim> >     &normal_vectors,
  std::vector<double>          &cell_JxW_values) const
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
  typename MappingQ1<dim,spacedim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    p_data=&data;

  const unsigned int n_q_points=q.size();
  this->compute_fill_face (cell, face_no, deal_II_numbers::invalid_unsigned_int,
                           n_q_points,
                           QProjector<dim>::DataSetDescriptor::
                           face (face_no,
				 cell->face_orientation(face_no),
				 cell->face_flip(face_no),
				 cell->face_rotation(face_no),
                                 n_q_points),
                           q.get_weights(),
                           *p_data,
                           quadrature_points, JxW_values,
                           exterior_forms, normal_vectors,
			   cell_JxW_values);

  // TODO: Verify implementation of cell_JxW_values
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
				       const unsigned int       face_no,
				       const unsigned int       sub_no,
				       const Quadrature<dim-1> &q,
				       typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
				       std::vector<Point<dim> >     &quadrature_points,
				       std::vector<double>          &JxW_values,
				       std::vector<Tensor<1,dim> >  &exterior_forms,
				       std::vector<Point<dim> >     &normal_vectors,
				       std::vector<double>          &cell_JxW_values) const
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
  typename MappingQ1<dim,spacedim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    p_data=&data;

  const unsigned int n_q_points=q.size();
  this->compute_fill_face (cell, face_no, sub_no,
                           n_q_points,
                           QProjector<dim>::DataSetDescriptor::
                           subface (face_no, sub_no,
                                     cell->face_orientation(face_no),
                                     cell->face_flip(face_no),
                                     cell->face_rotation(face_no),
                                     n_q_points),
                           q.get_weights(),
                           *p_data,
                           quadrature_points, JxW_values,
                           exterior_forms, normal_vectors,
			   cell_JxW_values);

  // TODO: Verify implementation of cell_JxW_values ...
}


#if deal_II_dimension==1

template <>
void
MappingQ<1>::set_laplace_on_quad_vector(Table<2,double> &) const
{
  Assert(false, ExcInternalError());
}

#else

template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::set_laplace_on_quad_vector(Table<2,double> &loqvs) const
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
	  ={80/1053., 1/81., 1/81., 11/1053., 25/117., 44/351.,
	    7/117., 16/351., 25/117., 44/351., 7/117., 16/351., 
	    1/81., 80/1053., 11/1053., 1/81., 7/117., 16/351.,
	    25/117., 44/351., 44/351., 25/117., 16/351., 7/117.,
	    1/81., 11/1053., 80/1053., 1/81., 44/351., 25/117.,
	    16/351., 7/117., 7/117., 16/351., 25/117., 44/351.,
	    11/1053., 1/81., 1/81., 80/1053., 16/351., 7/117.,
	    44/351., 25/117., 16/351., 7/117., 44/351., 25/117.};
        Assert (sizeof(loqv3)/sizeof(loqv3[0]) ==
                n_inner_2d * n_outer_2d,
                ExcInternalError());
        
	loqv_ptr=&loqv3[0];
	
	break;
      }


      case 4:
      {
	static const double loqv4[9*16]
	  ={0.07405921850311571, -0.001075744628905992,
	    -0.001075744628906007, 0.001914292239071463,
	    0.2231273865431892, 0.1346851306015187,
	    0.03812914216116724, 0.02913160002633252,
	    0.02200737428129396, 0.01600835564431224,
	    0.2231273865431891, 0.1346851306015187,
	    0.03812914216116723, 0.02913160002633253,
	    0.02200737428129391, 0.01600835564431222,
	    
	    0.00664803151334206, 0.006648031513342719,
	    0.002873452861657458, 0.002873452861657626,
	    0.07903572682584378, 0.05969238281250031,
	    0.03619864817415824, 0.07903572682584187,
	    0.0596923828124999, 0.03619864817415815,
	    0.1527716818820237, 0.2348152760709273,
	    0.152771681882024, 0.02496269311797778,
	    0.04081948955407129, 0.02496269311797789,

	    -0.001075744628906923, 0.07405921850311589,
	    0.001914292239071339, -0.001075744628905884,
	    0.02913160002633509, 0.02200737428129395,
	    0.01600835564431229, 0.2231273865431878,
	    0.1346851306015183, 0.0381291421611672,
	    0.03812914216116729, 0.1346851306015185,
	    0.2231273865431898, 0.01600835564431217,
	    0.02200737428129394, 0.02913160002633262,
	    
	    0.006648031513342073, 0.002873452861657473,
	    0.006648031513342726, 0.002873452861657636,
	    0.1527716818820238, 0.2348152760709273,
	    0.152771681882024, 0.02496269311797779,
	    0.04081948955407131, 0.0249626931179779,
	    0.07903572682584376, 0.05969238281250026,
	    0.03619864817415824, 0.07903572682584187,
	    0.0596923828124998, 0.0361986481741581,
	    
	    0.01106770833333302, 0.01106770833333336,
	    0.01106770833333337, 0.01106770833333374,
	    0.06770833333333424, 0.1035156250000011,
	    0.0677083333333344, 0.06770833333333376,
	    0.103515624999999, 0.06770833333333399,
	    0.06770833333333422, 0.1035156250000009,
	    0.06770833333333436, 0.0677083333333337,
	    0.1035156249999988, 0.0677083333333339,
	    
	    0.002873452861657185, 0.006648031513342362,
	    0.002873452861657334, 0.006648031513343038,
	    0.02496269311797779, 0.04081948955407401,
	    0.02496269311797788, 0.1527716818820234,
	    0.234815276070926, 0.1527716818820237,
	    0.03619864817415819, 0.05969238281250028,
	    0.07903572682584407, 0.03619864817415804,
	    0.05969238281249986, 0.0790357268258422,
	    
	    -0.001075744628906913, 0.00191429223907134,
	    0.07405921850311592, -0.001075744628905865,
	    0.03812914216116729, 0.1346851306015185,
	    0.2231273865431899, 0.01600835564431217,
	    0.02200737428129396, 0.02913160002633264,
	    0.02913160002633509, 0.02200737428129391,
	    0.01600835564431228, 0.2231273865431878,
	    0.1346851306015183, 0.03812914216116718,

	    0.002873452861657176, 0.002873452861657321,
	    0.006648031513342374, 0.006648031513343037,
	    0.03619864817415817, 0.05969238281250032,
	    0.07903572682584409, 0.03619864817415805,
	    0.05969238281249992, 0.07903572682584221,
	    0.02496269311797776, 0.04081948955407392,
	    0.02496269311797785, 0.1527716818820233,
	    0.2348152760709258, 0.1527716818820236,
	    
	    0.001914292239071237, -0.001075744628906803,
	    -0.001075744628906778, 0.07405921850311617,
	    0.01600835564431228, 0.02200737428129401,
	    0.02913160002633524, 0.03812914216116726,
	    0.1346851306015182, 0.2231273865431886,
	    0.01600835564431228, 0.02200737428129397,
	    0.02913160002633523, 0.03812914216116726,
	    0.1346851306015181, 0.2231273865431886,    
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


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::set_laplace_on_hex_vector(Table<2,double> &) const
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


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::compute_laplace_vector(Table<2,double> &lvs) const
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
  const unsigned int n_q_points=quadrature.size();
  
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



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::apply_laplace_vector(const Table<2,double> &lvs,
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


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
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
	a[i] = cell->vertex(i);
    }
}

  
template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::compute_support_points_laplace(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
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
	{
					   // in 3d also add the
					   // points located on
					   // the boundary faces
	  add_line_support_points (cell, a);
	  add_quad_support_points (cell, a);
	  apply_laplace_vector (laplace_on_hex_vector, a);
	  break;
	}
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


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::add_line_support_points (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
					std::vector<Point<dim> > &a) const
{
  static const StraightBoundary<dim> straight_boundary;
				   // if we only need the midpoint,
				   // then ask for it.
  if (degree==2)
    {
      for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
	{
	  const typename Triangulation<dim,spacedim>::line_iterator line = cell->line(line_no);
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
	  const typename Triangulation<dim,spacedim>::line_iterator line = cell->line(line_no);
	  
	  const Boundary<dim> * const boundary
	    = (line->at_boundary() ?
	       &line->get_triangulation().get_boundary(line->boundary_indicator()) :
	       &straight_boundary);
	  
	  boundary->get_intermediate_points_on_line (line, line_points);
	  if (dim==3)
	    {
					       // in 3D, lines might be in wrong
					       // orientation. if so, reverse
					       // the vector
	      if (cell->line_orientation(line_no))
		a.insert (a.end(), line_points.begin(), line_points.end());
 	      else
 		a.insert (a.end(), line_points.rbegin(), line_points.rend());
	    }
	  else
					     // in 2D, lines always have the
					     // correct orientation. simply
					     // append all points
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
      const bool face_orientation = cell->face_orientation(face_no),
		 face_flip        = cell->face_flip       (face_no),
		 face_rotation    = cell->face_rotation   (face_no);

#ifdef DEBUG      
                                       // some sanity checks up front
      for (unsigned int i=0; i<vertices_per_face; ++i)
        Assert(face->vertex_index(i)==cell->vertex_index(
	  GeometryInfo<3>::face_to_cell_vertices(face_no, i,
						 face_orientation,
						 face_flip,
						 face_rotation)),
	  ExcInternalError());

				       // indices of the lines that
				       // bound a face are given by
				       // GeometryInfo<3>::
				       // face_to_cell_lines
      for (unsigned int i=0; i<lines_per_face; ++i)
        Assert(face->line(i)==cell->line(GeometryInfo<3>::face_to_cell_lines(
	  face_no, i, face_orientation, face_flip, face_rotation)),
	       ExcInternalError());
#endif
      
				       // if face at boundary, then
				       // ask boundary object to
				       // return intermediate points
				       // on it
      if (face->at_boundary())
	{
	  face->get_triangulation().get_boundary(face->boundary_indicator())
	    .get_intermediate_points_on_quad (face, quad_points);
					   // in 3D, the orientation, flip and
					   // rotation of the face might not
					   // match what we expect here, namely
					   // the standard orientation. thus
					   // reorder points accordingly. since
					   // a Mapping uses the same shape
					   // function as an FEQ, we can ask a
					   // FEQ to do the reordering for us.
	  for (unsigned int i=0; i<quad_points.size(); ++i)
	    a.push_back(quad_points[feq.adjust_quad_dof_index_for_face_orientation(i,
										   face_orientation,
										   face_flip,
										   face_rotation)]);
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
	      
					       // sort the points into b. We
					       // used access from the cell (not
					       // from the face) to fill b, so
					       // we can assume a standard face
					       // orientation. Doing so, the
					       // calculated points will be in
					       // standard orientation as well.
              for (unsigned int i=0; i<vertices_per_face; ++i)
                b[i]=a[GeometryInfo<3>::face_to_cell_vertices(face_no, i)];
		      
              for (unsigned int i=0; i<lines_per_face; ++i)
                for (unsigned int j=0; j<degree-1; ++j)
                  b[vertices_per_face+i*(degree-1)+j]=
                    a[vertices_per_cell + GeometryInfo<3>::face_to_cell_lines(
		      face_no, i)*(degree-1)+j];

					       // Now b includes the support
					       // points on the quad and we can
					       // apply the laplace vector
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
					       // in 3D, the orientation, flip
					       // and rotation of the face might
					       // not match what we expect here,
					       // namely the standard
					       // orientation. thus reorder
					       // points accordingly. since a
					       // Mapping uses the same shape
					       // function as an FEQ, we can ask
					       // a FEQ to do the reordering for
					       // us.
	      for (unsigned int i=0; i<quad_points.size(); ++i)
		a.push_back(quad_points[feq.adjust_quad_dof_index_for_face_orientation(i,
										       face_orientation,
										       face_flip,
										       face_rotation)]);
	    }
	}
    }
}

#endif


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
add_quad_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &,
                        std::vector<Point<dim> > &) const
{
  Assert (dim > 2, ExcImpossibleInDim(dim));
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  VectorSlice<std::vector<Tensor<1,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &internal,
  const MappingType mapping_type) const
{
  switch (mapping_type)
    {
      case mapping_covariant:
	    transform_covariant(input, 0, output, internal);
	    return;
       case mapping_contravariant:
 	    transform_contravariant(input, 0, output, internal);
 	    return;
      default:
	    Assert(false, ExcNotImplemented());
    }
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<2,dim> > > input,
  VectorSlice<std::vector<Tensor<2,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &internal,
  const MappingType mapping_type) const
{
  switch (mapping_type)
    {
      case mapping_covariant:
	    transform_covariant(input, 0, output, internal);
	    return;
       case mapping_contravariant:
 	    transform_contravariant(input, 0, output, internal);
 	    return;
      default:
	    Assert(false, ExcNotImplemented());
    }
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::transform_covariant (
  const VectorSlice<const std::vector<Tensor<1,spacedim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<1,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data) const
{
  Assert (offset == 0, ExcInternalError());
  AssertDimension (input.size(), output.size());
  
  const typename MappingQ1<dim,spacedim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim,spacedim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());

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
  
  Tensor<1, spacedim> auxiliary;
  
  for (unsigned int i=0; i<output.size(); ++i)
    {
      for (unsigned int d=0;d<dim;++d)
	auxiliary[d] = input[i][d];
      contract (output[i], auxiliary, *(tensor++));
    }
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::transform_covariant (
  const VectorSlice<const std::vector<Tensor<2,spacedim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<2,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data) const
{
  Assert (offset == 0, ExcInternalError());
  AssertDimension (input.size(), output.size());
  
  const typename MappingQ1<dim,spacedim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim,spacedim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());

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



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::transform_contravariant (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<1,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data) const
{
  Assert (offset == 0, ExcInternalError());
  AssertDimension (input.size(), output.size());
  
  const typename MappingQ1<dim,spacedim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim,spacedim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());
  
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
  
  Tensor<1, spacedim> auxiliary;
  
  for (unsigned int i=0; i<output.size(); ++i)
    {
      for (unsigned int d=0;d<dim;++d)
	auxiliary[d] = input[i][d];
      contract (output[i], *(tensor++), auxiliary);
    }
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::transform_contravariant (
  const VectorSlice<const std::vector<Tensor<2,dim> > > input,
  const unsigned int                 offset,
  VectorSlice<std::vector<Tensor<2,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data) const
{
  Assert (offset == 0, ExcInternalError());
  AssertDimension (input.size(), output.size());
  
  const typename MappingQ1<dim,spacedim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim,spacedim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());
  
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



template<int dim, int spacedim>
Point<spacedim>
MappingQ<dim,spacedim>::
transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
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

  typename MappingQ1<dim,spacedim>::InternalData
    *p_data = (mdata->use_mapping_q1_on_current_cell ?
               &mdata->mapping_q1_data :
               &*mdata);

  compute_mapping_support_points(cell, p_data->mapping_support_points);
  
  return this->transform_unit_to_real_cell_internal(*p_data);
}



template<int dim, int spacedim>
Point<dim>
MappingQ<dim,spacedim>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<spacedim>                            &p) const
{
				   // first a Newton iteration based
				   // on a Q1 mapping
  Point<dim> p_unit = MappingQ1<dim,spacedim>::transform_real_to_unit_cell(cell, p);
  
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



template<int dim, int spacedim>
unsigned int
MappingQ<dim,spacedim>::get_degree() const
{
  return degree;
}



template<int dim, int spacedim>
Mapping<dim,spacedim> *
MappingQ<dim,spacedim>::clone () const
{
  return new MappingQ<dim,spacedim>(*this);
}

  
// explicit instantiation
template class MappingQ<deal_II_dimension>;

DEAL_II_NAMESPACE_CLOSE
