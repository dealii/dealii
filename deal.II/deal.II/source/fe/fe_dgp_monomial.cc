//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------

#include <base/quadrature.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_dgp_monomial.h>
#include <fe/fe_values.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


// TEST
#include <base/polynomial.h>



// namespace for some functions that are used in this file.
namespace 
{
				   // storage of hand-chosen support
				   // points
				   //
				   // For dim=2, dofs_per_cell of
				   // FE_DGPMonomial(k) is given by
				   // 0.5(k+1)(k+2), i.e.
				   // 
				   // k    0  1  2  3  4  5  6  7
				   // dofs 1  3  6 10 15 21 28 36
				   // 
				   // indirect access of unit points:
				   // the points for degree k are
				   // located at
				   //
				   // points[start_index[k]..start_index[k+1]-1]
  const unsigned int start_index2d[6]={0,1,4,10,20,35};
  const double points2d[35][2]=
  {{0,0},
   {0,0},{1,0},{0,1},
   {0,0},{1,0},{0,1},{1,1},{0.5,0},{0,0.5},
   {0,0},{1,0},{0,1},{1,1},{1./3.,0},{2./3.,0},{0,1./3.},{0,2./3.},{0.5,1},{1,0.5},
   {0,0},{1,0},{0,1},{1,1},{0.25,0},{0.5,0},{0.75,0},{0,0.25},{0,0.5},{0,0.75},{1./3.,1},{2./3.,1},{1,1./3.},{1,2./3.},{0.5,0.5}
  };

				   //
				   // For dim=3, dofs_per_cell of
				   // FE_DGPMonomial(k) is given by
				   // 1./6.(k+1)(k+2)(k+3), i.e.
				   // 
				   // k    0  1  2  3  4  5  6   7
				   // dofs 1  4 10 20 35 56 84 120
  const unsigned int start_index3d[6]={0,1,5,15/*,35*/};
  const double points3d[35][3]=
  {{0,0,0},
   {0,0,0},{1,0,0},{0,1,0},{0,0,1},
   {0,0,0},{1,0,0},{0,1,0},{0,0,1},{0.5,0,0},{0,0.5,0},{0,0,0.5},{1,1,0},{1,0,1},{0,1,1}
  };

  
  template<int dim>
  void generate_unit_points (const unsigned int,
			     std::vector<Point<dim> > &);

  template <>
  void generate_unit_points (const unsigned int,
			     std::vector<Point<1> > &)
  {
    Assert(false, ExcNotImplemented());
  }
  
  
  template <>
  void generate_unit_points (const unsigned int k,
			     std::vector<Point<2> > &p)
  {
    Assert(k<=4, ExcNotImplemented());
    Assert(p.size()==start_index2d[k+1]-start_index2d[k], ExcInternalError());
    for (unsigned int i=0; i<p.size(); ++i)
      {
	p[i](0)=points2d[start_index2d[k]+i][0];
	p[i](1)=points2d[start_index2d[k]+i][1];
      }
  }
  
  template <>
  void generate_unit_points (const unsigned int k,
			     std::vector<Point<3> > &p)
  {
    Assert(k<=2, ExcNotImplemented());
    Assert(p.size()==start_index3d[k+1]-start_index3d[k], ExcInternalError());
    for (unsigned int i=0; i<p.size(); ++i)
      {
	p[i](0)=points3d[start_index3d[k]+i][0];
	p[i](1)=points3d[start_index3d[k]+i][1];
	p[i](2)=points3d[start_index3d[k]+i][2];
      }
  }  
}



template <int dim>
FE_DGPMonomial<dim>::FE_DGPMonomial (const unsigned int degree)
		:
		FiniteElement<dim> (FiniteElementData<dim>(get_dpo_vector(degree), 1, degree),
				    std::vector<bool>(FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell,true),
				    std::vector<std::vector<bool> >(FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell,
								    std::vector<bool>(1,true))),
								      polynomial_space(degree)
{
  Assert(polynomial_space.n()==dofs_per_cell, ExcInternalError());
  
				   // DG doesn't have constraints, so
				   // leave them empty

				   // initialize the interpolation
				   // matrices
  initialize_embedding ();
//  initialize_restriction ();

                                   // note, that these elements have
                                   // neither support nor face-support
                                   // points, so leave these fields
                                   // empty
}



template <int dim>
std::string
FE_DGPMonomial<dim>::get_name () const
{
				   // note that the
				   // FETools::get_fe_from_name
				   // function depends on the
				   // particular format of the string
				   // this function returns, so they
				   // have to be kept in synch

#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream namebuf;
#else
  std::ostrstream namebuf;
#endif
  
  namebuf << "FE_DGPMonomial<" << dim << ">(" << get_degree() << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_DGPMonomial<dim>::clone() const
{
  return new FE_DGPMonomial<dim>(get_degree());
}



template <int dim>
double
FE_DGPMonomial<dim>::shape_value (const unsigned int i,
				  const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  return polynomial_space.compute_value(i, p);
}



template <int dim>
double
FE_DGPMonomial<dim>::shape_value_component (const unsigned int i,
					    const Point<dim> &p,
					    const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_value(i, p);
}



template <int dim>
Tensor<1,dim>
FE_DGPMonomial<dim>::shape_grad (const unsigned int i,
				 const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  return polynomial_space.compute_grad(i, p);
}


template <int dim>
Tensor<1,dim>
FE_DGPMonomial<dim>::shape_grad_component (const unsigned int i,
					   const Point<dim> &p,
					   const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_grad(i, p);
}



template <int dim>
Tensor<2,dim>
FE_DGPMonomial<dim>::shape_grad_grad (const unsigned int i,
				      const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  return polynomial_space.compute_grad_grad(i, p);
}



template <int dim>
Tensor<2,dim>
FE_DGPMonomial<dim>::shape_grad_grad_component (const unsigned int i,
						const Point<dim> &p,
						const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_grad_grad(i, p);
}



template <int dim>
void
FE_DGPMonomial<dim>::
get_interpolation_matrix (const FiniteElementBase<dim> &x_source_fe,
			  FullMatrix<double>           &interpolation_matrix) const
{
				   // this is only implemented, if the
				   // source FE is also a
				   // DGPMonomial element
  AssertThrow ((x_source_fe.get_name().find ("FE_DGPMonomial<") == 0)
               ||
               (dynamic_cast<const FE_DGPMonomial<dim>*>(&x_source_fe) != 0),
               typename FiniteElementBase<dim>::
               ExcInterpolationNotImplemented());
  
				   // ok, source is a Q element, so
				   // we will be able to do the work
  const FE_DGPMonomial<dim> &source_fe
    = dynamic_cast<const FE_DGPMonomial<dim>&>(x_source_fe);

  const unsigned int m=interpolation_matrix.m();
  const unsigned int n=interpolation_matrix.n();
  Assert (m == this->dofs_per_cell, ExcDimensionMismatch (m, this->dofs_per_cell));
  Assert (n == source_fe.dofs_per_cell, ExcDimensionMismatch (n, source_fe.dofs_per_cell));

  const unsigned int min_mn=
    interpolation_matrix.m()<interpolation_matrix.n() ?
    interpolation_matrix.m() : interpolation_matrix.n();

  for (unsigned int i=0; i<min_mn; ++i)
    interpolation_matrix(i,i)=1.;
}



template <int dim>
void
FE_DGPMonomial<dim>::initialize_embedding ()
{
  std::vector<Point<dim> > unit_points(this->dofs_per_cell);
  generate_unit_points(get_degree(), unit_points);
  
  FullMatrix<double> cell_interpolation (this->dofs_per_cell,
					 this->dofs_per_cell);
  FullMatrix<double> subcell_interpolation (this->dofs_per_cell,
					    this->dofs_per_cell);
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    this->prolongation[child].reinit (this->dofs_per_cell,
				      this->dofs_per_cell);
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
	{
	  const Point<dim> &p_subcell=unit_points[j];
	  
	  const Point<dim> p_cell =
	    GeometryInfo<dim>::child_to_cell_coordinates (p_subcell, child);
	  
	  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
	    {
	      cell_interpolation(j,i) = polynomial_space.compute_value (i, p_cell);
	      subcell_interpolation(j,i) = polynomial_space.compute_value (i, p_subcell);
	    }
	}
      
				       // then compute the embedding
				       // matrix for this child and
				       // this coordinate direction
      subcell_interpolation.gauss_jordan ();
      subcell_interpolation.mmult (this->prolongation[child], cell_interpolation);
      
				       // cut off very small values
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
	for (unsigned int j=0; j<this->dofs_per_cell; ++j)
	  if (std::fabs(this->prolongation[child](i,j)) < 2e-14*get_degree()*dim)
	    this->prolongation[child](i,j) = 0.;      
    }
}


template <int dim>
void
FE_DGPMonomial<dim>::initialize_restriction ()
{
  Assert(false, ExcNotImplemented());
}


//----------------------------------------------------------------------
// Auxiliary functions
//----------------------------------------------------------------------


template <int dim>
std::vector<unsigned int>
FE_DGPMonomial<dim>::get_dpo_vector(unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[dim] = ++deg;
  for (unsigned int i=1;i<dim;++i)
    {
      dpo[dim] *= deg+i;
      dpo[dim] /= i+1;
    }
  return dpo;
}


template <int dim>
UpdateFlags
FE_DGPMonomial<dim>::update_once (const UpdateFlags flags) const
{
				   // for this kind of elements, only
				   // the values can be precomputed
				   // once and for all. set this flag
				   // if the values are requested at
				   // all
  return (update_default | (flags & update_values));
}


template <int dim>
UpdateFlags
FE_DGPMonomial<dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;

  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_covariant_transformation;

  return out;
}


//----------------------------------------------------------------------
// Data field initialization
//----------------------------------------------------------------------

template <int dim>
typename Mapping<dim>::InternalDataBase *
FE_DGPMonomial<dim>::get_data (const UpdateFlags      update_flags,
			       const Mapping<dim>    &mapping,
			       const Quadrature<dim> &quadrature) const
{
				   // generate a new data object
  InternalData* data = new InternalData;
				   // check what needs to be
				   // initialized only once and what
				   // on every cell/face/subface we
				   // visit
  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  const UpdateFlags flags(data->update_flags);
  const unsigned int n_q_points = quadrature.n_quadrature_points;

				   // have some scratch arrays
  std::vector<double> values(0);
  std::vector<Tensor<1,dim> > grads(0);
  std::vector<Tensor<2,dim> > grad_grads(0);
  
				   // initialize fields only if really
				   // necessary. otherwise, don't
				   // allocate memory
  if (flags & update_values)
    {
      values.resize (this->dofs_per_cell);
      data->shape_values.reinit (this->dofs_per_cell,
				 n_q_points);
    }

  if (flags & update_gradients)
    {
      grads.resize (this->dofs_per_cell);
      data->shape_gradients.reinit (this->dofs_per_cell,
				    n_q_points);
    }

				   // if second derivatives through
				   // finite differencing is required,
				   // then initialize some objects for
				   // that
  if (flags & update_second_derivatives)
    data->initialize_2nd (this, mapping, quadrature);
  
				   // next already fill those fields
				   // of which we have information by
				   // now. note that the shape
				   // gradients are only those on the
				   // unit cell, and need to be
				   // transformed when visiting an
				   // actual cell  
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
	polynomial_space.compute(quadrature.point(i),
				 values, grads, grad_grads);
	for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	  {
	    if (flags & update_values)
	      data->shape_values[k][i] = values[k];
	    if (flags & update_gradients)
	      data->shape_gradients[k][i] = grads[k];
	  }
      }
  return data;
}



//----------------------------------------------------------------------
// Fill data of FEValues
//----------------------------------------------------------------------

template <int dim>
void
FE_DGPMonomial<dim>::fill_fe_values (const Mapping<dim>                   &mapping,
			     const typename DoFHandler<dim>::cell_iterator &cell,
			     const Quadrature<dim>                &quadrature,
			     typename Mapping<dim>::InternalDataBase       &mapping_data,
			     typename Mapping<dim>::InternalDataBase       &fedata,
			     FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
	for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());
	  mapping.transform_covariant(data.shape_gradients[k].begin(),
				      data.shape_gradients[k].end(),
				      fe_data.shape_gradients[k].begin(),
				      mapping_data);
	}
    }
  
  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell,
                       QProjector<dim>::DataSetDescriptor::cell(),
                       mapping_data, fe_data, data);
}



template <int dim>
void
FE_DGPMonomial<dim>::fill_fe_face_values (const Mapping<dim>                   &mapping,
				  const typename DoFHandler<dim>::cell_iterator &cell,
				  const unsigned int                    face,
				  const Quadrature<dim-1>              &quadrature,
				  typename Mapping<dim>::InternalDataBase       &mapping_data,
				  typename Mapping<dim>::InternalDataBase       &fedata,
				  FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);

				   // offset determines which data set
				   // to take (all data sets for all
				   // faces are stored contiguously)
  const typename QProjector<dim>::DataSetDescriptor offset
    = (QProjector<dim>::DataSetDescriptor::
       face (face, cell->face_orientation(face),
             quadrature.n_quadrature_points));
  
  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	if (flags & update_values)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() + offset <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());	  
	  mapping.transform_covariant(data.shape_gradients[k].begin(),
				      data.shape_gradients[k].end(),
				      fe_data.shape_gradients[k].begin()+offset,
				      mapping_data);
	};
    }

  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <int dim>
void
FE_DGPMonomial<dim>::fill_fe_subface_values (const Mapping<dim>                   &mapping,
				     const typename DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int                    face,
				     const unsigned int                    subface,
				     const Quadrature<dim-1>              &quadrature,
				     typename Mapping<dim>::InternalDataBase       &mapping_data,
				     typename Mapping<dim>::InternalDataBase       &fedata,
				     FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  
				   // offset determines which data set
				   // to take (all data sets for all
				   // sub-faces are stored contiguously)
  const typename QProjector<dim>::DataSetDescriptor offset
    = (QProjector<dim>::DataSetDescriptor::
       sub_face (face, subface, cell->face_orientation(face), 
                 quadrature.n_quadrature_points));

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	if (flags & update_values)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() + offset <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());	  
	  mapping.transform_covariant(data.shape_gradients[k].begin(),
				      data.shape_gradients[k].end(),
				      fe_data.shape_gradients[k].begin()+offset,
				      mapping_data);
	};
    }
  
  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <int dim>
unsigned int
FE_DGPMonomial<dim>::n_base_elements () const
{
  return 1;
}



template <int dim>
const FiniteElement<dim> &
FE_DGPMonomial<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <int dim>
unsigned int
FE_DGPMonomial<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}


#if deal_II_dimension == 1

template <>
bool
FE_DGPMonomial<1>::has_support_on_face (const unsigned int,
					const unsigned int face_index) const
{
  return face_index==1 || (face_index==0 && get_degree()==0);
}

#endif

#if deal_II_dimension == 2

template <>
bool
FE_DGPMonomial<2>::has_support_on_face (const unsigned int shape_index,
					const unsigned int face_index) const
{
  bool support_on_face=false;
  if (face_index==1 || face_index==2)
    support_on_face=true;
  else
    {
      unsigned int degrees[2];
      polynomial_space.directional_degrees(shape_index, degrees);
      if ((face_index==0 && degrees[1]==0) ||
	  (face_index==3 && degrees[0]==0))
	support_on_face=true;
    }
  return support_on_face;
}

#endif

#if deal_II_dimension == 3

template <>
bool
FE_DGPMonomial<3>::has_support_on_face (const unsigned int shape_index,
					const unsigned int face_index) const
{
  bool support_on_face=false;
  if (face_index==1 || face_index==3 || face_index==4)
    support_on_face=true;
  else
    {
      unsigned int degrees[3];
      polynomial_space.directional_degrees(shape_index, degrees);
      if ((face_index==0 && degrees[1]==0) ||
	  (face_index==2 && degrees[2]==0) ||
	  (face_index==5 && degrees[0]==0))
	support_on_face=true;
    }
  return support_on_face;
}

#endif


template <int dim>
unsigned int
FE_DGPMonomial<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template class FE_DGPMonomial<deal_II_dimension>;
