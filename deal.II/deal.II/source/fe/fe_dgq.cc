//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------

#include <base/quadrature.h>
#include <base/polynomial.h>
#include <base/tensor_product_polynomials.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_dgq.h>
#include <fe/fe_values.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif



// namespace for some functions that are used in this file. they are
// specific to numbering conventions used for the FE_DGQ element, and
// are thus not very interesting to the outside world
namespace
{
				   // auxiliary type to allow for some
				   // kind of explicit template
				   // specialization of the following
				   // functions
  template <int dim> struct int2type {};


				   // given an integer N, compute its
				   // integer square root (if it
				   // exists, otherwise give up)
  unsigned int int_sqrt (const unsigned int N)
  {
    for (unsigned int i=0; i<=N; ++i)
      if (i*i == N)
	return i;
    Assert (false, ExcInternalError());
    return static_cast<unsigned int>(-1);
  }


				   // given an integer N, compute its
				   // integer cube root (if it
				   // exists, otherwise give up)
  unsigned int int_cuberoot (const unsigned int N)
  {
    for (unsigned int i=0; i<=N; ++i)
      if (i*i*i == N)
	return i;
    Assert (false, ExcInternalError());
    return static_cast<unsigned int>(-1);
  }


				   // given N, generate i=0...N-1
				   // equidistant points in the
				   // interior of the interval [0,1]
  Point<1>
  generate_unit_point (const unsigned int i,
		       const unsigned int N,
		       const int2type<1>  )
  {
    Assert (i<N, ExcInternalError());
    if (N==1)
      return Point<1> (.5);
    else
      {
        const double h = 1./(N-1);
        return Point<1>(i*h);
      }
  }
  
    
				   // given N, generate i=0...N-1
				   // equidistant points in the domain
				   // [0,1]^2
  Point<2>
  generate_unit_point (const unsigned int i,
		       const unsigned int N,
		       const int2type<2>  )
  {
    Assert (i<N, ExcInternalError());
    
    if (N==1)
      return Point<2> (.5, .5);
    else
      {
        Assert (N>=4, ExcInternalError());
        const unsigned int N1d = int_sqrt(N);
        const double h = 1./(N1d-1);
        
        return Point<2> (i%N1d * h,
                         i/N1d * h);
      }
  }
  

  

				   // given N, generate i=0...N-1
				   // equidistant points in the domain
				   // [0,1]^3
  Point<3>
  generate_unit_point (const unsigned int i,
		       const unsigned int N,
		       const int2type<3>  )
  {
    Assert (i<N, ExcInternalError());
    if (N==1)
      return Point<3> (.5, .5, .5);
    else
      {
        Assert (N>=8, ExcInternalError());
    
        const unsigned int N1d = int_cuberoot(N);
        const double h = 1./(N1d-1);

        return Point<3> (i%N1d * h,
                         (i/N1d)%N1d * h,
                         i/(N1d*N1d) * h);
      }
  }
}




template <int dim>
FE_DGQ<dim>::FE_DGQ (const unsigned int degree)
		:
		FiniteElement<dim> (FiniteElementData<dim>(get_dpo_vector(degree),1),
				    std::vector<bool>(FiniteElementData<dim>(get_dpo_vector(degree),1).dofs_per_cell,
                                                      true),
				    std::vector<std::vector<bool> >(FiniteElementData<dim>(get_dpo_vector(degree),1).dofs_per_cell,
								    std::vector<bool>(1,true))),
                degree(degree),
                polynomial_space (Polynomials::LagrangeEquidistant::generate_complete_basis(degree))
{
				   // generate permutation/rotation
				   // index sets to generate some
				   // matrices from others
  std::vector<unsigned int> right;
  std::vector<unsigned int> top;
  rotate_indices (right, 'Z');
  if (dim>2)
    rotate_indices (top, 'X');

				   // if defined, copy over matrices
				   // from precomputed arrays and
				   // generate all other matrices by
				   // permutations
                                   //
                                   // (note: the matrix is defined if
                                   // something was entered into the
                                   // respective table, and what was
                                   // entered is not a NULL pointer --
                                   // this would allow for "holes")
  if ((degree < Matrices::n_embedding_matrices) &&
      (Matrices::embedding[degree] != 0))
    {
      for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
        this->prolongation[i].reinit (this->dofs_per_cell,
                                      this->dofs_per_cell);
      this->prolongation[0].fill (Matrices::embedding[degree]);
      switch (dim)
	{
	  case 1:
          {
	    this->prolongation[1].fill_permutation (this->prolongation[0],
						    right, right);
	    break;
          };
          
	  case 2:
          {
	    this->prolongation[1].fill_permutation (this->prolongation[0],
						    right, right);
	    this->prolongation[2].fill_permutation (this->prolongation[1],
						    right, right);
	    this->prolongation[3].fill_permutation (this->prolongation[2],
						    right, right);
	    break;
          };
          
	  case 3:
          {
	    this->prolongation[1].fill_permutation (this->prolongation[0],
						    right, right);
	    this->prolongation[5].fill_permutation (this->prolongation[1],
						    right, right);
	    this->prolongation[4].fill_permutation (this->prolongation[5],
						    right, right);
	    this->prolongation[7].fill_permutation (this->prolongation[4],
						    top, top);
	    this->prolongation[3].fill_permutation (this->prolongation[7],
						    top, top);
	    this->prolongation[6].fill_permutation (this->prolongation[5],
						    top, top);
	    this->prolongation[2].fill_permutation (this->prolongation[6],
						    top, top);
	    break;
          };
          
	  default:
	    Assert (false, ExcNotImplemented());
	}
    }
  else
				     // matrix undefined, leave matrix
				     // at size zero
    ;

				   // same as above: copy over matrix
				   // from predefined values and
				   // generate all others by rotation
  if ((degree < Matrices::n_projection_matrices) &&
      (Matrices::projection_matrices[degree] != 0))
    {
      for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
        this->restriction[i].reinit (this->dofs_per_cell,
                                     this->dofs_per_cell);
      this->restriction[0].fill (Matrices::projection_matrices[degree]);
      switch (dim)
	{
	  case 1:
          {
	    this->restriction[1].fill_permutation (this->restriction[0],
						   right, right);
	    break;
          };
          
	  case 2:
          {
	    this->restriction[1].fill_permutation (this->restriction[0],
						   right, right);
	    this->restriction[2].fill_permutation (this->restriction[1],
						   right, right);
	    this->restriction[3].fill_permutation (this->restriction[2],
						   right, right);
	    break;
          };
          
	  case 3:
          {
	    this->restriction[1].fill_permutation (this->restriction[0],
						   right, right);
	    this->restriction[5].fill_permutation (this->restriction[1],
						   right, right);
	    this->restriction[4].fill_permutation (this->restriction[5],
						   right, right);
	    this->restriction[7].fill_permutation (this->restriction[4],
						   top, top);
	    this->restriction[3].fill_permutation (this->restriction[7],
						   top, top);
	    this->restriction[6].fill_permutation (this->restriction[5],
						   top, top);
	    this->restriction[2].fill_permutation (this->restriction[6],
						   top, top);
	    break;
          };
          
	  default:
	    Assert (false, ExcNotImplemented());
	}
    }
  else
				     // matrix undefined, leave matrix
				     // at size zero
    ;
  
				   // finally fill in support points
  if (degree == 0)
    {
				       // constant elements, take
				       // midpoint
      this->unit_support_points.resize(1);
      for (unsigned int i=0; i<dim; ++i)
	this->unit_support_points[0](i) = 0.5;
    }
  else
    {
				       // number of points: (degree+1)^dim
      unsigned int n = degree+1;
      for (unsigned int i=1; i<dim; ++i)
	n *= degree+1;
      
      this->unit_support_points.resize(n);
      
      const double step = 1./degree;
      Point<dim> p;
      
      unsigned int k=0;
      for (unsigned int iz=0; iz <= ((dim>2) ? degree : 0) ; ++iz)
	for (unsigned int iy=0; iy <= ((dim>1) ? degree : 0) ; ++iy)
	  for (unsigned int ix=0; ix<=degree; ++ix)
	    {
	      p(0) = ix * step;
	      if (dim>1)
		p(1) = iy * step;
	      if (dim>2)
		p(2) = iz * step;
	      
	      this->unit_support_points[k++] = p;
	    };
    };

				   // note: no face support points for
				   // DG elements
}



template <int dim>
std::string
FE_DGQ<dim>::get_name () const
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
  
  namebuf << "FE_DGQ<" << dim << ">(" << degree << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_DGQ<dim>::clone() const
{
  return new FE_DGQ<dim>(degree);
}



template <int dim>
double
FE_DGQ<dim>::shape_value (const unsigned int i,
			  const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  return polynomial_space.compute_value(i, p);
}



template <int dim>
double
FE_DGQ<dim>::shape_value_component (const unsigned int i,
				    const Point<dim> &p,
				    const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_value(i, p);
}



template <int dim>
Tensor<1,dim>
FE_DGQ<dim>::shape_grad (const unsigned int i,
			 const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  return polynomial_space.compute_grad(i, p);
}


template <int dim>
Tensor<1,dim>
FE_DGQ<dim>::shape_grad_component (const unsigned int i,
				   const Point<dim> &p,
				   const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_grad(i, p);
}



template <int dim>
Tensor<2,dim>
FE_DGQ<dim>::shape_grad_grad (const unsigned int i,
			      const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  return polynomial_space.compute_grad_grad(i, p);
}



template <int dim>
Tensor<2,dim>
FE_DGQ<dim>::shape_grad_grad_component (const unsigned int i,
					const Point<dim> &p,
					const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i, 0, this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_grad_grad(i, p);
}


//----------------------------------------------------------------------
// Auxiliary functions
//----------------------------------------------------------------------


template <int dim>
std::vector<unsigned int>
FE_DGQ<dim>::get_dpo_vector(unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, static_cast<unsigned int>(0));
  dpo[dim] = ++deg;
  for (unsigned int i=1;i<dim;++i)
    dpo[dim] *= deg;
  return dpo;
}


template <int dim>
UpdateFlags
FE_DGQ<dim>::update_once (const UpdateFlags flags) const
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
FE_DGQ<dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;

  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_covariant_transformation;

  return out;
}



template <int dim>
void
FE_DGQ<dim>::rotate_indices (std::vector<unsigned int> &numbers,
			     const char                 direction) const
{
  const unsigned int n = degree+1;
  unsigned int s = n;
  for (unsigned int i=1;i<dim;++i)
    s *= n;
  numbers.resize (s);
  
  unsigned int l = 0;

  if (dim==1)
    {
				       // Mirror around midpoint
      for (unsigned int i=n;i>0;)
	numbers[l++]=--i;
    }
  else
    {
      switch (direction)
	{
					   // Rotate xy-plane
					   // counter-clockwise
	  case 'z':
	    for (unsigned int iz=0;iz<((dim>2) ? n:1);++iz)
	      for (unsigned int j=0;j<n;++j)
		for (unsigned int i=0;i<n;++i)
		  {
		    unsigned int k = n*i-j+n-1 + n*n*iz;
		    numbers[l++] = k;
		  }
	    break;
					     // Rotate xy-plane
					     // clockwise
	  case 'Z':
	    for (unsigned int iz=0;iz<((dim>2) ? n:1);++iz)
	      for (unsigned int iy=0;iy<n;++iy)
		for (unsigned int ix=0;ix<n;++ix)
		  {
		    unsigned int k = n*ix-iy+n-1 + n*n*iz;
		    numbers[k] = l++;
		  }
	    break;
					     // Rotate yz-plane
					     // counter-clockwise
	  case 'x':
	    Assert (dim>2,
		    typename FiniteElementData<dim>::ExcSpaceDimensionMismatch (dim,3));
	    for (unsigned int iz=0;iz<n;++iz)
	      for (unsigned int iy=0;iy<n;++iy)
		for (unsigned int ix=0;ix<n;++ix)
		  {
		    unsigned int k = n*(n*iy-iz+n-1) + ix;
		    numbers[l++] = k;
		  }
	    break;
					     // Rotate yz-plane
					     // clockwise
	  case 'X':
	    Assert (dim>2, 
		    typename FiniteElementData<dim>::ExcSpaceDimensionMismatch (dim,3));
	    for (unsigned int iz=0;iz<n;++iz)
	      for (unsigned int iy=0;iy<n;++iy)
		for (unsigned int ix=0;ix<n;++ix)
		  {
		    unsigned int k = n*(n*iy-iz+n-1) + ix;
		    numbers[k] = l++;
		  }
	    break;
	  default:
	    Assert (false, ExcNotImplemented ());
	}
    }
}



template <int dim>
void
FE_DGQ<dim>::
get_interpolation_matrix (const FiniteElementBase<dim> &x_source_fe,
			  FullMatrix<double>           &interpolation_matrix) const
{
				   // this is only implemented, if the
				   // source FE is also a
				   // DGQ element
  AssertThrow ((x_source_fe.get_name().find ("FE_DGQ<") == 0)
               ||
               (dynamic_cast<const FE_DGQ<dim>*>(&x_source_fe) != 0),
               typename FiniteElementBase<dim>::
               ExcInterpolationNotImplemented());
  
				   // ok, source is a Q element, so
				   // we will be able to do the work
  const FE_DGQ<dim> &source_fe
    = dynamic_cast<const FE_DGQ<dim>&>(x_source_fe);

  Assert (interpolation_matrix.m() == this->dofs_per_cell,
	  ExcDimensionMismatch (interpolation_matrix.m(),
				this->dofs_per_cell));
  Assert (interpolation_matrix.n() == source_fe.dofs_per_cell,
	  ExcDimensionMismatch (interpolation_matrix.m(),
				source_fe.dofs_per_cell));
  
  
				   // compute the interpolation
				   // matrices in much the same way as
				   // we do for the embedding matrices
				   // from mother to child.
  FullMatrix<double> cell_interpolation (this->dofs_per_cell,
					 this->dofs_per_cell);
  FullMatrix<double> source_interpolation (this->dofs_per_cell,
					   source_fe.dofs_per_cell);
  FullMatrix<double> tmp (this->dofs_per_cell,
			  source_fe.dofs_per_cell);
  for (unsigned int j=0; j<this->dofs_per_cell; ++j)
    {
                                       // generate a point on this
                                       // cell and evaluate the
                                       // shape functions there
      const Point<dim> p = generate_unit_point (j, this->dofs_per_cell,
                                                int2type<dim>());
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        cell_interpolation(j,i)
          = polynomial_space.compute_value (i, p);

      for (unsigned int i=0; i<source_fe.dofs_per_cell; ++i)
        source_interpolation(j,i)
          = source_fe.polynomial_space.compute_value (i, p);
    }

                                   // then compute the
                                   // interpolation matrix matrix
                                   // for this coordinate
                                   // direction
  cell_interpolation.gauss_jordan ();
  cell_interpolation.mmult (interpolation_matrix,
                            source_interpolation);

                                   // cut off very small values
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    for (unsigned int j=0; j<source_fe.dofs_per_cell; ++j)
      if (std::fabs(interpolation_matrix(i,j)) < 1e-15)
        interpolation_matrix(i,j) = 0.;

				   // make sure that the row sum of
				   // each of the matrices is 1 at
				   // this point. this must be so
				   // since the shape functions sum up
				   // to 1
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    {
      double sum = 0.;
      for (unsigned int j=0; j<source_fe.dofs_per_cell; ++j)
        sum += interpolation_matrix(i,j);

      Assert (std::fabs(sum-1) < 5e-14*std::max(degree,1U)*dim,
              ExcInternalError());
    }
}





//----------------------------------------------------------------------
// Data field initialization
//----------------------------------------------------------------------

template <int dim>
typename Mapping<dim>::InternalDataBase *
FE_DGQ<dim>::get_data (const UpdateFlags      update_flags,
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
FE_DGQ<dim>::fill_fe_values (const Mapping<dim>                   &mapping,
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
	};
    }
  
  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, 0, mapping_data, fe_data, data);
}



template <int dim>
void
FE_DGQ<dim>::fill_fe_face_values (const Mapping<dim>                   &mapping,
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
  const unsigned int offset = face * quadrature.n_quadrature_points;
  
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
FE_DGQ<dim>::fill_fe_subface_values (const Mapping<dim>                   &mapping,
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
  const unsigned int offset = (face * GeometryInfo<dim>::subfaces_per_face + subface)
			      * quadrature.n_quadrature_points;

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
FE_DGQ<dim>::n_base_elements () const
{
  return 1;
}



template <int dim>
const FiniteElement<dim> &
FE_DGQ<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <int dim>
unsigned int
FE_DGQ<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}



template <int dim>
bool
FE_DGQ<dim>::has_support_on_face (const unsigned int shape_index,
				  const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

  unsigned int n = degree+1;
  unsigned int n2 = n*n;
  
  switch (dim)
    {
      case 1:
    {
				       // in 1d, things are simple. since
				       // there is only one degree of
				       // freedom per vertex in this
				       // class, the first is on vertex 0
				       // (==face 0 in some sense), the
				       // second on face 1:
      return (((shape_index == 0) && (face_index == 0)) ||
	      ((shape_index == 1) && (face_index == 1)));
    };
      
      case 2:
    {
      if (face_index==0 && shape_index < n)
	return true;
      if (face_index==1 && (shape_index % n) == degree)
	return true;
      if (face_index==2 && shape_index >= this->dofs_per_cell-n)
	return true;
      if (face_index==3 && (shape_index % n) == 0)
	return true;
      return false;
    };
      
      case 3:
    {
      const unsigned int in2 = shape_index % n2;
	
				       // y=0
      if (face_index==0 && in2 < n )
	return true;
				       // y=1
      if (face_index==1 && in2 >= n2-n)
	return true;
				       // z=0
      if (face_index==2 && shape_index < n2)
	return true;
				       // x=1
      if (face_index==3 && (shape_index % n) == n-1)
	return true;
				       // z=1
      if (face_index==4 && shape_index >= this->dofs_per_cell - n2)
	return true;
				       // x=0
      if (face_index==5 && (shape_index % n) == 0)
	return true;
      return false;
    };

      default:
	Assert (false, ExcNotImplemented());
    }
  return true;
}



template <int dim>
unsigned int
FE_DGQ<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template <int dim>
unsigned int
FE_DGQ<dim>::get_degree () const
{
  return degree;
}




template class FE_DGQ<deal_II_dimension>;
