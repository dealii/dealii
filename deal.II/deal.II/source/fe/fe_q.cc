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
#include <fe/fe_q.h>
#include <fe/fe_values.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


// namespace for some functions that are used in this file. they are
// specific to numbering conventions used for the FE_Q element, and
// are thus not very interesting to the outside world
namespace
{
				   // auxiliary type to allow for some
				   // kind of explicit template
				   // specialization of the following
				   // functions
  template <int dim> struct int2type {};

				   // given a permutation array,
				   // compute and return the inverse
				   // permutation
#ifdef DEAL_II_ANON_NAMESPACE_BUG
  static
#endif
  std::vector<unsigned int>
  invert_numbering (const std::vector<unsigned int> &in)
  {
    std::vector<unsigned int> out (in.size());
    for (unsigned int i=0; i<in.size(); ++i)
      out[in[i]]=i;
    return out;
  }


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
    const double h = 1./(N-1);
    return Point<1>(i*h);
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
    Assert (N>=4, ExcInternalError());
    
    const unsigned int N1d = int_sqrt(N);
    const double h = 1./(N1d-1);

    return Point<2> (i%N1d * h,
		     i/N1d * h);
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
    Assert (N>=8, ExcInternalError());
    
    const unsigned int N1d = int_cuberoot(N);
    const double h = 1./(N1d-1);

    return Point<3> (i%N1d * h,
		     (i/N1d)%N1d * h,
		     i/(N1d*N1d) * h);
  }

  

				   // given N, generate i=0...N-1
				   // equidistant points in the
				   // interior of the interval [0,1]
  Point<1>
  generate_face_unit_point (const unsigned int i,
			    const unsigned int N,
			    const int2type<1>  )
  {
    Assert (i<N, ExcInternalError());
    const double h = 1./(N+1);
    return Point<1>((1+i)*h);
  }
  
    
				   // given N, generate i=0...N-1
				   // equidistant points in the domain
				   // [0,1]^2, but excluding the four
				   // vertices (since we don't have to
				   // consider shape functions on
				   // child cells that are located on
				   // existing vertices)
  Point<2>
  generate_face_unit_point (const unsigned int i,
			    const unsigned int N,
			    const int2type<2>  )
  {
    Assert (i<N, ExcInternalError());
    
    const unsigned int N1d = int_sqrt(N+4);    
    const double h = 1./(N1d+1);
    
				     // i gives the index in the list
				     // of points excluding the four
				     // vertices. convert this into an
				     // index for all N1d**2 points
				     //
				     // we do so by
				     // - adding one if the point is
				     // beyond the lower left vertex
				     // (actually, all points are)

				     // - adding one if the point is
				     // beyond the lower right one
				     // - adding one if it is beyond
				     // the upper left one
				     // - not adding one for the upper
				     // right vertex, since no point
				     // can be beyond that one anyway
				     // :-)
    const unsigned int true_i = (1
				 +
				 (i >= N1d-2 ? 1 : 0)
				 +
				 (i >= N1d*(N1d-1)-2 ? 1 : 0));
    return Point<2> ((true_i%N1d)*h,
		     (true_i/N1d)*h);
  }
  


				   // return whether shape function j,
				   // as given in the numbering
				   // specific to the computation of
				   // the constraint matrix, is active
				   // on the given subface
  bool
  constraint_function_is_active_on_child (const unsigned int          j,
					  const unsigned int          subface,
					  const FiniteElementData<2> &fe_data)
  {
				     // note that in our weird
				     // numbering, the zeroth function
				     // is the one associated with the
				     // center node, then come the
				     // ones on subface 0, then those
				     // on subface 1. the initial one
				     // is active on both subfaces,
				     // all other ones only on one of
				     // the subfaces
    return !(((j>=1) && (j<1+fe_data.dofs_per_line) && (subface == 1)) ||
	     ((j>=1+fe_data.dofs_per_line) && (subface == 0)));
  }



  bool
  constraint_function_is_active_on_child (const unsigned int          j,
					  const unsigned int          subface,
					  const FiniteElementData<3> &fe_data)
  {
				     // in 3d: in our weird numbering,
				     // the zeroth function is the one
				     // associated with the center
				     // node, then come the four edge
				     // midpoints, then the ones on
				     // the 12 edges then those on
				     // subfaces. some are active on
				     // more than one child

    if (j < 5)
				       // one one of the five vertices
      {
	switch (j)
	  {
	    case 0: return true;
	    case 1: return (subface == 0) || (subface == 1);
	    case 2: return (subface == 1) || (subface == 2);
	    case 3: return (subface == 2) || (subface == 3);
	    case 4: return (subface == 3) || (subface == 0);
	  }
      }
    else if (j < 5 + 12*fe_data.dofs_per_line)
				       // one one of the 12 lines
      {
	const unsigned int line = (j-5)/fe_data.dofs_per_line;
	Assert (line<12, ExcInternalError());

	switch (line)
	  {
	    case 0: return (subface == 0) || (subface == 1);
	    case 1: return (subface == 1) || (subface == 2);
	    case 2: return (subface == 2) || (subface == 3);
	    case 3: return (subface == 3) || (subface == 0);
	    case 4: return (subface == 0);
	    case 5: return (subface == 1);
	    case 6: return (subface == 1);
	    case 7: return (subface == 2);
	    case 8: return (subface == 3);
	    case 9: return (subface == 2);
	    case 10: return (subface == 0);
	    case 11: return (subface == 2);
	  }
      }
    else
				       // interior
      {
	const unsigned int quad = (j-5-12*fe_data.dofs_per_line)/fe_data.dofs_per_quad;
	Assert (quad<4, ExcInternalError());
	return quad == subface;
      }
    
    Assert (false, ExcInternalError());
    return static_cast<unsigned int>(-1);
  }


				   // given index j in the weird
				   // constraint numbering, compute
				   // its index in the polynomials
				   // space of a given subface
  unsigned int
  constraint_get_local_j (const unsigned int          j,
			  const unsigned int          subface,
			  const FiniteElementData<2> &fe_data)
  {
				     // the zeroth shape function is a
				     // little special, since it has
				     // index N on subface 0 and index
				     // 0 on subface 1
	  
    return (subface == 0 ?
	    (j == 0 ? 1+fe_data.dofs_per_line : j) :
	    (j == 0 ? 0 : j-fe_data.dofs_per_line));
  }
  

  unsigned int
  constraint_get_local_j (const unsigned int          /*j*/,
			  const unsigned int          /*subface*/,
			  const FiniteElementData<3> &/*fe_data*/)
  {
    Assert (false, ExcNotImplemented());
//    const unsigned int N1d = 2+fe_data.dofs_per_line;
    return static_cast<unsigned int>(-1);
  }



				   // in the constraint numbering:
				   // return true if the support point
				   // of shape function j and
				   // evaluation point i coincide. to
				   // make things simpler, also pass
				   // the subface on which j is
				   // located
  bool
  constraint_is_support_point (const unsigned int          i,
			       const unsigned int          j,
			       const unsigned int          subface,
			       const FiniteElementData<2> &fe_data)
  {
    return ((subface == 0) && (((j==0) && (i==fe_data.dofs_per_line))
			       ||
			       ((j!=0) && (i==j-1))))
		      ||
	   ((subface == 1) && (((j==0) && (i==fe_data.dofs_per_line))
			       ||
			       ((j!=0) && (i==j))));
  }
  

  bool
  constraint_is_support_point (const unsigned int          /*i*/,
			       const unsigned int          /*j*/,
			       const unsigned int          /*subface*/,
			       const FiniteElementData<3> &/*fe_data*/)
  {
    Assert (false, ExcNotImplemented());
    return false;
  }  
}



template <int dim>
FE_Q<dim>::FE_Q (const unsigned int degree)
		:
		FiniteElement<dim> (FiniteElementData<dim>(get_dpo_vector(degree),1),
				    std::vector<bool> (FiniteElementData<dim>(get_dpo_vector(degree),1).dofs_per_cell,
                                                       false),
				    std::vector<std::vector<bool> >(FiniteElementData<dim>(get_dpo_vector(degree),1).dofs_per_cell,
								    std::vector<bool>(1,true))),
                degree(degree),
                renumber(lexicographic_to_hierarchic_numbering (*this, degree)),
		renumber_inverse(invert_numbering (renumber)),
		face_renumber(face_lexicographic_to_hierarchic_numbering (degree)),
		polynomial_space(Polynomials::LagrangeEquidistant::generate_complete_basis(degree))
{
  
				   // copy constraint and embedding
				   // matrices if they are
				   // defined. otherwise leave them at
				   // invalid size
  initialize_constraints ();
  initialize_embedding ();
  initialize_restriction ();

				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();
}



template <int dim>
std::string
FE_Q<dim>::get_name () const
{
#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream namebuf;
#else
  std::ostrstream namebuf;
#endif
  
  namebuf << "FE_Q<" << dim << ">(" << degree << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_Q<dim>::clone() const
{
  return new FE_Q<dim>(degree);
}



template <int dim>
double
FE_Q<dim>::shape_value (const unsigned int i,
			const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return polynomial_space.compute_value(renumber_inverse[i], p);
}


template <int dim>
double
FE_Q<dim>::shape_value_component (const unsigned int i,
				  const Point<dim> &p,
				  const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_value(renumber_inverse[i], p);
}



template <int dim>
Tensor<1,dim>
FE_Q<dim>::shape_grad (const unsigned int i,
		       const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return polynomial_space.compute_grad(renumber_inverse[i], p);
}



template <int dim>
Tensor<1,dim>
FE_Q<dim>::shape_grad_component (const unsigned int i,
				 const Point<dim> &p,
				 const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_grad(renumber_inverse[i], p);
}



template <int dim>
Tensor<2,dim>
FE_Q<dim>::shape_grad_grad (const unsigned int i,
			    const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return polynomial_space.compute_grad_grad(renumber_inverse[i], p);
}



template <int dim>
Tensor<2,dim>
FE_Q<dim>::shape_grad_grad_component (const unsigned int i,
				      const Point<dim> &p,
				      const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_grad_grad(renumber_inverse[i], p);
}



template <int dim>
void
FE_Q<dim>::
get_interpolation_matrix (const FiniteElementBase<dim> &x_source_fe,
			  FullMatrix<double>           &interpolation_matrix) const
{
				   // this is only implemented, if the
				   // source FE is also a
				   // Q element
  AssertThrow ((x_source_fe.get_name().find ("FE_Q<") == 0)
               ||
               (dynamic_cast<const FE_Q<dim>*>(&x_source_fe) != 0),
               typename FiniteElementBase<dim>::
               ExcInterpolationNotImplemented());
  
				   // ok, source is a Q element, so
				   // we will be able to do the work
  const FE_Q<dim> &source_fe
    = dynamic_cast<const FE_Q<dim>&>(x_source_fe);

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
        cell_interpolation(renumber[j],renumber[i])
          = polynomial_space.compute_value (i, p);

      for (unsigned int i=0; i<source_fe.dofs_per_cell; ++i)
        source_interpolation(renumber[j],source_fe.renumber[i])
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

      Assert (std::fabs(sum-1) < 2e-14*degree*dim,
              ExcInternalError());
    }
}


//----------------------------------------------------------------------
// Auxiliary functions
//----------------------------------------------------------------------



template <int dim>
void FE_Q<dim>::initialize_unit_support_points ()
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
	  
	  this->unit_support_points[renumber[k++]] = p;
	};
}


#if deal_II_dimension == 1

template <>
void FE_Q<1>::initialize_unit_face_support_points ()
{
				   // no faces in 1d, so nothing to do
}

#endif


template <int dim>
void FE_Q<dim>::initialize_unit_face_support_points ()
{
  const unsigned int codim = dim-1;
  
				   // number of points: (degree+1)^codim
  unsigned int n = degree+1;
  for (unsigned int i=1; i<codim; ++i)
    n *= degree+1;
  
  this->unit_face_support_points.resize(n);
  
  const double step = 1./degree;
  Point<codim> p;
  
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((codim>2) ? degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((codim>1) ? degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=degree; ++ix)
	{
	  p(0) = ix * step;
	  if (codim>1)
	    p(1) = iy * step;
	  if (codim>2)
	    p(2) = iz * step;
	  
	  this->unit_face_support_points[face_renumber[k++]] = p;
	};
}



template <int dim>
std::vector<unsigned int>
FE_Q<dim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, static_cast<unsigned int>(1));
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}



template <int dim>
std::vector<unsigned int>
FE_Q<dim>::lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe_data,
						  const unsigned int            degree)
{
  std::vector<unsigned int> renumber (fe_data.dofs_per_cell);
  
  const unsigned int n = degree+1;

  if (degree == 0)
    {
      Assert ((fe_data.dofs_per_vertex == 0) &&
	      ((fe_data.dofs_per_line == 0) || (dim == 1)) &&
	      ((fe_data.dofs_per_quad == 0) || (dim == 2)) &&
	      ((fe_data.dofs_per_hex == 0)  || (dim == 3)),
	      ExcInternalError());
      renumber[0] = 0;
    };

  if (degree > 0)
    {
      Assert (fe_data.dofs_per_vertex == 1, ExcInternalError());
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	{
	  unsigned int index = 0;
					   // Find indices of vertices.
					   // Unfortunately, somebody
					   // switched the upper corner
					   // points of a quad. The same
					   // person decided to find a very
					   // creative numbering of the
					   // vertices of a hexahedron.
					   // Therefore, this looks quite
					   // sophisticated.
					   //
					   // NB: This same person
					   // claims to have had good
					   // reasons then, but seems to
					   // have forgotten about
					   // them. At least, the
					   // numbering was discussed
					   // with the complaining
					   // person back then when all
					   // began :-)
	  switch (dim)
	    {
	      case 1:
	      {
		const unsigned int values[GeometryInfo<1>::vertices_per_cell]
		  = { 0, degree };
		index = values[i];
		break;
	      };
	     
	      case 2:
	      {
		const unsigned int values[GeometryInfo<2>::vertices_per_cell]
		  = { 0, degree, n*degree+degree, n*degree };
		index = values[i];
		break;
	      };
	     
	      case 3:
	      {
		const unsigned int values[GeometryInfo<3>::vertices_per_cell]
		  = { 0, degree,
		      n*n*degree + degree, n*n*degree,
		      n*degree, n*degree+degree,
		      n*n*degree + n*degree+degree, n*n*degree + n*degree};
		index = values[i];
		break;
	      };
	     
	      default:
		    Assert(false, ExcNotImplemented());
	    }

	  Assert (index<renumber.size(), ExcInternalError());
	  renumber[index] = i;
	}
    };
  
				   // for degree 2 and higher: Lines,
				   // quads, hexes etc also carry
				   // degrees of freedom
  if (degree > 1)
    {
      Assert (fe_data.dofs_per_line == degree-1, ExcInternalError());
      Assert ((fe_data.dofs_per_quad == (degree-1)*(degree-1)) ||
	      (dim < 2), ExcInternalError());
      Assert ((fe_data.dofs_per_hex == (degree-1)*(degree-1)*(degree-1)) ||
	      (dim < 3), ExcInternalError());
	    
      for (int i=0; i<static_cast<signed int>(GeometryInfo<dim>::lines_per_cell); ++i)
	{
	  unsigned int index = fe_data.first_line_index + i*fe_data.dofs_per_line;
	  unsigned int incr = 0;
	  unsigned int tensorstart = 0;
					   // This again looks quite
					   // strange because of the odd
					   // numbering scheme.
	  switch (i+100*dim)
	    {
					       // lines in x-direction
	      case 100:
	      case 200: case 202:
	      case 300: case 302: case 304: case 306:
		    incr = 1;
		    break;
						     // lines in y-direction
	      case 201: case 203:
	      case 308: case 309: case 310: case 311:
		    incr = n;
		    break;
						     // lines in z-direction
	      case 301: case 303: case 305: case 307:
		    incr = n*n;
		    break;
	      default:
		    Assert(false, ExcNotImplemented());
	    }
	  switch (i+100*dim)
	    {
					       // x=y=z=0
	      case 100:
	      case 200: case 203:
	      case 300: case 303: case 308:
		    tensorstart = 0;
		    break;
						     // x=1 y=z=0
	      case 201:
	      case 301: case 309:
		    tensorstart = degree;
		    break;
						     // y=1 x=z=0
	      case 202:
	      case 304: case 307:
		    tensorstart = n*degree;
		    break;
						     // x=z=1 y=0
	      case 310:
		    tensorstart = n*n*degree+degree;
		    break;
						     // z=1 x=y=0
	      case 302: case 311:
		    tensorstart = n*n*degree;
		    break;
						     // x=y=1 z=0
	      case 305:
		    tensorstart = n*degree+degree;
		    break;
						     // y=z=1 x=0
	      case 306:
		    tensorstart = n*n*n-n;
		    break;
	      default:
		    Assert(false, ExcNotImplemented());	      
	    }
	  
	  for (unsigned int jx = 1; jx<degree ;++jx)
	    {
	      unsigned int tensorindex = tensorstart + jx * incr;
	      Assert (tensorindex<renumber.size(), ExcInternalError());
	      renumber[tensorindex] = index++;
	    }
	}

      for (int i=0; i<static_cast<signed int>(GeometryInfo<dim>::quads_per_cell); ++i)
	{
	  unsigned int index = fe_data.first_quad_index+i*fe_data.dofs_per_quad;
	  unsigned int tensorstart = 0;
	  unsigned int incx = 0;
	  unsigned int incy = 0;
	  switch (i)
	    {
	      case 0:
		    tensorstart = 0; incx = 1;
		    if (dim==2)
		      incy = n;
		    else
		      incy = n*n;
		    break;
	      case 1:
		    tensorstart = n*degree; incx = 1; incy = n*n;
		    break;
	      case 2:
		    tensorstart = 0; incx = 1; incy = n;
		    break;
	      case 3:
		    tensorstart = degree; incx = n; incy = n*n;
		    break;
	      case 4:
		    tensorstart = n*n*degree; incx = 1; incy = n;
		    break;
	      case 5:
		    tensorstart = 0; incx = n; incy = n*n;
		    break;
	      default:
		    Assert(false, ExcNotImplemented());	      
	    }
	  
	  for (unsigned int jy = 1; jy<degree; jy++)
	    for (unsigned int jx = 1; jx<degree ;++jx)
	      {
		unsigned int tensorindex = tensorstart
					   + jx * incx + jy * incy;
		Assert (tensorindex<renumber.size(), ExcInternalError());
		renumber[tensorindex] = index++;
	      }
	}

      if (GeometryInfo<dim>::hexes_per_cell > 0)
	for (int i=0; i<static_cast<signed int>(GeometryInfo<dim>::hexes_per_cell); ++i)
	  {
	    unsigned int index = fe_data.first_hex_index;
	    
	    for (unsigned int jz = 1; jz<degree; jz++)
	      for (unsigned int jy = 1; jy<degree; jy++)
		for (unsigned int jx = 1; jx<degree; jx++)
		  {
		    const unsigned int tensorindex = jx + jy*n + jz*n*n;
		    Assert (tensorindex<renumber.size(), ExcInternalError());
		    renumber[tensorindex]=index++;
		  }  
	  } 
    }

  return renumber;
}



template <int dim>
std::vector<unsigned int>
FE_Q<dim>::face_lexicographic_to_hierarchic_numbering (const unsigned int degree)
{
  const FiniteElementData<dim-1> fe_data(FE_Q<dim-1>::get_dpo_vector(degree),1);
  return FE_Q<dim-1>::lexicographic_to_hierarchic_numbering (fe_data, degree); 
}


#if deal_II_dimension == 1

template <>
std::vector<unsigned int>
FE_Q<1>::face_lexicographic_to_hierarchic_numbering (const unsigned int)
{
  return std::vector<unsigned int>();
}

#endif

#if deal_II_dimension == 1

template <>
void
FE_Q<1>::initialize_constraints ()
{
				   // no constraints in 1d
}

#endif


#if deal_II_dimension == 2

template <>
void
FE_Q<2>::initialize_constraints ()
{
  const unsigned int dim = 2;

				   // restricted to each face, the
				   // traces of the shape functions is
				   // an element of P_{k} (in 2d), or
				   // Q_{k} (in 3d), where k is the
				   // degree of the element
				   //
				   // from this, we interpolate
				   // between mother and cell
				   // face. for the general case, this
				   // may be a little complicated if
				   // we don't use Lagrange
				   // interpolation polynomials, since
				   // then we can't just use point
				   // interpolation. what we do
				   // instead is to evaluate at a
				   // number of points and then invert
				   // the interpolation matrix. here,
				   // for the FE_Q elements, we
				   // actually do have Lagrange
				   // polynomials, but we still follow
				   // the general scheme since this
				   // code here is the master copy for
				   // what we use in other elements as
				   // well. however, there are places
				   // where we make use of the fact
				   // that we have Lagrange
				   // interpolation polynomials.

				   // mathematically speaking, the
				   // interpolation process works in
				   // the following way: on each
				   // subface, we want that finite
				   // element solututions from both
				   // sides coincide. i.e. if a and b
				   // are expansion coefficients for
				   // the shape functions from both
				   // sides, we seek a relation
				   // between x and y such that
				   //   sum_i a_i phi^c_i(x)
				   //   == sum_j b_j phi_j(x)
				   // for all points x on the
				   // interface. here, phi^c_i are the
				   // shape functions on the small
				   // cell on one side of the face,
				   // and phi_j those on the big cell
				   // on the other side. To get this
				   // relation, it suffices to look at
				   // a sufficient number of points
				   // for which this has to hold. if
				   // there are n functions, then we
				   // need n evaluation points, and we
				   // choose them equidistantly.
				   //
				   // what one then gets is a matrix
				   // system
				   //    a A  ==  b B
				   // where
				   //    A_ij = phi^c_i(x_j)
  				   //    B_ij = phi_i(x_j)
				   // and the relation we are looking for
				   // is
				   //    a = (A^T)^-1 B^T b
				   //
				   // below, we build up these
				   // matrices, but rather than
				   // transposing them after the
				   // fact, we do so while building
				   // them. A will be
				   // subface_interpolation, B will be
				   // face_interpolation. note that we
				   // build up these matrices for all
				   // faces at once, rather than
				   // considering them separately. the
				   // reason is that we finally will
				   // want to have them in this order
				   // anyway, as this is the format we
				   // need inside deal.II
  TensorProductPolynomials<dim-1>
    face_polynomials (Polynomials::LagrangeEquidistant::
		      generate_complete_basis (degree));
  Assert (face_polynomials.n() == this->dofs_per_face, ExcInternalError());
  
  const unsigned int n_small_functions = this->interface_constraints_size()[0];
  
  FullMatrix<double> face_interpolation (n_small_functions, this->dofs_per_face);
  FullMatrix<double> subface_interpolation (n_small_functions, n_small_functions);

  const std::vector<unsigned int>
    face_renumber_inverse (invert_numbering(face_renumber));
  
  for (unsigned int i=0; i<n_small_functions; ++i)
    {
				       // generate a quadrature point
				       // xi. it is actually not so
				       // important where this point
				       // lies, as long as we make
				       // sure that they are not
				       // equal. however, we will want
				       // them to be the (equidistant)
				       // Lagrange points, since then
				       // the subface_interpolation
				       // matrix has a most positive
				       // property: it is a
				       // permutation of the identity
				       // matrix. so create an
				       // equidistant mesh of points
				       // in the interior of the face
				       // (in 2d). for 3d, things are
				       // somewhat more convoluted as
				       // usual, since the new (child)
				       // shape functions are not only
				       // located in the interior of
				       // the face, but also on the
				       // edges, with the exception of
				       // the four vertices of the
				       // face. the function we call
				       // takes care of all this
      const Point<dim-1> p_face = generate_face_unit_point (i, n_small_functions,
							    int2type<dim-1>());

				       // evaluate the big face
				       // shape function at this
				       // point. note that the
				       // numbering of our shape
				       // functions is different
				       // from that of the
				       // polynomial, which orders
				       // them in the order of
				       // interpolation points.
				       //
				       // face_renumber_inverse will
				       // get us over this little
				       // conversion
      for (unsigned int j=0; j<this->dofs_per_face; ++j)
	{
	  face_interpolation(i,j)
	    = face_polynomials.compute_value(face_renumber_inverse[j], p_face);
					   // if the value is small up
					   // to round-off, then
					   // simply set it to zero to
					   // avoid unwanted fill-in
					   // of the constraint
					   // matrices (which would
					   // then increase the number
					   // of other DoFs a
					   // constrained DoF would
					   // couple to)
	  if (std::fabs(face_interpolation(i,j)) < 1e-14)
	    face_interpolation(i,j) = 0;
	}
	
				       // then evaluate all the
				       // small shape functions at
				       // this point.
      for (unsigned int j=0; j<n_small_functions; ++j)
	{
					   // first thing is to check
					   // which face the present
					   // point is on,
					   // i.e. whether it is left
					   // or right of the middle
					   // vertex in 2d, or
					   // something more complex
					   // in 3d (note that we
					   // might actually be
					   // sitting on top of the
					   // center vertex, or on on
					   // interface between
					   // children, but that
					   // doesn't really bother
					   // us: the shape functions
					   // associated with that
					   // have the same value
					   // whether we consider the
					   // left or the right
					   // subface, and all other
					   // shape functions should
					   // be zero there as well,
					   // so it doesn't really
					   // matter whether we
					   // account for this fact or
					   // not...)
	  const unsigned int subface
	    = GeometryInfo<dim-1>::child_cell_from_point (p_face);

					   // then check whether small
					   // shape function number j
					   // is nonzero on this
					   // face. as usual with our
					   // numbering of shape
					   // functions in constraint
					   // matrices, this is messy,
					   // so have a function that
					   // does this for us
					   //
					   // if not active, then the
					   // entry in the matrix will
					   // remain zero, and we
					   // simply go on with the
					   // next entry
	  if (! constraint_function_is_active_on_child (j, subface, *this))
	    continue;
	    
					   // otherwise: compute the
					   // coordinates of this
					   // evaluation point on
					   // the small face
	  const Point<dim-1> p_subface
	    = GeometryInfo<dim-1>::cell_to_child_coordinates (p_face, subface);
	
					   // then get the index of
					   // small shape function j
					   // on this subface. again,
					   // divert to a function
					   // that is specialized for
					   // this
	  const unsigned int local_j
	    = constraint_get_local_j (j, subface, *this);

					   // so evaluate this shape
					   // function there. now,
					   // since we have been
					   // careful with our choice
					   // of evaluation points,
					   // this is not actually
					   // necessary: the values of
					   // the small shape
					   // functions at these
					   // points should be either
					   // zero, and we can
					   // precompute which they
					   // are. However, we double
					   // check just to be sure we
					   // didn't do something
					   // stupid...
					   //
					   // (we could just set the
					   // evaluated value, but
					   // we'd end up with a lot
					   // of almost-zero entries,
					   // which will then carry
					   // over to the final
					   // result. this clutters up
					   // the constraint matrices,
					   // which we want to keep as
					   // small as possible.)
	  if (constraint_is_support_point (i, j, subface, *this))
	    subface_interpolation(i, j) = 1.;
	  else
	    subface_interpolation(i, j) = 0.;
	  Assert (std::fabs (subface_interpolation(i, j) -
			     face_polynomials.compute_value(local_j, p_subface))
		  < 1e-12,
		  ExcInternalError());
	}
    }

				   // what we now want to do is to
				   // compute
				   //   (subface_intp)^-1 face_intp
				   // which should give us the
				   // desired hanging node constraints.
				   // rather than actually doing this,
				   // we note that we have constructed
				   // subface_interpolation to be a
				   // permutation of the unit matrix.
				   // rather than doing a gauss jordan
				   // inversion, we note that the
				   // inverse is actually given by the
				   // transpose of the matrix. This has
				   // the additional benefit of being
				   // more stable and in particular of
				   // not adding almost-zeros
  this->interface_constraints
    .TableBase<2,double>::reinit (this->interface_constraints_size());
  subface_interpolation.Tmmult (this->interface_constraints,
				face_interpolation);

				   // in 3d we still have the
				   // constraint matrices, so make the
				   // check
  if (dim == 3)
    if (degree < Matrices::n_constraint_matrices+1)
      {
	FullMatrix<double> x;
	x.TableBase<2,double>::reinit (this->interface_constraints_size());
	x.fill (Matrices::constraint_matrices[degree-1]);

	for (unsigned int i=0; i<x.m(); ++i)
	  for (unsigned int j=0; j<x.n(); ++j)
	    Assert (std::fabs (x(i,j) - this->interface_constraints(i,j))
                    <
                    1e-14,
		    ExcInternalError());
      }    
}

#endif

#if deal_II_dimension == 3

template <>
void
FE_Q<3>::initialize_constraints ()
{
                                   // the algorithm for 2d is written
                                   // in a way so that it can be
                                   // extended to 3d as well. however,
                                   // the weird numbering convention
                                   // makes this really really hard,
                                   // so we abandoned this project at
                                   // one point. the plan is to change
                                   // the numbering convention for the
                                   // constraint matrices, and then
                                   // the approach for 2d will be
                                   // readily extendable to 3d as
                                   // well, but until this happens we
                                   // rather prefer to go back to the
                                   // precomputed matrices in 3d
  if (degree < Matrices::n_constraint_matrices+1)
    {
      this->interface_constraints
        .TableBase<2,double>::reinit (this->interface_constraints_size());
      this->interface_constraints.fill (Matrices::constraint_matrices[degree-1]);
    }
}

#endif


template <int dim>
void
FE_Q<dim>::initialize_embedding ()
{
				   // compute the interpolation
				   // matrices in much the same way as
				   // we do for the constraints. it's
				   // actually simpler here, since we
				   // don't have this weird
				   // renumbering stuff going on
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
					   // generate a point on
					   // the child cell and
					   // evaluate the shape
					   // functions there
	  const Point<dim> p_subcell = generate_unit_point (j, this->dofs_per_cell,
							    int2type<dim>());
	  const Point<dim> p_cell =
	    GeometryInfo<dim>::child_to_cell_coordinates (p_subcell, child);

	  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
	    {
	      const double
		cell_value    = polynomial_space.compute_value (i, p_cell),
		subcell_value = polynomial_space.compute_value (i, p_subcell);

					       // cut off values that
					       // are too small. note
					       // that we have here
					       // Lagrange
					       // interpolation
					       // functions, so they
					       // should be zero at
					       // almost all points,
					       // and one at the
					       // others, at least on
					       // the subcells. so set
					       // them to their exact
					       // values
                                               //
                                               // the actual cut-off
                                               // value is somewhat
                                               // fuzzy, but it works
                                               // for
                                               // 1e-14*degree*dim,
                                               // which is kind of
                                               // reasonable given
                                               // that we compute the
                                               // values of the
                                               // polynomials via an
                                               // degree-step
                                               // recursion and then
                                               // multiply the
                                               // 1d-values. this
                                               // gives us a linear
                                               // growth in
                                               // degree*dim, times a
                                               // small constant.
	      if (std::fabs(cell_value) < 2e-14*degree*dim)
		cell_interpolation(renumber[j], renumber[i]) = 0.;
	      else
		cell_interpolation(renumber[j], renumber[i]) = cell_value;

	      if (std::fabs(subcell_value) < 2e-14*degree*dim)
		subcell_interpolation(renumber[j], renumber[i]) = 0.;
	      else
		if (std::fabs(subcell_value-1) < 2e-14*degree*dim)
		  subcell_interpolation(renumber[j], renumber[i]) = 1.;
		else
						   // we have put our
						   // evaluation
						   // points onto the
						   // interpolation
						   // points, so we
						   // should either
						   // get zeros or
						   // ones!
		  Assert (false, ExcInternalError());
	    }
	}

				       // then compute the embedding
				       // matrix for this child and
				       // this coordinate
				       // direction. by the same trick
				       // as with the constraint
				       // matrices, don't compute the
				       // inverse of
				       // subcell_interpolation, but
				       // use the fact that we have
				       // put our interpolation points
				       // onto the interpolation
				       // points of the Lagrange
				       // polynomials used here. then,
				       // the subcell_interpolation
				       // matrix is just a permutation
				       // of the identity matrix and
				       // its inverse is also its
				       // transpose
      subcell_interpolation.Tmmult (this->prolongation[child],
                                    cell_interpolation);

					 // cut off very small values
					 // here
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
	for (unsigned int j=0; j<this->dofs_per_cell; ++j)
	  if (std::fabs(this->prolongation[child](i,j)) < 2e-14*degree*dim)
	    this->prolongation[child](i,j) = 0.;

				       // and make sure that the row
				       // sum is 1. this must be so
				       // since for this element, the
				       // shape functions add up to on
      for (unsigned int row=0; row<this->dofs_per_cell; ++row)
	{
	  double sum = 0;
	  for (unsigned int col=0; col<this->dofs_per_cell; ++col)
	    sum += this->prolongation[child](row,col);
	  Assert (std::fabs(sum-1.) < 2e-14*degree*dim,
		  ExcInternalError());
	};
    }
}



template <int dim>
void
FE_Q<dim>::initialize_restriction ()
{
                                   // for these Lagrange interpolation
                                   // polynomials, construction of the
                                   // restriction matrices is
                                   // relatively simple. the reason is
                                   // that the interpolation points on
                                   // the mother cell are always also
                                   // interpolation points for some
                                   // shape function on one or the
                                   // other child, because we have
                                   // chosen equidistant Lagrange
                                   // interpolation points for the
                                   // polynomials
                                   //
                                   // so the only thing we have to
                                   // find out is: for each shape
                                   // function on the mother cell,
                                   // which is the child cell
                                   // (possibly more than one) on
                                   // which it is located, and which
                                   // is the corresponding shape
                                   // function there. rather than
                                   // doing it for the shape functions
                                   // on the mother cell, we take the
                                   // interpolation points there are
                                   // also search which shape function
                                   // corresponds to it (too lazy to
                                   // do this mapping by hand)
                                   //
                                   // note that the interpolation
                                   // point of a shape function can be
                                   // on the boundary between
                                   // subcells. in that case,
                                   // restriction from children to
                                   // mother may produce two or more
                                   // entries for a dof on the mother
                                   // cell. however, this doesn't
                                   // hurt: since the element is
                                   // continuous, the contribution
                                   // from each child should yield the
                                   // same result, and since the
                                   // element is non-additive we just
                                   // overwrite one value (compute one
                                   // one child) by the same value
                                   // (compute on a later child), so
                                   // we don't have to care about this
  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    this->restriction[c].reinit (this->dofs_per_cell, this->dofs_per_cell);
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    {
      const Point<dim> p_cell = generate_unit_point (i, this->dofs_per_cell,
                                                     int2type<dim>());
      unsigned int mother_dof = 0;
      for (; mother_dof<this->dofs_per_cell; ++mother_dof)
        {
          const double val
            = polynomial_space.compute_value(renumber_inverse[mother_dof],
                                             p_cell);
          if (std::fabs (val-1.) < 2e-14*degree*dim)
                                             // ok, this is the right
                                             // dof
            break;
          else
                                             // make sure that all
                                             // other shape functions
                                             // are zero there
            Assert (std::fabs(val) < 2e-14*degree*dim,
                    ExcInternalError());
        }
                                       // check also the shape
                                       // functions after tat
      for (unsigned int j=mother_dof+1; j<this->dofs_per_cell; ++j)
        Assert (std::fabs (polynomial_space.compute_value(renumber_inverse[j],
                                                          p_cell))
                < 2e-14*degree*dim,
                ExcInternalError());

                                       // then find the children on
                                       // which the interpolation
                                       // point is located
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell;
           ++child)
        {
                                           // first initialize this
                                           // column of the matrix
          for (unsigned int j=0; j<this->dofs_per_cell; ++j)
            this->restriction[child](mother_dof, j) = 0.;

                                           // then check whether this
                                           // interpolation point is
                                           // inside this child cell
          const Point<dim> p_subcell
            = GeometryInfo<dim>::cell_to_child_coordinates (p_cell, child);
          if (GeometryInfo<dim>::is_inside_unit_cell (p_subcell))
            {
                                               // find the one child
                                               // shape function
                                               // corresponding to
                                               // this point. do it in
                                               // the same way as
                                               // above
              unsigned int child_dof = 0;
              for (; child_dof<this->dofs_per_cell; ++child_dof)
                {
                  const double val
                    = polynomial_space.compute_value(renumber_inverse[child_dof],
                                                     p_subcell);
                  if (std::fabs (val-1.) < 2e-14*degree*dim)
                    break;
                  else
                    Assert (std::fabs(val) < 2e-14*degree*dim,
                            ExcInternalError());
                }
              for (unsigned int j=child_dof+1; j<this->dofs_per_cell; ++j)
                Assert (std::fabs (polynomial_space.compute_value(renumber_inverse[j],
                                                                  p_subcell))
                        < 2e-14*degree*dim,
                        ExcInternalError());

                                               // so now that we have
                                               // it, set the
                                               // corresponding value
                                               // in the matrix
              this->restriction[child](mother_dof, child_dof) = 1.;
            }
        }
    }
}



template <int dim>
UpdateFlags
FE_Q<dim>::update_once (const UpdateFlags flags) const
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
FE_Q<dim>::update_each (const UpdateFlags flags) const
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
FE_Q<dim>::get_data (const UpdateFlags      update_flags,
		     const Mapping<dim>    &mapping,
		     const Quadrature<dim> &quadrature) const
{
				   // generate a new data object and
				   // initialize some fields
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

				   // some scratch arrays
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
	
	if (flags & update_values)
	  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	    data->shape_values[renumber[k]][i] = values[k];
	
	if (flags & update_gradients)
	  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	    data->shape_gradients[renumber[k]][i] = grads[k];
      }
  return data;
}




//----------------------------------------------------------------------
// Fill data of FEValues
//----------------------------------------------------------------------

template <int dim>
void
FE_Q<dim>::fill_fe_values (const Mapping<dim>                   &mapping,
			   const typename DoFHandler<dim>::cell_iterator &cell,
			   const Quadrature<dim>                &quadrature,
			   typename Mapping<dim>::InternalDataBase &mapping_data,
			   typename Mapping<dim>::InternalDataBase &fedata,
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
	for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
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
FE_Q<dim>::fill_fe_face_values (const Mapping<dim>                   &mapping,
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
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
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
FE_Q<dim>::fill_fe_subface_values (const Mapping<dim>                   &mapping,
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
  const unsigned int offset = ((face * GeometryInfo<dim>::subfaces_per_face + subface)
                               * quadrature.n_quadrature_points);

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
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
FE_Q<dim>::n_base_elements () const
{
  return 1;
}



template <int dim>
const FiniteElement<dim> &
FE_Q<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <int dim>
unsigned int
FE_Q<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}



template <int dim>
bool
FE_Q<dim>::has_support_on_face (const unsigned int shape_index,
				const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));


				   // in 1d, things are simple. since
				   // there is only one degree of
				   // freedom per vertex in this
				   // class, the first is on vertex 0
				   // (==face 0 in some sense), the
				   // second on face 1:
  if (dim==1)
    return (((shape_index == 0) && (face_index == 0)) ||
	    ((shape_index == 1) && (face_index == 1)));
  else
				     // more dimensions
    {
                                       // first, special-case interior
                                       // shape functions, since they
                                       // have no support no-where on
                                       // the boundary
      if (((dim==2) && (shape_index>=this->first_quad_index))
          ||
          ((dim==3) && (shape_index>=this->first_hex_index)))
        return false;
                                       
                                       // let's see whether this is a
                                       // vertex
      if (shape_index < this->first_line_index) 
        {
                                           // for Q elements, there is
                                           // one dof per vertex, so
                                           // shape_index==vertex_number. check
                                           // whether this vertex is
                                           // on the given face. thus,
                                           // for each face, give a
                                           // list of vertices
          const unsigned int vertex_no = shape_index;
          Assert (vertex_no < GeometryInfo<dim>::vertices_per_cell,
                  ExcInternalError());
          switch (dim)
            {
              case 2:
              {
                static const unsigned int face_vertices[4][2] =
                  { {0,1},{1,2},{2,3},{0,3} };
                return ((face_vertices[face_index][0] == vertex_no)
                        ||
                        (face_vertices[face_index][1] == vertex_no));
              };

              case 3:
              {
                static const unsigned int face_vertices[6][4] =
                  { {0,1,2,3},{4,5,6,7},{0,1,5,4},
                    {1,5,6,2},{3,2,6,7},{0,4,7,3} };
                return ((face_vertices[face_index][0] == vertex_no)||
                        (face_vertices[face_index][1] == vertex_no)||
                        (face_vertices[face_index][2] == vertex_no)||
                        (face_vertices[face_index][3] == vertex_no));
              };

              default:
                    Assert (false, ExcNotImplemented());
            };
        }
      else if (shape_index < this->first_quad_index)
                                         // ok, dof is on a line
        {
          const unsigned int line_index
            = (shape_index - this->first_line_index) / this->dofs_per_line;
          Assert (line_index < GeometryInfo<dim>::lines_per_cell,
                  ExcInternalError());

                                           // in 2d, the line is the
                                           // face, so get the line
                                           // index
          if (dim == 2)
            return (line_index == face_index);
          else if (dim == 3)
            {
                                               // see whether the
                                               // given line is on the
                                               // given face. use
                                               // table technique
                                               // again
              static const unsigned int face_lines[6][4] =
                { {0,1,2,3},{4,5,6,7},{0,8,9,4},
                  {1,9,5,10},{2,10,6,11},{3,8,7,11} };
              return ((face_lines[face_index][0] == line_index)||
                      (face_lines[face_index][1] == line_index)||
                      (face_lines[face_index][2] == line_index)||
                      (face_lines[face_index][3] == line_index));
            }
          else
            Assert (false, ExcNotImplemented());
        }
      else if (shape_index < this->first_hex_index)
                                         // dof is on a quad
        {
          const unsigned int quad_index 
            = (shape_index - this->first_quad_index) / this->dofs_per_quad;
          Assert (static_cast<signed int>(quad_index) <
                  static_cast<signed int>(GeometryInfo<dim>::quads_per_cell),
                  ExcInternalError());
          
                                           // in 2d, cell bubble are
                                           // zero on all faces. but
                                           // we have treated this
                                           // case above already
          Assert (dim != 2, ExcInternalError());

                                           // in 3d,
                                           // quad_index=face_index
          if (dim == 3)
            return (quad_index == face_index);
          else
            Assert (false, ExcNotImplemented());
        }
      else
                                         // dof on hex
        {
                                           // can only happen in 3d,
                                           // but this case has
                                           // already been covered
                                           // above
          Assert (false, ExcNotImplemented());
          return false;
        };
    };

                                   // we should not have gotten here
  Assert (false, ExcInternalError());
  return false;
}



template <int dim>
unsigned int
FE_Q<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template <int dim>
unsigned int
FE_Q<dim>::get_degree () const
{
  return degree;
}



template class FE_Q<deal_II_dimension>;
