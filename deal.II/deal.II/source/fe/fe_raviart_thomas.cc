//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------

#include <base/quadrature.h>
#include <base/table.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_values.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


// namespace for some functions that are used in this file. they are
// specific to numbering conventions used for the FE_RT element, and
// are thus not very interesting to the outside world
namespace 
{
  				   // auxiliary type to allow for some
				   // kind of explicit template
				   // specialization of the following
				   // functions
  template <int dim> struct int2type {};

  
				   // generate the j-th out of a total
				   // of N points on the unit square
				   // in 2d. N needs not be a square
				   // number, but must be the product
				   // of two integers
				   //
				   // there is one complication: we
				   // want to generate interpolation
				   // points on the unit square for
				   // the shape functions for this
				   // element, but for that we need to
				   // make sure that these
				   // interpolation points make the
				   // resulting matrix rows linearly
				   // independent. this is a problem
				   // since we have anisotropic
				   // polynomials, so for example for
				   // the lowest order elements, we
				   // have as polynomials in for the
				   // x-component of the shape
				   // functions only "x" and "1-x",
				   // i.e. no y-dependence. if we
				   // select as interpolation points
				   // the points (.5,0) and (.5,1),
				   // we're hosed!
				   //
				   // thus, the third parameter gives
				   // the coordinate direction in
				   // which the polynomial degree is
				   // highest. we use this to select
				   // interpolation points primarily
				   // in this direction then
  Point<2> generate_unit_point (const unsigned int j,
				const unsigned int N,
				const unsigned int d,
				const int2type<2>  &)
  {
    Assert (d<2, ExcInternalError());
    
				     // factorize N int N1*N2. note
				     // that we always have N1<=N2,
				     // since the square root is
				     // rounded down
    const unsigned int N1 = static_cast<unsigned int>(std::sqrt(1.*N));
    const unsigned int N2 = N/N1;
    Assert (N1*N2 == N, ExcInternalError());

    const unsigned int Nx = (d==0 ? N2 : N1),
		       Ny = (d==1 ? N2 : N1);
    
    return Point<2> (Nx == 1 ? .5 : 1.*(j%Nx)/(Nx-1),
		     Ny == 1 ? .5 : 1.*(j/Nx)/(Ny-1));
  }
  

				   // generate the j-th out of a total
				   // of N points on the unit cube
				   // in 3d. N needs not be a cube
				   // number, but must be the product
				   // of three integers
				   //
				   // the same applies as above for
				   // the meaning of the parameter "d"
  Point<3> generate_unit_point (const unsigned int /*j*/,
				const unsigned int N,
				const unsigned int d,
				const int2type<3>  &)
  {
    Assert (d<3, ExcInternalError());

    const unsigned int N1 = static_cast<unsigned int>(std::pow(1.*N, 1./3.));
    const unsigned int N2 = static_cast<unsigned int>(std::sqrt(1.*N/N1));
    const unsigned int N3 = N/(N1*N2);
    Assert (N1*N2*N3 == N, ExcInternalError());

    Assert (false, ExcNotImplemented());

    return Point<3> ();
  }
  
}



template <int dim>
FE_RaviartThomas<dim>::FE_RaviartThomas (const unsigned int degree)
		:
		FiniteElement<dim> (FiniteElementData<dim>(get_dpo_vector(degree),
							   dim, degree+1),
				    get_ria_vector (degree),
				    std::vector<std::vector<bool> >(FiniteElementData<dim>(get_dpo_vector(degree),dim,degree+1).dofs_per_cell,
								    std::vector<bool>(dim,true))),
		degree(degree),
                polynomials (create_polynomials(degree)),
                renumber (compute_renumber(degree))
{
  Assert (dim >= 2, ExcNotUsefulInThisDimension());

                                   // check formula (III.3.22) in the
                                   // book by Brezzi & Fortin about
                                   // the number of degrees of freedom
                                   // per cell
  Assert (((dim==2) &&
           (this->dofs_per_cell == 2*(degree+1)*(degree+2)))
          ||
          ((dim==3) &&
           (this->dofs_per_cell == 3*(degree+1)*(degree+1)*(degree+2))),
          ExcInternalError());
  Assert (renumber.size() == this->dofs_per_cell,
          ExcInternalError());

				   // initialize the various matrices
  initialize_constraints ();
  initialize_embedding ();
  initialize_restriction ();

				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();

                                   // then make
                                   // system_to_component_table
                                   // invalid, since this has no
                                   // meaning for the present element
  std::vector<std::pair<unsigned,unsigned> > tmp1, tmp2;
  this->system_to_component_table.swap (tmp1);
  this->face_system_to_component_table.swap (tmp2);
}



template <int dim>
std::string
FE_RaviartThomas<dim>::get_name () const
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
  
  namebuf << "FE_RaviartThomas<" << dim << ">(" << degree << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_RaviartThomas<dim>::clone() const
{
  return new FE_RaviartThomas<dim>(degree);
}


template <int dim>
double
FE_RaviartThomas<dim>::shape_value_component (const unsigned int i,
                                              const Point<dim>    &p,
                                              const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));

                                   // check whether this shape
                                   // function has a contribution in
                                   // this component at all, and if so
                                   // delegate to the respective
                                   // polynomial
  if (component == renumber[i].first)
    return polynomials[component].compute_value(renumber[i].second, p);
  else
    return 0;
}



template <int dim>
Tensor<1,dim>
FE_RaviartThomas<dim>::shape_grad_component (const unsigned int  i,
                                             const Point<dim>   &p,
                                             const unsigned int  component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));

                                   // check whether this shape
                                   // function has a contribution in
                                   // this component at all, and if so
                                   // delegate to the respective
                                   // polynomial
  if (component == renumber[i].first)
    return polynomials[component].compute_grad(renumber[i].second, p);
  else
    return Tensor<1,dim>();
}



template <int dim>
Tensor<2,dim>
FE_RaviartThomas<dim>::shape_grad_grad_component (const unsigned int i,
                                                  const Point<dim>  &p,
                                                  const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));

                                   // check whether this shape
                                   // function has a contribution in
                                   // this component at all, and if so
                                   // delegate to the respective
                                   // polynomial
  if (component == renumber[i].first)
    return polynomials[component].compute_grad_grad(renumber[i].second, p);
  else
    return Tensor<2,dim>();
}



#if deal_II_dimension == 1

template <>
void
FE_RaviartThomas<1>::
get_interpolation_matrix (const FiniteElementBase<1> &,
			  FullMatrix<double>         &) const
{
  Assert (false, ExcNotUsefulInThisDimension());
}

#endif


template <int dim>
void
FE_RaviartThomas<dim>::
get_interpolation_matrix (const FiniteElementBase<dim> &x_source_fe,
			  FullMatrix<double>           &interpolation_matrix) const
{
				   // this is only implemented, if the
				   // source FE is also a
				   // Raviart-Thomas element,
				   // otherwise throw an exception, as
				   // the documentation says
  AssertThrow ((x_source_fe.get_name().find ("FE_RaviartThomas<") == 0)
               ||
               (dynamic_cast<const FE_RaviartThomas<dim>*>(&x_source_fe) != 0),
               typename FiniteElementBase<dim>::
               ExcInterpolationNotImplemented());
  
				   // ok, source is a RT element, so
				   // we will be able to do the work
  const FE_RaviartThomas<dim> &source_fe
    = dynamic_cast<const FE_RaviartThomas<dim>&>(x_source_fe);

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
  const unsigned int dofs_per_coordinate = this->dofs_per_cell/dim;
  Assert (dofs_per_coordinate*dim == this->dofs_per_cell,
	  ExcInternalError());
  for (unsigned int d=0; d<dim; ++d)
    Assert (polynomials[d].n() == dofs_per_coordinate, ExcInternalError());

  const unsigned int source_dofs_per_coordinate = source_fe.dofs_per_cell/dim;
  Assert (source_dofs_per_coordinate*dim == source_fe.dofs_per_cell,
	  ExcInternalError());
  for (unsigned int d=0; d<dim; ++d)
    Assert (source_fe.polynomials[d].n() == source_dofs_per_coordinate, ExcInternalError());
  
  FullMatrix<double> cell_interpolation (dofs_per_coordinate,
					 dofs_per_coordinate);
  FullMatrix<double> source_interpolation (dofs_per_coordinate,
					   source_dofs_per_coordinate);
  FullMatrix<double> tmp (dofs_per_coordinate,
			  source_dofs_per_coordinate);
  for (unsigned int d=0; d<dim; ++d)
    {
      for (unsigned int j=0; j<dofs_per_coordinate; ++j)
	{
					   // generate a point on this
					   // cell and evaluate the
					   // shape functions there
					   //
					   // see the comment for that
					   // function to see why the
					   // third parameter is
					   // necessary
	  const Point<dim> p = generate_unit_point (j, dofs_per_coordinate,
						    d, int2type<dim>());
	  for (unsigned int i=0; i<dofs_per_coordinate; ++i)
	    cell_interpolation(j,i) = polynomials[d].compute_value (i, p);

	  for (unsigned int i=0; i<source_dofs_per_coordinate; ++i)
	    source_interpolation(j,i) = source_fe.polynomials[d].compute_value (i, p);
	}

				       // then compute the
				       // interpolation matrix matrix
				       // for this coordinate
				       // direction
      cell_interpolation.gauss_jordan ();
      cell_interpolation.mmult (tmp, source_interpolation);

				       // finally transfer the
				       // results for this
				       // coordinate into the matrix
				       // corresponding to the
				       // entire space on this
				       // cell. cut off very small
				       // values here
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
	if (renumber[i].first == d)
	  for (unsigned int j=0; j<source_fe.dofs_per_cell; ++j)
	    if (source_fe.renumber[j].first == d)
	      if (std::fabs(tmp(renumber[i].second,
				source_fe.renumber[j].second)) > 1e-15)
		interpolation_matrix(i,j) = tmp(renumber[i].second,
						source_fe.renumber[j].second);
    }

				   // if this were a Lagrange
				   // interpolation element, we could
				   // make sure that the row sum of
				   // each of the matrices is 1 at
				   // this point. note that this won't
				   // work here, since we are working
				   // with hierarchical elements for
				   // which the shape functions don't
				   // sum up to 1
				   //
				   // however, we can make sure that
				   // only components couple that have
				   // the same vector component
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    for (unsigned int j=0; j<source_fe.dofs_per_cell; ++j)
      Assert ((interpolation_matrix(i,j) == 0.) ||
	      (renumber[i].first == source_fe.renumber[j].first),
	      ExcInternalError());
}



//----------------------------------------------------------------------
// Auxiliary and internal functions
//----------------------------------------------------------------------




#if deal_II_dimension == 1

template <>
void
FE_RaviartThomas<1>::initialize_constraints ()
{
  Assert (false, ExcNotUsefulInThisDimension());
}

#endif

#if deal_II_dimension == 2

template <>
void
FE_RaviartThomas<2>::initialize_constraints ()
{
  const unsigned int dim = 2;
  
  this->interface_constraints.
    TableBase<2,double>::reinit (this->interface_constraints_size());

				   // this case is too easy, so
				   // special case it
  if (degree == 0)
    {
      this->interface_constraints(0,0) = this->interface_constraints(1,0) = .5;
      return;
    }

				   // for higher orders of the
				   // Raviart-Thomas element:
    
				   // restricted to each face, the
				   // normal component of the shape
				   // functions is an element of P_{k}
				   // (in 2d), or Q_{k} (in 3d), where
				   // k is the degree of the element
				   //
				   // from this, we interpolate
				   // between mother and cell
				   // face. this is slightly
				   // complicated by the fact that we
				   // don't use Lagrange interpolation
				   // polynomials, but rather
				   // hierarchical polynomials, so we
				   // can't just use point
				   // interpolation. what we do
				   // instead is to evaluate at a
				   // number of points and then invert
				   // the interpolation matrix

				   // mathematically speaking, this
				   // works in the following way: on
				   // each subface, we want that
				   // finite element solututions from
				   // both sides coincide. i.e. if a
				   // and b are expansion coefficients
				   // for the shape functions from
				   // both sides, we seek a relation
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
  const std::vector<Polynomials::Polynomial<double> >
    face_polynomials (Polynomials::Hierarchical::
		      generate_complete_basis (degree));
  Assert (face_polynomials.size() == this->dofs_per_face, ExcInternalError());
  
  FullMatrix<double> face_interpolation (2*this->dofs_per_face, this->dofs_per_face);
  FullMatrix<double> subface_interpolation (2*this->dofs_per_face, 2*this->dofs_per_face);
  
				   // generate the matrix for the
				   // evaluation points on the big
				   // face, and the corresponding
				   // points in the coordinate system
				   // of the small face. order the
				   // shape functions in the same way
				   // we want to have them in the
				   // final matrix. extend shape
				   // functions on the small faces by
				   // zero to the other face on which
				   // they are not defined (we do this
				   // by simply not considering these
				   // entries in the matrix)
				   //
				   // note the agreeable fact that for
				   // this element, all the shape
				   // functions we presently care for
				   // are face-based (i.e. not vertex
				   // shape functions); thus, for this
				   // element, we can skip the
				   // annoying index shifting for the
				   // constraints matrix due to its
				   // weird format
  for (unsigned int subface=0; subface<GeometryInfo<dim>::subfaces_per_face; ++subface)
    for (unsigned int i=0; i<this->dofs_per_face; ++i)
      {
	const double p_face (1.*i/degree/2 + (subface == 0 ? 0. : .5));
	const double p_subface (1.*i/degree);

	for (unsigned int j=0; j<this->dofs_per_face; ++j)
	  {
	    face_interpolation(subface*this->dofs_per_face+i,
			       j)
	      = face_polynomials[j].value(p_face);
	    subface_interpolation(subface*this->dofs_per_face+i,
				  subface*this->dofs_per_face+j)
	      = face_polynomials[j].value(p_subface);
	  }
      }
      
  subface_interpolation.gauss_jordan ();
  subface_interpolation.mmult (this->interface_constraints,
			       face_interpolation);
    
				   // there is one additional thing to
				   // be considered: since the shape
				   // functions on the real cell
				   // contain the Jacobian (actually,
				   // the determinant of the inverse),
				   // there is an additional factor of
				   // 2 when going from the big to the
				   // small cell:
  this->interface_constraints *= 1./2;

				   // finally: constraints become
				   // really messy if the matrix in
				   // question has some entries that
				   // are almost zero, but not
				   // quite. this will happen in the
				   // above procedure due to
				   // round-off. let us simply delete
				   // these entries
  for (unsigned int i=0; i<this->interface_constraints.m(); ++i)
    for (unsigned int j=0; j<this->interface_constraints.n(); ++j)
      if (std::fabs(this->interface_constraints(i,j)) < 1e-14)
	this->interface_constraints(i,j) = 0.;
}

#endif

#if deal_II_dimension == 3

template <>
void
FE_RaviartThomas<3>::initialize_constraints ()
{
  Assert (false, ExcNotImplemented());
}

#endif


#if deal_II_dimension == 1

template <>
void
FE_RaviartThomas<1>::initialize_embedding ()
{
  Assert (false, ExcNotUsefulInThisDimension());
}

#endif


template <int dim>
void
FE_RaviartThomas<dim>::initialize_embedding ()
{
				   // compute the interpolation
				   // matrices in much the same way as
				   // we do for the constraints. it's
				   // actually simpler here, since we
				   // don't have this weird
				   // renumbering stuff going on
				   //
				   // it is, however, slightly
				   // complicated by the fact that we
				   // have vector-valued elements
				   // here, so we do all the stuff for
				   // the degrees of freedom
				   // corresponding to each coordinate
				   // direction separately
  const unsigned int dofs_per_coordinate = this->dofs_per_cell/dim;
  Assert (dofs_per_coordinate*dim == this->dofs_per_cell,
	  ExcInternalError());
  for (unsigned int d=0; d<dim; ++d)
    Assert (polynomials[d].n() == dofs_per_coordinate, ExcInternalError());
  
  FullMatrix<double> cell_interpolation (dofs_per_coordinate,
					 dofs_per_coordinate);
  FullMatrix<double> subcell_interpolation (dofs_per_coordinate,
					    dofs_per_coordinate);
  FullMatrix<double> tmp (dofs_per_coordinate,
			  dofs_per_coordinate);
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    this->prolongation[child].reinit (this->dofs_per_cell,
				      this->dofs_per_cell);
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    for (unsigned int d=0; d<dim; ++d)
      {
	for (unsigned int j=0; j<dofs_per_coordinate; ++j)
	  {
					     // generate a point on
					     // the child cell and
					     // evaluate the shape
					     // functions there
					     //
					     // see the comment for
					     // that function to see
					     // why the third
					     // parameter is necessary
	    const Point<dim> p_subcell = generate_unit_point (j, dofs_per_coordinate,
							      d, int2type<dim>());
	    const Point<dim> p_cell =
	      GeometryInfo<dim>::child_to_cell_coordinates (p_subcell, child);

	    for (unsigned int i=0; i<dofs_per_coordinate; ++i)
	      {
		cell_interpolation(j,i) = polynomials[d].compute_value (i, p_cell);
		subcell_interpolation(j,i) = polynomials[d].compute_value (i, p_subcell);
	      }
	  }

					 // then compute the embedding
					 // matrix for this child and
					 // this coordinate direction
	subcell_interpolation.gauss_jordan ();
	subcell_interpolation.mmult (tmp, cell_interpolation);

					 // finally transfer the
					 // results for this
					 // coordinate into the matrix
					 // corresponding to the
					 // entire space on this
					 // cell. cut off very small
					 // values here
	for (unsigned int i=0; i<this->dofs_per_cell; ++i)
	  if (renumber[i].first == d)
	    for (unsigned int j=0; j<this->dofs_per_cell; ++j)
	      if (renumber[j].first == d)
		if (std::fabs(tmp(renumber[i].second,renumber[j].second)) > 1e-15)
		  this->prolongation[child](i,j) = tmp(renumber[i].second,
						       renumber[j].second);
      }

				   // if this were a Lagrange
				   // interpolation element, we could
				   // make sure that the row sum of
				   // each of the matrices is 1 at
				   // this point. note that this won't
				   // work here, since we are working
				   // with hierarchical elements for
				   // which the shape functions don't
				   // sum up to 1
				   //
				   // however, we can make sure that
				   // only components couple that have
				   // the same vector component
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
	Assert ((this->prolongation[child](i,j) == 0.) ||
		(renumber[i].first == renumber[j].first),
		ExcInternalError());
      
  
				   // there is one additional thing to
				   // be considered: since the shape
				   // functions on the real cell
				   // contain the Jacobian (actually,
				   // the determinant of the inverse),
				   // there is an additional factor of
				   // 2 when going from the big to the
				   // small cell:
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    this->prolongation[child] *= 1./2;
}


#if deal_II_dimension == 1

template <>
void
FE_RaviartThomas<1>::initialize_restriction ()
{}

#endif


#if deal_II_dimension == 2

template <>
void
FE_RaviartThomas<2>::initialize_restriction ()
{
  const unsigned int dim = 2;
  switch (degree)
    {
      case 0:
      {
					 // this is a strange element,
					 // since it is both additive
					 // and then it is also
					 // not. ideally, we would
					 // like to have the value of
					 // the shape function on the
					 // coarse line to be the mean
					 // value of that on the two
					 // child ones. thus, one
					 // should make it
					 // additive. however,
					 // additivity only works if
					 // an element does not have
					 // any continuity
					 // requirements, since
					 // otherwise degrees of
					 // freedom are shared between
					 // adjacent elements, and
					 // when we make the element
					 // additive, that would mean
					 // that we end up adding up
					 // contributions not only
					 // from the child cells of
					 // this cell, but also from
					 // the child cells of the
					 // neighbor, and since we
					 // cannot know whether there
					 // even exists a neighbor we
					 // cannot simply make the
					 // element additive.
					 //
					 // so, until someone comes
					 // along with a better
					 // alternative, we do the
					 // following: make the
					 // element non-additive, and
					 // simply pick the value of
					 // one of the child lines for
					 // the value of the mother
					 // line (note that we have to
					 // multiply by two, since the
					 // shape functions scale with
					 // the inverse Jacobian). we
					 // thus throw away the
					 // information of one of the
					 // child lines, but there
					 // seems to be no other way
					 // than that...
					 //
					 // note: to make things
					 // consistent, and
					 // restriction independent of
					 // the order in which we
					 // travel across the cells of
					 // the coarse grid, we have
					 // to make sure that we take
					 // the same small line when
					 // visiting its two
					 // neighbors, to get the
					 // value for the mother
					 // line. we take the first
					 // line always, in the
					 // canonical direction of
					 // lines
	for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	  this->restriction[c].reinit (this->dofs_per_cell,
				       this->dofs_per_cell);
              
	this->restriction[0](0,0) = 2.;
	this->restriction[1](1,1) = 2.;
	this->restriction[3](2,2) = 2.;
	this->restriction[0](3,3) = 2.;

	break;
      };


      case 1:
      {
	for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	  this->restriction[c].reinit (this->dofs_per_cell,
				       this->dofs_per_cell);

					 // first set the corner
					 // nodes. note that they are
					 // non-additive
	this->restriction[0](0,0) = 2.;
	this->restriction[0](6,6) = 2.;

	this->restriction[1](1,1) = 2.;
	this->restriction[1](2,2) = 2.;

	this->restriction[2](3,3) = 2.;
	this->restriction[2](5,5) = 2.;

	this->restriction[3](4,4) = 2.;
	this->restriction[3](7,7) = 2.;

					 // then also set the bubble
					 // nodes. they _are_
					 // additive. to understand
					 // what's going on, recall
					 // that the bubble shape
					 // functions have value -1
					 // (!) at the center point,
					 // by construction of the
					 // polynomials, and that the
					 // corner nodes have values
					 // 1/2 there since they are
					 // just the linears, and not
					 // some interpolating
					 // polynomial
					 //
					 // (actually, the
					 // additive/non-additive
					 // business shouldn't make
					 // that much of a difference:
					 // node 4 on cell 0 and node
					 // 0 on cell 3 must have the
					 // same value, since normal
					 // components are
					 // continuous. so we could
					 // pick either and make these
					 // shape functions
					 // non-additive as well. we
					 // choose to take the mean
					 // value, which should be the
					 // same as either value, and
					 // make the shape function
					 // additive)
	this->restriction[0](10,0) = 1.;
	this->restriction[0](10,4) = -1.;
	this->restriction[3](10,0) = -1.;
	this->restriction[3](10,4) = 1.;

	this->restriction[1](11,1) = 1.;
	this->restriction[1](11,5) = -1.;
	this->restriction[2](11,1) = -1.;
	this->restriction[2](11,5) = 1.;
	
	this->restriction[0](8,6) = 1.;
	this->restriction[0](8,2) = -1.;
	this->restriction[1](8,6) = -1.;
	this->restriction[1](8,2) = 1.;
	
	this->restriction[3](9,7) = 1.;
	this->restriction[3](9,3) = -1.;
	this->restriction[2](9,7) = -1.;
	this->restriction[2](9,3) = 1.;
	
	break;
      };
	
					// in case we don't have the
					// matrices (yet), leave them
					// empty. this does not
					// prevent the use of this FE,
					// but will prevent the use of
					// these matrices
    };
}

#endif

#if deal_II_dimension == 3

template <>
void
FE_RaviartThomas<3>::initialize_restriction ()
{
  Assert (false, ExcNotImplemented());
}

#endif


template <int dim>
void
FE_RaviartThomas<dim>::initialize_unit_support_points ()
{
  this->unit_support_points.resize (this->dofs_per_cell);
  switch (dim) 
    {
      case 2:
      {
        Assert (degree+1 == this->dofs_per_face, ExcInternalError());
        
                                         // associate support points
                                         // with mid-face points if a
                                         // shape function has a
                                         // non-zero normal component
                                         // there, otherwise with the
                                         // cell center. the reason
                                         // for this non-unique
                                         // support point is that we
                                         // use hierarchical shape
                                         // functions, rather than
                                         // Lagrange functions, for
                                         // which we get into the same
                                         // trouble as in the
                                         // FE_Q_Hierarchical element;
                                         // see the respective
                                         // function there

                                         // start with the face shape
                                         // functions
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[0*this->dofs_per_face+i] = Point<dim>(.5, .0);
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[1*this->dofs_per_face+i] = Point<dim>(1., .5);
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[2*this->dofs_per_face+i] = Point<dim>(.5, 1.);
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[3*this->dofs_per_face+i] = Point<dim>(.0, .5);

                                         // associate the rest with
                                         // the cell center
        for (unsigned int i=4*this->dofs_per_face; i<this->dofs_per_cell; ++i)
          this->unit_support_points[i] = Point<dim>(.5, .5);
        
        break;
      }

      case 3:
      {
                                         // same as in 2d
        Assert ((degree+1)*(degree+1) == this->dofs_per_face, ExcInternalError());
        
                                         // start with the face shape
                                         // functions
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[0*this->dofs_per_face+i] = Point<dim>(.5, .0, .5);
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[1*this->dofs_per_face+i] = Point<dim>(.5, 1., .5);
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[2*this->dofs_per_face+i] = Point<dim>(.5, .5, 0.);
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[3*this->dofs_per_face+i] = Point<dim>(1., .5, .5);
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[4*this->dofs_per_face+i] = Point<dim>(.5, .5, 1.);
        for (unsigned int i=0; i<this->dofs_per_face; ++i)
          this->unit_support_points[5*this->dofs_per_face+i] = Point<dim>(.0, .5, .5);

                                         // associate the rest with
                                         // the cell center
        for (unsigned int i=6*this->dofs_per_face; i<this->dofs_per_cell; ++i)
          this->unit_support_points[i] = Point<dim>(.5, .5, .5);
        
        break;
      }

      default:
	    Assert (false, ExcNotImplemented());
    };
}


#if deal_II_dimension == 1

template <>
void FE_RaviartThomas<1>::initialize_unit_face_support_points ()
{
				   // no faces in 1d, so nothing to do
}

#endif


template <int dim>
void FE_RaviartThomas<dim>::initialize_unit_face_support_points ()
{
  this->unit_face_support_points.resize (this->dofs_per_face);

                                         // like with cell
                                         // unit_support_points:
                                         // associate all of the in
                                         // the face mid-point, since
                                         // there is no other useful
                                         // way
  for (unsigned int i=0; i<this->dofs_per_face; ++i)
    this->unit_face_support_points[i] = (dim == 2 ?
                                         Point<dim-1>(.5) :
                                         Point<dim-1>(.5,.5));
}


#if deal_II_dimension == 1

template <>
std::vector<unsigned int>
FE_RaviartThomas<1>::get_dpo_vector (const unsigned int)
{
  Assert (false, ExcNotUsefulInThisDimension());
  return std::vector<unsigned int>();
}

#endif


template <int dim>
std::vector<unsigned int>
FE_RaviartThomas<dim>::get_dpo_vector (const unsigned int degree)
{
                                   // the element is face-based (not
                                   // to be confused with George
                                   // W. Bush's Faith Based
                                   // Initiative...), and we have
                                   // (degree+1)^(dim-1) DoFs per face
  unsigned int dofs_per_face = 1;
  for (unsigned int d=0; d<dim-1; ++d)
    dofs_per_face *= degree+1;

                                   // and then there are interior dofs
  const unsigned int
    interior_dofs = dim*degree*dofs_per_face;
  
  std::vector<unsigned int> dpo(dim+1);
  dpo[dim-1] = dofs_per_face;
  dpo[dim]   = interior_dofs;
  
  return dpo;
}



#if deal_II_dimension == 1

template <>
std::vector<bool>
FE_RaviartThomas<1>::get_ria_vector (const unsigned int)
{
  Assert (false, ExcNotUsefulInThisDimension());
  return std::vector<bool>();
}

#endif


template <int dim>
std::vector<bool>
FE_RaviartThomas<dim>::get_ria_vector (const unsigned int degree)
{
  unsigned int dofs_per_cell, dofs_per_face;
  switch (dim)
    {
      case 2:
	    dofs_per_face = degree+1;
	    dofs_per_cell = 2*(degree+1)*(degree+2);
	    break;
      case 3:
	    dofs_per_face = (degree+1)*(degree+1);
	    dofs_per_cell = 3*(degree+1)*(degree+1)*(degree+2);
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    }
  Assert (FiniteElementData<dim>(get_dpo_vector(degree),dim).dofs_per_cell ==
	  dofs_per_cell,
	  ExcInternalError());
  Assert (FiniteElementData<dim>(get_dpo_vector(degree),dim).dofs_per_face ==
	  dofs_per_face,
	  ExcInternalError());
  
				   // all face dofs need to be
				   // non-additive, since they have
				   // continuity requirements.
				   // however, the interior dofs are
				   // made additive
  std::vector<bool> ret_val(dofs_per_cell,false);
  for (unsigned int i=GeometryInfo<dim>::faces_per_cell*dofs_per_face;
       i < dofs_per_cell; ++i)
    ret_val[i] = true;

  return ret_val;
}


#if deal_II_dimension == 1

template <>
std::vector<AnisotropicPolynomials<1> >
FE_RaviartThomas<1>::create_polynomials (const unsigned int)
{
  Assert (false, ExcNotUsefulInThisDimension());
  return std::vector<AnisotropicPolynomials<1> > ();
}

#endif


#if deal_II_dimension == 2

template <>
std::vector<AnisotropicPolynomials<2> >
FE_RaviartThomas<2>::create_polynomials (const unsigned int degree)
{
  const unsigned int dim = 2;
  
                                   // use the fact that the RT(k)
                                   // spaces are spanned by the
                                   // functions
                                   // P_{k+1,k} \times P_{k,k+1},
                                   // see the book by Brezzi and
                                   // Fortin
  const std::vector<Polynomials::Polynomial<double> > pols[2]
    = { Polynomials::Hierarchical::generate_complete_basis (degree+1),
        Polynomials::Hierarchical::generate_complete_basis (degree)};

                                   // create spaces (k+1,k) and (k,k+1)
  std::vector<std::vector<Polynomials::Polynomial<double> > >
    pols_vector_1(dim), pols_vector_2(dim);
  pols_vector_1[0] = pols[0];
  pols_vector_1[1] = pols[1];

  pols_vector_2[0] = pols[1];
  pols_vector_2[1] = pols[0];
  
  const AnisotropicPolynomials<dim> anisotropic[dim]
    = { AnisotropicPolynomials<dim> (pols_vector_1),
        AnisotropicPolynomials<dim> (pols_vector_2) };

                                   // work around a stupid bug in
                                   // gcc2.95 where the compiler
                                   // complains about reaching the end
                                   // of a non-void function when we
                                   // simply return the following
                                   // object unnamed, rather than
                                   // first creating a named object
                                   // and then returning it...
  const std::vector<AnisotropicPolynomials<dim> >
    ret_val (&anisotropic[0], &anisotropic[dim]);
  return ret_val;
}

#endif


#if deal_II_dimension == 3

template <>
std::vector<AnisotropicPolynomials<3> >
FE_RaviartThomas<3>::create_polynomials (const unsigned int degree)
{
  const unsigned int dim = 3;
  
                                   // use the fact that the RT(k)
                                   // spaces are spanned by the
                                   // functions
                                   // P_{k+1,k,k} \times P_{k,k+1,k}
                                   // \times P_{k,k,k+1},
                                   // see the book by Brezzi and
                                   // Fortin
  const std::vector<Polynomials::Polynomial<double> > pols[2]
    = { Polynomials::Hierarchical::generate_complete_basis (degree+1),
        Polynomials::Hierarchical::generate_complete_basis (degree)};

                                   // create spaces (k+1,k,k),
                                   // (k,k+1,k) and (k,k,k+1)
  std::vector<std::vector<Polynomials::Polynomial<double> > >
    pols_vector_1(dim), pols_vector_2(dim), pols_vector_3(dim);
  pols_vector_1[0] = pols[0];
  pols_vector_1[1] = pols[1];
  pols_vector_1[2] = pols[1];

  pols_vector_2[0] = pols[1];
  pols_vector_2[1] = pols[0];
  pols_vector_2[2] = pols[1];
  
  pols_vector_3[0] = pols[1];
  pols_vector_3[1] = pols[1];
  pols_vector_3[2] = pols[0];
  
  const AnisotropicPolynomials<dim> anisotropic[dim]
    = { AnisotropicPolynomials<dim> (pols_vector_1),
        AnisotropicPolynomials<dim> (pols_vector_2),
        AnisotropicPolynomials<dim> (pols_vector_3) };

                                   // work around a stupid bug in
                                   // gcc2.95 where the compiler
                                   // complains about reaching the end
                                   // of a non-void function when we
                                   // simply return the following
                                   // object unnamed, rather than
                                   // first creating a named object
                                   // and then returning it...
  const std::vector<AnisotropicPolynomials<dim> >
    ret_val (&anisotropic[0], &anisotropic[dim]);
  return ret_val;
}

#endif



#if deal_II_dimension == 1

template <>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomas<1>::compute_renumber (const unsigned int)
{
  Assert (false, ExcNotUsefulInThisDimension());
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}

#endif


#if deal_II_dimension == 2

template <>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomas<2>::compute_renumber (const unsigned int degree)
{
  const unsigned int dim = 2;
  
  std::vector<std::pair<unsigned int, unsigned int> > ret_val;
  
                                   // to explain the following: the
                                   // first (degree+1) shape functions
                                   // are on face 0, and point in
                                   // y-direction, so are for the
                                   // second vector component. then
                                   // there are (degree+1) shape
                                   // functions on face 1, which is
                                   // for the x vector component, and
                                   // so on. since the order of face
                                   // degrees of freedom is arbitrary,
                                   // we simply use the same order as
                                   // that provided by the 1d
                                   // polynomial class on which this
                                   // element is based. after
                                   // 4*(degree+1), the remaining
                                   // shape functions are all bubbles,
                                   // so we can number them in any way
                                   // we want. we do so by first
                                   // numbering the x-vectors, then
                                   // the y-vectors
                                   //
                                   // now, we have to find a mapping
                                   // from the above ordering to:
                                   // first which vector component
                                   // they belong to (easy), and
                                   // second the index within this
                                   // component as provided by the
                                   // AnisotropicPolynomials class
                                   //
                                   // this is mostly a counting
                                   // argument, tedious and error
                                   // prone, and so boring to explain
                                   // that we rather not try to do so
                                   // here (it's simple, but boring,
                                   // as said), aside from a few
                                   // comments below

                                   // face 0
  for (unsigned int i=0; i<degree+1; ++i)
    ret_val.push_back (std::make_pair (1, i));
  
                                   // face 1
  for (unsigned int i=0; i<degree+1; ++i)
    ret_val.push_back (std::make_pair (0, (degree+2)*i+1));
  
                                   // face 2
  for (unsigned int i=0; i<degree+1; ++i)
    ret_val.push_back (std::make_pair (1, (degree+1)+i));
  
                                   // face 3
  for (unsigned int i=0; i<degree+1; ++i)
    ret_val.push_back (std::make_pair (0, (degree+2)*i));

                                   // then go on with interior bubble
                                   // functions, first for the
                                   // x-direction, then for the
                                   // y-direction
  for (unsigned int x=0; x<degree; ++x)
    for (unsigned int y=0; y<degree+1; ++y)
      {
        const unsigned int index_in_component = (x+2) + y*(degree+2);
        Assert (index_in_component < (degree+1)*(degree+2),
                ExcInternalError());
        ret_val.push_back (std::make_pair(0, index_in_component));
      }
  for (unsigned int x=0; x<degree+1; ++x)
    for (unsigned int y=0; y<degree; ++y)
      {
        const unsigned int index_in_component = 2*(degree+1) + y + x*degree;
        Assert (index_in_component < (degree+1)*(degree+2),
                ExcInternalError());
        ret_val.push_back (std::make_pair(1, index_in_component));
      }

#ifdef DEBUG  
                                   // make sure we have actually used
                                   // up all elements of the tensor
                                   // product polynomial
  Assert (ret_val.size() == 2*(degree+1)*(degree+2),
          ExcInternalError());
  std::vector<bool> test[dim] = { std::vector<bool>(ret_val.size()/dim, false),
                                  std::vector<bool>(ret_val.size()/dim, false) };
  for (unsigned int i=0; i<ret_val.size(); ++i)
    {
      Assert (ret_val[i].first < dim, ExcInternalError());
      Assert (ret_val[i].second < test[ret_val[i].first].size(),
	      ExcInternalError());
      Assert (test[ret_val[i].first][ret_val[i].second] == false,
              ExcInternalError());
      
      test[ret_val[i].first][ret_val[i].second] = true;
    }
  for (unsigned int d=0; d<dim; ++d)
    Assert (std::find (test[d].begin(), test[d].end(), false) == test[d].end(),
            ExcInternalError());
#endif
  
  return ret_val;
}

#endif


#if deal_II_dimension == 3

template <>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomas<3>::compute_renumber (const unsigned int /*degree*/)
{
  Assert (false, ExcNotImplemented());
  return std::vector<std::pair<unsigned int, unsigned int> > ();  
}

#endif




template <int dim>
UpdateFlags
FE_RaviartThomas<dim>::update_once (const UpdateFlags) const
{
				   // even the values have to be
				   // computed on the real cell, so
				   // nothing can be done in advance
  return update_default;
}



template <int dim>
UpdateFlags
FE_RaviartThomas<dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_values)
    out |= update_values             | update_covariant_transformation;
  if (flags & update_gradients)
    out |= update_gradients          | update_covariant_transformation;
  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_covariant_transformation;

  return out;
}



//----------------------------------------------------------------------
// Data field initialization
//----------------------------------------------------------------------

template <int dim>
typename Mapping<dim>::InternalDataBase *
FE_RaviartThomas<dim>::get_data (const UpdateFlags      update_flags,
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

 				   // initialize fields only if really
 				   // necessary. otherwise, don't
 				   // allocate memory
  if (flags & update_values)
    data->shape_values.resize (this->dofs_per_cell,
                               std::vector<Tensor<1,dim> > (n_q_points));

  if (flags & update_gradients)
    data->shape_gradients.resize (this->dofs_per_cell,
                                  std::vector<Tensor<2,dim> > (n_q_points));

 				   // if second derivatives through
 				   // finite differencing is required,
 				   // then initialize some objects for
 				   // that
  if (flags & update_second_derivatives)
    data->initialize_2nd (this, mapping, quadrature);

 				   // next already fill those fields
 				   // of which we have information by
 				   // now. note that the shape values
 				   // and gradients are only those on
 				   // the unit cell, and need to be
 				   // transformed when visiting an
 				   // actual cell
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    for (unsigned int q=0; q<n_q_points; ++q)
      {
        if (flags & update_values)
          for (unsigned int c=0; c<dim; ++c)
            data->shape_values[i][q][c]
              = shape_value_component(i,quadrature.point(q),c);
	
        if (flags & update_gradients)
          for (unsigned int c=0; c<dim; ++c)
            data->shape_gradients[i][q][c]
              = shape_grad_component(i,quadrature.point(q),c);
      }
   
  return data;
}




//----------------------------------------------------------------------
// Fill data of FEValues
//----------------------------------------------------------------------

template <int dim>
void
FE_RaviartThomas<dim>::fill_fe_values (const Mapping<dim>                   &mapping,
                                       const typename Triangulation<dim>::cell_iterator &cell,
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

				   // get the flags indicating the
				   // fields that have to be filled
  const UpdateFlags flags(fe_data.current_update_flags());

  const unsigned int n_q_points = quadrature.n_quadrature_points;
				  
				   // fill shape function
				   // values. these are vector-valued,
				   // so we have to transform
				   // them. since the output format
				   // (in data.shape_values) is a
				   // sequence of doubles (one for
				   // each non-zero shape function
				   // value, and for each quadrature
				   // point, rather than a sequence of
				   // small vectors, we have to use a
				   // number of conversions
  if (flags & update_values)
    {
      std::vector<Tensor<1,dim> > shape_values (n_q_points);

      Assert (data.shape_values.n_rows() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_values.n_cols() == n_q_points,
	      ExcInternalError());
      
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
					   // first transform shape
					   // values...
	  Assert (fe_data.shape_values[k].size() == n_q_points,
		  ExcInternalError());
	  mapping.transform_covariant(fe_data.shape_values[k], 0,
                                      shape_values,
                                      mapping_data);

					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_values[k*dim+d][q] = shape_values[q][d];
	};
    };
  
      
  if (flags & update_gradients)
    {
      std::vector<Tensor<2,dim> > shape_grads1 (n_q_points);
      std::vector<Tensor<2,dim> > shape_grads2 (n_q_points);

      Assert (data.shape_gradients.size() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_gradients[0].size() == n_q_points,
	      ExcInternalError());

                                       // loop over all shape
                                       // functions, and treat the
                                       // gradients of each shape
                                       // function at all quadrature
                                       // points
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
                                           // treat the gradients of
                                           // this particular shape
                                           // function at all
                                           // q-points. if Dv is the
                                           // gradient of the shape
                                           // function on the unit
                                           // cell, then
                                           // (J^-T)Dv(J^-1) is the
                                           // value we want to have on
                                           // the real cell. so, we
                                           // will have to apply a
                                           // covariant transformation
                                           // to Dv twice. since the
                                           // interface only allows
                                           // multiplication with
                                           // (J^-1) from the right,
                                           // we have to trick a
                                           // little in between
	  Assert (fe_data.shape_gradients[k].size() == n_q_points,
		  ExcInternalError());
                                           // do first transformation
	  mapping.transform_covariant(fe_data.shape_gradients[k], 0,
                                      shape_grads1,
                                      mapping_data);
                                           // transpose matrix
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
                                           // do second transformation
	  mapping.transform_covariant(shape_grads2, 0, shape_grads1,
                                      mapping_data);
                                           // transpose back
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
          
					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_gradients[k*dim+d][q] = shape_grads2[q][d];
	};
    }

  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell,
                       QProjector<dim>::DataSetDescriptor::cell(),
                       mapping_data, fe_data, data);
}



template <int dim>
void
FE_RaviartThomas<dim>::fill_fe_face_values (const Mapping<dim>                   &mapping,
                                            const typename Triangulation<dim>::cell_iterator &cell,
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

  				   // get the flags indicating the
				   // fields that have to be filled
  const UpdateFlags flags(fe_data.current_update_flags());

  const unsigned int n_q_points = quadrature.n_quadrature_points;
				  
				   // fill shape function
				   // values. these are vector-valued,
				   // so we have to transform
				   // them. since the output format
				   // (in data.shape_values) is a
				   // sequence of doubles (one for
				   // each non-zero shape function
				   // value, and for each quadrature
				   // point, rather than a sequence of
				   // small vectors, we have to use a
				   // number of conversions
  if (flags & update_values)
    {
      Assert (fe_data.shape_values.size() == this->dofs_per_cell,
              ExcInternalError());
      Assert (fe_data.shape_values[0].size() ==
              GeometryInfo<dim>::faces_per_cell * n_q_points *
	      (dim == 3 ? 2 : 1),
              ExcInternalError());
      
      std::vector<Tensor<1,dim> > shape_values (n_q_points);

      Assert (data.shape_values.n_rows() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_values.n_cols() == n_q_points,
	      ExcInternalError());
      
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
					   // first transform shape
					   // values...
	  mapping.transform_covariant(fe_data.shape_values[k], offset,
                                      shape_values,
                                      mapping_data);

					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_values[k*dim+d][q] = shape_values[q][d];
	};
    };
  
      
  if (flags & update_gradients)
    {
      Assert (fe_data.shape_values.size() == this->dofs_per_cell,
              ExcInternalError());
      Assert (fe_data.shape_gradients[0].size() ==
              GeometryInfo<dim>::faces_per_cell * n_q_points *
	      (dim == 3 ? 2 : 1),
              ExcInternalError());

      std::vector<Tensor<2,dim> > shape_grads1 (n_q_points);
      std::vector<Tensor<2,dim> > shape_grads2 (n_q_points);

      Assert (data.shape_gradients.size() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_gradients[0].size() == n_q_points,
	      ExcInternalError());

                                       // loop over all shape
                                       // functions, and treat the
                                       // gradients of each shape
                                       // function at all quadrature
                                       // points
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
                                           // treat the gradients of
                                           // this particular shape
                                           // function at all
                                           // q-points. if Dv is the
                                           // gradient of the shape
                                           // function on the unit
                                           // cell, then
                                           // (J^-T)Dv(J^-1) is the
                                           // value we want to have on
                                           // the real cell. so, we
                                           // will have to apply a
                                           // covariant transformation
                                           // to Dv twice. since the
                                           // interface only allows
                                           // multiplication with
                                           // (J^-1) from the right,
                                           // we have to trick a
                                           // little in between
                                           // 
                                           // do first transformation
	  mapping.transform_covariant(fe_data.shape_gradients[k], offset,
                                      shape_grads1,
                                      mapping_data);
                                           // transpose matrix
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
                                           // do second transformation
	  mapping.transform_covariant(shape_grads2, 0, shape_grads1,
                                      mapping_data);
                                           // transpose back
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
          
					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_gradients[k*dim+d][q] = shape_grads2[q][d];
	};
    }

  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <int dim>
void
FE_RaviartThomas<dim>::fill_fe_subface_values (const Mapping<dim>                   &mapping,
                                               const typename Triangulation<dim>::cell_iterator &cell,
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
				   // faces are stored contiguously)
  const typename QProjector<dim>::DataSetDescriptor offset
    = (QProjector<dim>::DataSetDescriptor::
       sub_face (face, subface, cell->face_orientation(face),
                 quadrature.n_quadrature_points));

  				   // get the flags indicating the
				   // fields that have to be filled
  const UpdateFlags flags(fe_data.current_update_flags());

  const unsigned int n_q_points = quadrature.n_quadrature_points;
				  
				   // fill shape function
				   // values. these are vector-valued,
				   // so we have to transform
				   // them. since the output format
				   // (in data.shape_values) is a
				   // sequence of doubles (one for
				   // each non-zero shape function
				   // value, and for each quadrature
				   // point, rather than a sequence of
				   // small vectors, we have to use a
				   // number of conversions
  if (flags & update_values)
    {
      Assert (fe_data.shape_values[0].size() ==
              GeometryInfo<dim>::faces_per_cell *
	      GeometryInfo<dim>::subfaces_per_face *
	      n_q_points,
              ExcInternalError());
      
      std::vector<Tensor<1,dim> > shape_values (n_q_points);

      Assert (data.shape_values.n_rows() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_values.n_cols() == n_q_points,
	      ExcInternalError());
      
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
					   // first transform shape
					   // values...
	  mapping.transform_covariant(fe_data.shape_values[k], offset,
                                      shape_values,
                                      mapping_data);

					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_values[k*dim+d][q] = shape_values[q][d];
	};
    };
  
      
  if (flags & update_gradients)
    {
      Assert (fe_data.shape_gradients.size() ==
              GeometryInfo<dim>::faces_per_cell *
	      GeometryInfo<dim>::subfaces_per_face *
	      n_q_points,
              ExcInternalError());

      std::vector<Tensor<2,dim> > shape_grads1 (n_q_points);
      std::vector<Tensor<2,dim> > shape_grads2 (n_q_points);

      Assert (data.shape_gradients.size() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_gradients[0].size() == n_q_points,
	      ExcInternalError());

                                       // loop over all shape
                                       // functions, and treat the
                                       // gradients of each shape
                                       // function at all quadrature
                                       // points
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
                                           // treat the gradients of
                                           // this particular shape
                                           // function at all
                                           // q-points. if Dv is the
                                           // gradient of the shape
                                           // function on the unit
                                           // cell, then
                                           // (J^-T)Dv(J^-1) is the
                                           // value we want to have on
                                           // the real cell. so, we
                                           // will have to apply a
                                           // covariant transformation
                                           // to Dv twice. since the
                                           // interface only allows
                                           // multiplication with
                                           // (J^-1) from the right,
                                           // we have to trick a
                                           // little in between
                                           // 
                                           // do first transformation
	  mapping.transform_covariant(fe_data.shape_gradients[k], offset,
                                      shape_grads1,
                                      mapping_data);
                                           // transpose matrix
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
                                           // do second transformation
	  mapping.transform_covariant(shape_grads2, 0, shape_grads1,
                                      mapping_data);
                                           // transpose back
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
          
					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_gradients[k*dim+d][q] = shape_grads2[q][d];
	};
    }

  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <int dim>
unsigned int
FE_RaviartThomas<dim>::n_base_elements () const
{
  return 1;
}



template <int dim>
const FiniteElement<dim> &
FE_RaviartThomas<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <int dim>
unsigned int
FE_RaviartThomas<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}



template <int dim>
bool
FE_RaviartThomas<dim>::has_support_on_face (const unsigned int shape_index,
                                            const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

				   // Return computed values if we
				   // know them easily. Otherwise, it
				   // is always safe to return true.
  switch (degree)
    {
      case 0:
      {
        switch (dim)
          {
            case 2:
            {
                                               // only on the one
                                               // non-adjacent face
                                               // are the values
                                               // actually zero. list
                                               // these in a table
              const unsigned int
                opposite_faces[GeometryInfo<2>::faces_per_cell]
                = { 2, 3, 0, 1};
              
              return (face_index != opposite_faces[shape_index]);
            }
            
            default:
	      return true;
          };
      };
      
      default:  // other degree
	return true;
    };
  
  return true;
}



template <int dim>
unsigned int
FE_RaviartThomas<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template <int dim>
unsigned int
FE_RaviartThomas<dim>::get_degree () const
{
  return degree;
}



template class FE_RaviartThomas<deal_II_dimension>;
