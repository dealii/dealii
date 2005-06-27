//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/quadrature.h>
#include <base/qprojector.h>
#include <fe/fe_q.h>
#include <fe/fe_tools.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


// namespace for some functions that are used in this file. they are
// specific to numbering conventions used for the FE_Q element, and
// are thus not very interesting to the outside world. we'd like to
// simply put them into an anonymous namespace, but that triggers an
// odd error with icc which can't compile this small snippet if the
// function is static:
// --------------------
// template <int> struct int2type {};
// 
// namespace {    
//   static void SYMBOL (const int2type<1> & ) {}
// }
// 
// template <int dim> void g() {
//   SYMBOL(int2type<dim>());
// }
//
// template void g<1>();
// --------------------
// the function needs to be static because of another icc bug, though.
// work around this by packing everything into a namespace of its own
// and have the anonymous namespace inside
//
// this is now intel icc issue 216082
namespace FE_Q_Helper
{
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
    inline
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
#if defined(DEAL_II_ANON_NAMESPACE_BUG) && defined(DEAL_II_ANON_NAMESPACE_LINKAGE_BUG)
    static
#endif
    inline
    unsigned int int_sqrt (const unsigned int N)
    {
      for (unsigned int i=0; i<=N; ++i)
	if (i*i == N)
	  return i;
      Assert (false, ExcInternalError());
      return deal_II_numbers::invalid_unsigned_int;
    }


				     // given an integer N, compute its
				     // integer cube root (if it
				     // exists, otherwise give up)
#if defined(DEAL_II_ANON_NAMESPACE_BUG) && defined(DEAL_II_ANON_NAMESPACE_LINKAGE_BUG)
    static
#endif
    inline
    unsigned int int_cuberoot (const unsigned int N)
    {
      for (unsigned int i=0; i<=N; ++i)
	if (i*i*i == N)
	  return i;
      Assert (false, ExcInternalError());
      return deal_II_numbers::invalid_unsigned_int;
    }


				     // given N, generate i=0...N-1
				     // equidistant points in the
				     // interior of the interval [0,1]
#if defined(DEAL_II_ANON_NAMESPACE_BUG) && defined(DEAL_II_ANON_NAMESPACE_LINKAGE_BUG)
    static
#endif
    inline
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
#if defined(DEAL_II_ANON_NAMESPACE_BUG) && defined(DEAL_II_ANON_NAMESPACE_LINKAGE_BUG)
    static
#endif
    inline
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
#if defined(DEAL_II_ANON_NAMESPACE_BUG) && defined(DEAL_II_ANON_NAMESPACE_LINKAGE_BUG)
    static
#endif
    inline
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
    
  }
}



template <int dim>
FE_Q<dim>::FE_Q (const unsigned int degree)
		:
		FE_Poly<TensorProductPolynomials<dim>, dim> (
		  TensorProductPolynomials<dim>(Polynomials::LagrangeEquidistant::generate_complete_basis(degree)),
		  FiniteElementData<dim>(get_dpo_vector(degree),1, degree),
		  std::vector<bool> (FiniteElementData<dim>(
		    get_dpo_vector(degree),1, degree).dofs_per_cell, false),
		  std::vector<std::vector<bool> >(FiniteElementData<dim>(
		    get_dpo_vector(degree),1, degree).dofs_per_cell, std::vector<bool>(1,true))),
						    face_index_map(FE_Q_Helper::invert_numbering(face_lexicographic_to_hierarchic_numbering (degree)))
{
  std::vector<unsigned int> renumber (this->dofs_per_cell);
  FETools::hierarchic_to_lexicographic_numbering (*this, renumber);
  this->poly_space.set_numbering(renumber);
  
				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();

				   // compute constraint, embedding
				   // and restriction matrices
  initialize_constraints ();
  initialize_embedding ();
  initialize_restriction ();
}



template <int dim>
std::string
FE_Q<dim>::get_name () const
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
  
  namebuf << "FE_Q<" << dim << ">(" << this->degree << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_Q<dim>::clone() const
{
  return new FE_Q<dim>(this->degree);
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

  const std::vector<unsigned int> &index_map=
    this->poly_space.get_numbering();
  
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
      const Point<dim>
	p = FE_Q_Helper::generate_unit_point (index_map[j], this->dofs_per_cell,
					      FE_Q_Helper::int2type<dim>());
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        cell_interpolation(j,i) = this->poly_space.compute_value (i, p);

      for (unsigned int i=0; i<source_fe.dofs_per_cell; ++i)
        source_interpolation(j,i) = source_fe.poly_space.compute_value (i, p);
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

      Assert (std::fabs(sum-1) < 2e-14*this->degree*dim,
              ExcInternalError());
    }
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------



template <int dim>
void FE_Q<dim>::initialize_unit_support_points ()
{
				   // number of points: (degree+1)^dim
  unsigned int n = this->degree+1;
  for (unsigned int i=1; i<dim; ++i)
    n *= this->degree+1;
  
  this->unit_support_points.resize(n);

  const std::vector<unsigned int> &index_map_inverse=
    this->poly_space.get_numbering_inverse();
  
  const double step = 1./this->degree;
  Point<dim> p;
  
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((dim>2) ? this->degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((dim>1) ? this->degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=this->degree; ++ix)
	{
	  p(0) = ix * step;
	  if (dim>1)
	    p(1) = iy * step;
	  if (dim>2)
	    p(2) = iz * step;
	  
	  this->unit_support_points[index_map_inverse[k++]] = p;
	}
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
  unsigned int n = this->degree+1;
  for (unsigned int i=1; i<codim; ++i)
    n *= this->degree+1;
  
  this->unit_face_support_points.resize(n);

  const std::vector<unsigned int> &face_index_map_inverse=
    FE_Q_Helper::invert_numbering(face_index_map);
  
  const double step = 1./this->degree;
  Point<codim> p;
  
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((codim>2) ? this->degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((codim>1) ? this->degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=this->degree; ++ix)
	{
	  p(0) = ix * step;
	  if (codim>1)
	    p(1) = iy * step;
	  if (codim>2)
	    p(2) = iz * step;
	  
	  this->unit_face_support_points[face_index_map_inverse[k++]] = p;
	}
}



template <int dim>
std::vector<unsigned int>
FE_Q<dim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}


template <int dim>
std::vector<unsigned int>
FE_Q<dim>::face_lexicographic_to_hierarchic_numbering (const unsigned int degree)
{
  const FiniteElementData<dim-1> face_data(FE_Q<dim-1>::get_dpo_vector(degree),1);
  std::vector<unsigned int> face_renumber (face_data.dofs_per_cell);  
  FETools::lexicographic_to_hierarchic_numbering (face_data, face_renumber);
  return face_renumber;
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
				   // degree of the element.  from
				   // this, we interpolate between
				   // mother and cell face.

				   // the interpolation process works
				   // as followings: on each subface,
				   // we want that finite element
				   // solutions from both sides
				   // coincide. i.e. if a and b are
				   // expansion coefficients for the
				   // shape functions from both sides,
				   // we seek a relation between x and
				   // y such that
				   //   sum_j a_j phi^c_j(x)
				   //   == sum_j b_j phi_j(x)  
				   // for all points x on the
				   // interface. here, phi^c_j are the
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
				   // we obtain the matrix system
				   //    A a  ==  B b
				   // where
				   //    A_ij = phi^c_j(x_i)
  				   //    B_ij = phi_j(x_i)
				   // and the relation we are looking
				   // for is
				   //    a = A^-1 B b
				   //
				   // for the special case of Lagrange
				   // interpolation polynomials, A_ij
				   // reduces to delta_ij, and
				   //    a_i = B_ij b_j
				   // Hence,
				   // interface_constraints(i,j)=B_ij.
				   //
				   // for the general case, where we
				   // don't have Lagrange
				   // interpolation polynomials, this
				   // is a little more
				   // complicated. Then we would
				   // evaluate at a number of points
				   // and invert the interpolation
				   // matrix A.
				   //
				   // Note, that we build up these
				   // matrices for all subfaces at
				   // once, rather than considering
				   // them separately. the reason is
				   // that we finally will want to
				   // have them in this order anyway,
				   // as this is the format we need
				   // inside deal.II

				   // In the following the points x_i
				   // are constructed in following
				   // order (n=degree-1)
				   // *----------*---------*
				   //     1..n   0  n+1..2n
				   // i.e. first the midpoint of the
				   // line, then the support points on
				   // subface 0 and on subface 1
  std::vector<Point<dim-1> > constraint_points;
                                   // Add midpoint
  constraint_points.push_back (Point<dim-1> (0.5));

  if (this->degree>1)
    {
      const unsigned int n=this->degree-1;
      const double step=1./this->degree;
				       // subface 0
      for (unsigned int i=1; i<=n; ++i)
	constraint_points.push_back (
	  GeometryInfo<dim-1>::child_to_cell_coordinates(Point<dim-1>(i*step),0));
				       // subface 1
      for (unsigned int i=1; i<=n; ++i)
	constraint_points.push_back (
	  GeometryInfo<dim-1>::child_to_cell_coordinates(Point<dim-1>(i*step),1));
    }

                                   // Now construct relation between
                                   // destination (child) and source (mother)
                                   // dofs.
  const std::vector<Polynomials::Polynomial<double> > polynomials=
    Polynomials::LagrangeEquidistant::generate_complete_basis(this->degree);

  this->interface_constraints
    .TableBase<2,double>::reinit (this->interface_constraints_size());

  for (unsigned int i=0; i<constraint_points.size(); ++i)
    for (unsigned j=0; j<this->degree+1; ++j)
      {
	this->interface_constraints(i,j) = 
	  polynomials[face_index_map[j]].value (constraint_points[i](0));
		    
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
          if (std::fabs(this->interface_constraints(i,j)) < 1e-14)
            this->interface_constraints(i,j) = 0;
      }
}

#endif

#if deal_II_dimension == 3
template <>
void
FE_Q<3>::initialize_constraints ()
{
  const unsigned int dim = 3;

				   // For a detailed documentation of
				   // the interpolation see the
				   // FE_Q<2>::initialize_constraints
				   // function.

				   // In the following the points x_i
				   // are constructed in the order as
				   // described in the documentation
				   // of the FiniteElementBase class
				   // (fe_base.h), i.e.
				   //   *--13--3--14--*
				   //   |      |      |
				   //   16 20  7  19  12
				   //   |      |      |
				   //   4--8---0--6---2
				   //   |      |      |
				   //   15 17  5  18  11
				   //   |      |      |
				   //   *--9---1--10--*
  std::vector<Point<dim-1> > constraint_points;

                                   // Add midpoint
  constraint_points.push_back (Point<dim-1> (0.5, 0.5));

                                   // Add midpoints of lines of
                                   // "mother-face"
  constraint_points.push_back (Point<dim-1> (0.5, 0));
  constraint_points.push_back (Point<dim-1> (1, 0.5));
  constraint_points.push_back (Point<dim-1> (0.5, 1));
  constraint_points.push_back (Point<dim-1> (0, 0.5));

  if (this->degree>1)
    {
      const unsigned int n=this->degree-1;
      const double step=1./this->degree;
      vector<Point<dim-2> > line_support_points(n);
      for (unsigned int i=0; i<n; ++i)
	line_support_points[i](0)=(i+1)*step;
      Quadrature<dim-2> qline(line_support_points);

				       // auxiliary points in 2d
      vector<Point<dim-1> > p_line(n);
  
				       // Add nodes of lines interior
				       // in the "mother-face"

				       // line 5: use line 15
      QProjector<dim-1>::project_to_subface(qline, 3, 0, p_line);
      for (unsigned int i=0; i<n; ++i)
	constraint_points.push_back (p_line[i] + Point<dim-1> (0.5, 0));
				       // line 6: use line 10
      QProjector<dim-1>::project_to_subface(qline, 0, 1, p_line);
      for (unsigned int i=0; i<n; ++i)
	constraint_points.push_back (p_line[i] + Point<dim-1> (0, 0.5));
				       // line 7: use line 16
      QProjector<dim-1>::project_to_subface(qline, 3, 1, p_line);
      for (unsigned int i=0; i<n; ++i)
	constraint_points.push_back (p_line[i] + Point<dim-1> (0.5, 0));
				       // line 8: use line 9
      QProjector<dim-1>::project_to_subface(qline, 0, 0, p_line);
      for (unsigned int i=0; i<n; ++i)
	constraint_points.push_back (p_line[i] + Point<dim-1> (0, 0.5));      
      
				       // DoFs on bordering lines
				       // lines 9-16
      for (unsigned int face=0; face<GeometryInfo<dim-1>::faces_per_cell; ++face)
	for (unsigned int subface=0;
	     subface<GeometryInfo<dim-1>::subfaces_per_face; ++subface)
	  {
	    QProjector<dim-1>::project_to_subface(qline, face, subface, p_line);
	    constraint_points.insert(constraint_points.end(),
				     p_line.begin(), p_line.end());
	  }
      
				       // Create constraints for
				       // interior nodes
      vector<Point<dim-1> > inner_points(n*n);
      for (unsigned int i=0, iy=1; iy<=n; ++iy)
	for (unsigned int ix=1; ix<=n; ++ix)
	  inner_points[i++] = Point<dim-1> (ix*step, iy*step);
      
      for (unsigned int child=0; 
	   child<GeometryInfo<dim-1>::children_per_cell; ++child)
	for (unsigned int i=0; i<inner_points.size(); ++i)
	  constraint_points.push_back (
	    GeometryInfo<dim-1>::child_to_cell_coordinates(inner_points[i], child));
    }

                                   // Now construct relation between
                                   // destination (child) and source (mother)
                                   // dofs.
  const unsigned int pnts=(this->degree+1)*(this->degree+1);  
  const std::vector<Polynomials::Polynomial<double> > polynomial_basis=
    Polynomials::LagrangeEquidistant::generate_complete_basis(this->degree);
 
  const TensorProductPolynomials<dim-1> face_polynomials(polynomial_basis);

  this->interface_constraints
    .TableBase<2,double>::reinit (this->interface_constraints_size());

  for (unsigned int i=0; i<constraint_points.size(); ++i)
    {
      const double interval = (double) (this->degree * 2);
      bool mirror[dim - 1];
      Point<dim-1> constraint_point;

                                       // Eliminate FP errors in constraint
                                       // points. Due to their origin, they
                                       // must all be fractions of the unit
                                       // interval. If we have polynomial
                                       // degree 4, the refined element has 8
                                       // intervals.  Hence the coordinates
                                       // must be 0, 0.125, 0.25, 0.375 etc.
                                       // Now the coordinates of the
                                       // constraint points will be multiplied
                                       // by the inverse of the interval size
                                       // (in the example by 8).  After that
                                       // the coordinates must be integral
                                       // numbers. Hence a normal truncation
                                       // is performed and the coordinates
                                       // will be scaled back. The equal
                                       // treatment of all coordinates should
                                       // eliminate any FP errors.
      for (unsigned int k=0; k<dim-1; ++k)
        {
          const int coord_int =
            static_cast<int> (constraint_points[i](k) * interval + 0.25);
          constraint_point(k) = 1.*coord_int / interval;

                                           // The following lines of code
                                           // should eliminate the problems
                                           // with the Constraint-Matrix,
                                           // which appeared for P>=4. The
                                           // Constraint-Matrix class
                                           // complained about different
                                           // constraints for the same entry
                                           // of the Constraint-Matrix.
                                           // Actually this difference could
                                           // be attributed to FP errors, as
                                           // it was in the range of
                                           // 1.0e-16. These errors originate
                                           // in the loss of symmetry in the
                                           // FP approximation of the
                                           // shape-functions.  Considering a
                                           // 3rd order shape function in 1D,
                                           // we have N0(x)=N3(1-x) and
                                           // N1(x)=N2(1-x).  For higher order
                                           // polynomials the FP
                                           // approximations of the shape
                                           // functions do not satisfy these
                                           // equations any more!  Thus in the
                                           // following code everything is
                                           // computed in the interval x \in
                                           // [0..0.5], which is sufficient to
                                           // express all values that could
                                           // come out from a computation of
                                           // any shape function in the full
                                           // interval [0..1]. If x > 0.5 the
                                           // computation is done for 1-x with
                                           // the shape function N_{p-n}
                                           // instead of N_n.  Hence symmetry
                                           // is preserved and everything
                                           // works fine...
                                           //
                                           // For a different explanation of
                                           // the problem, see the discussion
                                           // in the FiniteElementBase class
                                           // for constraint matrices in 3d.
          mirror[k] = (constraint_point(k) > 0.5);
          if (mirror[k])
            constraint_point(k) = 1.0 - constraint_point(k);
        }

      for (unsigned j=0; j<pnts; ++j)
        {
          unsigned int indices[2]
            = { face_index_map[j] % (this->degree + 1),
                face_index_map[j] / (this->degree + 1) };
          
          for (unsigned int k = 0; k<2; ++k)
            if (mirror[k])
              indices[k] = this->degree - indices[k];
          
          const unsigned int
            new_index = indices[1] * (this->degree + 1) + indices[0];

          this->interface_constraints(i,j) = 
            face_polynomials.compute_value (new_index, constraint_point);
	    
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
          if (std::fabs(this->interface_constraints(i,j)) < 1e-14)
            this->interface_constraints(i,j) = 0;
        }
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
  const std::vector<unsigned int> &index_map=
    this->poly_space.get_numbering();
  
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
	  const Point<dim> p_subcell
	    = FE_Q_Helper::generate_unit_point (index_map[j], this->dofs_per_cell,
						FE_Q_Helper::int2type<dim>());
	  const Point<dim> p_cell =
	    GeometryInfo<dim>::child_to_cell_coordinates (p_subcell, child);

	  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
	    {
	      const double
		cell_value    = this->poly_space.compute_value (i, p_cell),
		subcell_value = this->poly_space.compute_value (i, p_subcell);

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
	      if (std::fabs(cell_value) < 2e-14*this->degree*dim)
		cell_interpolation(j, i) = 0.;
	      else
		cell_interpolation(j, i) = cell_value;

	      if (std::fabs(subcell_value) < 2e-14*this->degree*dim)
		subcell_interpolation(j, i) = 0.;
	      else
		if (std::fabs(subcell_value-1) < 2e-14*this->degree*dim)
		  subcell_interpolation(j, i) = 1.;
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
	  if (std::fabs(this->prolongation[child](i,j)) < 2e-14*this->degree*dim)
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
	  Assert (std::fabs(sum-1.) < 2e-14*this->degree*dim,
		  ExcInternalError());
	}
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
      const Point<dim> p_cell
	= FE_Q_Helper::generate_unit_point (i, this->dofs_per_cell,
					    FE_Q_Helper::int2type<dim>());
      unsigned int mother_dof = 0;
      for (; mother_dof<this->dofs_per_cell; ++mother_dof)
        {
          const double val
            = this->poly_space.compute_value(mother_dof, p_cell);
          if (std::fabs (val-1.) < 2e-14*this->degree*dim)
                                             // ok, this is the right
                                             // dof
            break;
          else
                                             // make sure that all
                                             // other shape functions
                                             // are zero there
            Assert (std::fabs(val) < 2e-14*this->degree*dim,
                    ExcInternalError());
        }
                                       // check also the shape
                                       // functions after tat
      for (unsigned int j=mother_dof+1; j<this->dofs_per_cell; ++j)
        Assert (std::fabs (this->poly_space.compute_value(j, p_cell))
                < 2e-14*this->degree*dim,
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
                    = this->poly_space.compute_value(child_dof, p_subcell);
                  if (std::fabs (val-1.) < 2e-14*this->degree*dim)
                    break;
                  else
                    Assert (std::fabs(val) < 2e-14*this->degree*dim,
                            ExcInternalError());
                }
              for (unsigned int j=child_dof+1; j<this->dofs_per_cell; ++j)
                Assert (std::fabs (this->poly_space.compute_value(j, p_subcell))
                        < 2e-14*this->degree*dim,
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


//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------


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

	  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
	    if (GeometryInfo<dim>::face_to_cell_vertices(face_index, v) == vertex_no)
	      return true;

	  return false;
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
                                               // given face.
	      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
		if (GeometryInfo<3>::face_to_cell_lines(face_index, l) == line_index)
		  return true;

	      return false;
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
        }
    }

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



template class FE_Q<deal_II_dimension>;
