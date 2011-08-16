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

#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/base/quadrature_lib.h>

#include <vector>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


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
// template <int dim, int spacedim> void g() {
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



				// in initialize_embedding() and
				// initialize_restriction(), want to undo
				// tensorization on inner loops for
				// performance reasons. this clears a
				// dim-array
    template <int dim>
#ifdef DEAL_II_ANON_NAMESPACE_BUG
    static
#endif
    inline
    void
    zero_indices (unsigned int indices[dim])
    {
      for (unsigned int d=0; d<dim; ++d)
	indices[d] = 0;
    }



				// in initialize_embedding() and
				// initialize_restriction(), want to undo
				// tensorization on inner loops for
				// performance reasons. this increments tensor
				// product indices
    template <int dim>
#ifdef DEAL_II_ANON_NAMESPACE_BUG
    static
#endif
    inline
    void
    increment_indices (unsigned int       indices[dim],
		       const unsigned int dofs1d)
    {
      ++indices[0];
      for (int d=0; d<dim-1; ++d)
	if (indices[d]==dofs1d)
	  {
	    indices[d] = 0;
	    indices[d+1]++;
	  }
    }



				// in initialize_embedding() and
				// initialize_restriction(), want to undo
				// tensorization on inner loops for
				// performance reasons, and we need to again
				// access 1D polynomials. This function
				// creates them from dim-dimensional support
				// points.
    template <int dim>
#ifdef DEAL_II_ANON_NAMESPACE_BUG
    static
#endif
    inline
    std::vector<Polynomials::Polynomial<double> >
    generate_poly_space1d (const std::vector<Point<dim> >  &unit_support_points,
			   const std::vector<unsigned int> &index_map_inverse,
			   const unsigned int               dofs1d)
    {
      AssertDimension (Utilities::fixed_power<dim> (dofs1d),
		       unit_support_points.size());
      std::vector<Point<1> > points1d (dofs1d);
      for (unsigned int i=0; i<dofs1d; ++i)
	{
	  const unsigned int j = index_map_inverse[i];
	  points1d[i] = Point<1>(unit_support_points[j](0));
	  for (unsigned int d=1; d<dim; ++d)
	    Assert (unit_support_points[j][d] == 0.,
		    ExcInternalError());
	}
      return Polynomials::generate_complete_Lagrange_basis (points1d);
    }
  }
}



/**
 * A class with the same purpose as the similarly named class of the
 * Triangulation class. See there for more information.
 */
template <int xdim, int xspacedim>
struct FE_Q<xdim,xspacedim>::Implementation
{
				     /**
				      * Initialize the hanging node
				      * constraints matrices. Called from the
				      * constructor in case the finite element
				      * is based on quadrature points.
				      */
    template <int spacedim>
    static
    void initialize_constraints (const Quadrature<1> &,
				 FE_Q<1,spacedim>    &)
      {
					 // no constraints in 1d
      }


    template <int spacedim>
    static
    void initialize_constraints (const Quadrature<1> &points,
				 FE_Q<2,spacedim>    &fe)
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
					 // as follows: on each subface,
					 // we want that finite element
					 // solutions from both sides
					 // coincide. i.e. if a and b are
					 // expansion coefficients for the
					 // shape functions from both sides,
					 // we seek a relation between a and
					 // b such that
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

	if (fe.degree>1)
	  {
	    const unsigned int n=fe.degree-1;
	    const double step=1./fe.degree;
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
	  Polynomials::generate_complete_Lagrange_basis(points.get_points());

	fe.interface_constraints
	  .TableBase<2,double>::reinit (fe.interface_constraints_size());

	for (unsigned int i=0; i<constraint_points.size(); ++i)
	  for (unsigned j=0; j<fe.degree+1; ++j)
	    {
	      fe.interface_constraints(i,j) =
		polynomials[fe.face_index_map[j]].value (constraint_points[i](0));

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
	      if (std::fabs(fe.interface_constraints(i,j)) < 1e-14)
		fe.interface_constraints(i,j) = 0;
	    }
      }


    template <int spacedim>
    static
    void initialize_constraints (const Quadrature<1> &points,
				 FE_Q<3,spacedim>    &fe)
      {
	const unsigned int dim = 3;

					 // For a detailed documentation of
					 // the interpolation see the
					 // FE_Q<2>::initialize_constraints
					 // function.

					 // In the following the points x_i
					 // are constructed in the order as
					 // described in the documentation
					 // of the FiniteElement class
					 // (fe_base.h), i.e.
					 //   *--15--4--16--*
					 //   |      |      |
					 //   10 19  6  20  12
					 //   |      |      |
					 //   1--7---0--8---2
					 //   |      |      |
					 //   9  17  5  18  11
					 //   |      |      |
					 //   *--13--3--14--*
	std::vector<Point<dim-1> > constraint_points;

					 // Add midpoint
	constraint_points.push_back (Point<dim-1> (0.5, 0.5));

					 // Add midpoints of lines of
					 // "mother-face"
	constraint_points.push_back (Point<dim-1> (0, 0.5));
	constraint_points.push_back (Point<dim-1> (1, 0.5));
	constraint_points.push_back (Point<dim-1> (0.5, 0));
	constraint_points.push_back (Point<dim-1> (0.5, 1));

	if (fe.degree>1)
	  {
	    const unsigned int n=fe.degree-1;
	    const double step=1./fe.degree;
	    std::vector<Point<dim-2> > line_support_points(n);
	    for (unsigned int i=0; i<n; ++i)
	      line_support_points[i](0)=(i+1)*step;
	    Quadrature<dim-2> qline(line_support_points);

					     // auxiliary points in 2d
	    std::vector<Point<dim-1> > p_line(n);

					     // Add nodes of lines interior
					     // in the "mother-face"

					     // line 5: use line 9
	    QProjector<dim-1>::project_to_subface(qline, 0, 0, p_line);
	    for (unsigned int i=0; i<n; ++i)
	      constraint_points.push_back (p_line[i] + Point<dim-1> (0.5, 0));
					     // line 6: use line 10
	    QProjector<dim-1>::project_to_subface(qline, 0, 1, p_line);
	    for (unsigned int i=0; i<n; ++i)
	      constraint_points.push_back (p_line[i] + Point<dim-1> (0.5, 0));
					     // line 7: use line 13
	    QProjector<dim-1>::project_to_subface(qline, 2, 0, p_line);
	    for (unsigned int i=0; i<n; ++i)
	      constraint_points.push_back (p_line[i] + Point<dim-1> (0, 0.5));
					     // line 8: use line 14
	    QProjector<dim-1>::project_to_subface(qline, 2, 1, p_line);
	    for (unsigned int i=0; i<n; ++i)
	      constraint_points.push_back (p_line[i] + Point<dim-1> (0, 0.5));

					     // DoFs on bordering lines
					     // lines 9-16
	    for (unsigned int face=0; face<GeometryInfo<dim-1>::faces_per_cell; ++face)
	      for (unsigned int subface=0;
		   subface<GeometryInfo<dim-1>::max_children_per_face; ++subface)
		{
		  QProjector<dim-1>::project_to_subface(qline, face, subface, p_line);
		  constraint_points.insert(constraint_points.end(),
					   p_line.begin(), p_line.end());
		}

					     // Create constraints for
					     // interior nodes
	    std::vector<Point<dim-1> > inner_points(n*n);
	    for (unsigned int i=0, iy=1; iy<=n; ++iy)
	      for (unsigned int ix=1; ix<=n; ++ix)
		inner_points[i++] = Point<dim-1> (ix*step, iy*step);

					     // at the moment do this for
					     // isotropic face refinement only
	    for (unsigned int child=0;
		 child<GeometryInfo<dim-1>::max_children_per_cell; ++child)
	      for (unsigned int i=0; i<inner_points.size(); ++i)
		constraint_points.push_back (
		  GeometryInfo<dim-1>::child_to_cell_coordinates(inner_points[i], child));
	  }

					 // Now construct relation between
					 // destination (child) and source (mother)
					 // dofs.
	const unsigned int pnts=(fe.degree+1)*(fe.degree+1);
	const std::vector<Polynomials::Polynomial<double> > polynomial_basis=
	  Polynomials::generate_complete_Lagrange_basis(points.get_points());

	const TensorProductPolynomials<dim-1> face_polynomials(polynomial_basis);

	fe.interface_constraints
	  .TableBase<2,double>::reinit (fe.interface_constraints_size());

	for (unsigned int i=0; i<constraint_points.size(); ++i)
	  {
	    const double interval = (double) (fe.degree * 2);
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
						 // in the FiniteElement class
						 // for constraint matrices in 3d.
		mirror[k] = (constraint_point(k) > 0.5);
		if (mirror[k])
		  constraint_point(k) = 1.0 - constraint_point(k);
	      }

	    for (unsigned j=0; j<pnts; ++j)
	      {
		unsigned int indices[2]
		  = { fe.face_index_map[j] % (fe.degree + 1),
		      fe.face_index_map[j] / (fe.degree + 1) };

		for (unsigned int k = 0; k<2; ++k)
		  if (mirror[k])
		    indices[k] = fe.degree - indices[k];

		const unsigned int
		  new_index = indices[1] * (fe.degree + 1) + indices[0];

		fe.interface_constraints(i,j) =
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
		if (std::fabs(fe.interface_constraints(i,j)) < 1e-14)
		  fe.interface_constraints(i,j) = 0;
	      }
	  }
      }
};





template <int dim, int spacedim>
FE_Q<dim,spacedim>::FE_Q (const unsigned int degree)
		:
		FE_Poly<TensorProductPolynomials<dim>, dim, spacedim> (
		  TensorProductPolynomials<dim>(Polynomials::LagrangeEquidistant::generate_complete_basis(degree)),
		  FiniteElementData<dim>(get_dpo_vector(degree),
					 1, degree,
					 FiniteElementData<dim>::H1),
		  std::vector<bool> (1, false),
		  std::vector<std::vector<bool> >(1, std::vector<bool>(1,true))),
		face_index_map(FE_Q_Helper::invert_numbering(face_lexicographic_to_hierarchic_numbering (degree)))
{
  Assert (degree > 0,
          ExcMessage ("This element can only be used for polynomial degrees "
                      "at least zero"));

  std::vector<unsigned int> renumber (this->dofs_per_cell);
  FETools::hierarchic_to_lexicographic_numbering (*this, renumber);
  this->poly_space.set_numbering(renumber);

				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();

				   // reinit constraints
  initialize_constraints ();

				   // Reinit the vectors of restriction and
				   // prolongation matrices to the right sizes
				   // and compute the matrices
  this->reinit_restriction_and_prolongation_matrices();
  initialize_embedding();
  initialize_restriction();

  initialize_quad_dof_index_permutation();
}



template <int dim, int spacedim>
FE_Q<dim,spacedim>::FE_Q (const Quadrature<1> &points)
		:
		FE_Poly<TensorProductPolynomials<dim>, dim, spacedim> (
		  TensorProductPolynomials<dim>(Polynomials::generate_complete_Lagrange_basis(points.get_points())),
		  FiniteElementData<dim>(get_dpo_vector(points.size()-1),
					 1, points.size()-1,
					 FiniteElementData<dim>::H1),
		  std::vector<bool> (1, false),
		  std::vector<std::vector<bool> >(1, std::vector<bool>(1,true))),
		face_index_map(FE_Q_Helper::invert_numbering(face_lexicographic_to_hierarchic_numbering (points.size()-1)))
{
  const int degree = points.size()-1;

  Assert (degree > 0,
          ExcMessage ("This element can only be used for polynomial degrees "
                      "at least zero"));
  Assert (points.point(0)(0) == 0,
	  ExcMessage ("The first support point has to be zero."));
  Assert (points.point(degree)(0) == 1,
	  ExcMessage ("The last support point has to be one."));

  std::vector<unsigned int> renumber (this->dofs_per_cell);
  FETools::hierarchic_to_lexicographic_numbering (*this, renumber);
  this->poly_space.set_numbering(renumber);

				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points (points);
  initialize_unit_face_support_points (points);

				   // reinit constraints
  Implementation::initialize_constraints (points, *this);

				   // Reinit the vectors of restriction and
				   // prolongation matrices to the right sizes
				   // and compute the matrices
  this->reinit_restriction_and_prolongation_matrices();
  initialize_embedding();
  initialize_restriction();

  initialize_quad_dof_index_permutation();
}



template <int dim, int spacedim>
std::string
FE_Q<dim,spacedim>::get_name () const
{
				   // note that the
				   // FETools::get_fe_from_name
				   // function depends on the
				   // particular format of the string
				   // this function returns, so they
				   // have to be kept in synch

  std::ostringstream namebuf;
  bool type = true;
  const unsigned int n_points = this->degree +1;
  std::vector<double> points(n_points);
  const unsigned int dofs_per_cell = this->dofs_per_cell;
  const std::vector<Point<dim> > &unit_support_points = this->unit_support_points;
  unsigned int index = 0;

				   // Decode the support points
				   // in one coordinate direction.
  for (unsigned int j=0;j<dofs_per_cell;j++)
    {
      if ((dim>1) ? (unit_support_points[j](1)==0 &&
	   ((dim>2) ? unit_support_points[j](2)==0: true)) : true)
	{
	  if (index == 0)
	    points[index] = unit_support_points[j](0);
          else if (index == 1)
	    points[n_points-1] = unit_support_points[j](0);
          else
	    points[index-1] = unit_support_points[j](0);

	  index++;
	}
    }
  Assert (index == n_points,
	  ExcMessage ("Could not decode support points in one coordinate direction."));

				  // Check whether the support
				  // points are equidistant.
  for(unsigned int j=0;j<n_points;j++)
    if (std::fabs(points[j] - (double)j/this->degree) > 1e-15)
      {
	type = false;
	break;
      }

  if (type == true)
    namebuf << "FE_Q<" << dim << ">(" << this->degree << ")";
  else
    {

				  // Check whether the support
				  // points come from QGaussLobatto.
      const QGaussLobatto<1> points_gl(n_points);
      type = true;
      for(unsigned int j=0;j<n_points;j++)
	if (points[j] != points_gl.point(j)(0))
	  {
	    type = false;
	    break;
	  }
      if(type == true)
	namebuf << "FE_Q<" << dim << ">(QGaussLobatto(" << this->degree+1 << "))";
      else
	namebuf << "FE_Q<" << dim << ">(QUnknownNodes(" << this->degree << "))";
    }
  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_Q<dim,spacedim>::clone() const
{
  return new FE_Q<dim,spacedim>(*this);
}



template <int dim, int spacedim>
void
FE_Q<dim,spacedim>::
get_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
			  FullMatrix<double>       &interpolation_matrix) const
{
				   // this is only implemented, if the
				   // source FE is also a
				   // Q element
  typedef FE_Q<dim,spacedim> FEQ;
  typedef FiniteElement<dim,spacedim> FEL;

  AssertThrow ((x_source_fe.get_name().find ("FE_Q<") == 0)
               ||
               (dynamic_cast<const FEQ*>(&x_source_fe) != 0),
               typename FEL::
               ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.m() == this->dofs_per_cell,
	  ExcDimensionMismatch (interpolation_matrix.m(),
				this->dofs_per_cell));
  Assert (interpolation_matrix.n() == x_source_fe.dofs_per_cell,
	  ExcDimensionMismatch (interpolation_matrix.m(),
				x_source_fe.dofs_per_cell));

				   // ok, source is a Q element, so
				   // we will be able to do the work
  const FE_Q<dim,spacedim> &source_fe
    = dynamic_cast<const FE_Q<dim,spacedim>&>(x_source_fe);

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
                                   // read in a point on this
                                   // cell and evaluate the
                                   // shape functions there
      const Point<dim> p = this->unit_support_points[j];
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        cell_interpolation(j,i) = this->poly_space.compute_value (i, p);

      for (unsigned int i=0; i<source_fe.dofs_per_cell; ++i)
        source_interpolation(j,i) = source_fe.poly_space.compute_value (i, p);
    }

                                   // then compute the
                                   // interpolation matrix
                                   // for this coordinate
                                   // direction
  cell_interpolation.gauss_jordan ();
  cell_interpolation.mmult (interpolation_matrix,
                            source_interpolation);

  const double eps = 2e-13*this->degree*dim;

                                   // cut off very small values
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    for (unsigned int j=0; j<source_fe.dofs_per_cell; ++j)
      if (std::fabs(interpolation_matrix(i,j)) < eps)
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

      Assert (std::fabs(sum-1) < eps, ExcInternalError());
    }
}



template <>
void
FE_Q<1>::
get_face_interpolation_matrix (const FiniteElement<1> &/*x_source_fe*/,
			       FullMatrix<double>     &/*interpolation_matrix*/) const
{
  Assert (false,
	  FiniteElement<1>::
	  ExcInterpolationNotImplemented ());
}



template <>
void
FE_Q<1>::
get_subface_interpolation_matrix (const FiniteElement<1> &/*x_source_fe*/,
				  const unsigned int      /*subface*/,
				  FullMatrix<double>     &/*interpolation_matrix*/) const
{
  Assert (false,
	  FiniteElement<1>::
	  ExcInterpolationNotImplemented ());
}


template <>
void
FE_Q<1,2>::
get_face_interpolation_matrix (const FiniteElement<1,2> &/*x_source_fe*/,
			       FullMatrix<double>     &/*interpolation_matrix*/) const
{
  typedef FiniteElement<1,2> FEL;
  Assert (false,
	  FEL::
	   ExcInterpolationNotImplemented ());
}


template <>
void
FE_Q<1,2>::
get_subface_interpolation_matrix (const FiniteElement<1,2> &/*x_source_fe*/,
				  const unsigned int      /*subface*/,
				  FullMatrix<double>     &/*interpolation_matrix*/) const
{
  typedef FiniteElement<1,2> FEL;
  Assert (false,
	  FEL::
	  ExcInterpolationNotImplemented ());
}



template <int dim, int spacedim>
void
FE_Q<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
			       FullMatrix<double>       &interpolation_matrix) const
{
				   // this is only implemented, if the
				   // source FE is also a
				   // Q element
  typedef FE_Q<dim,spacedim> FEQ;
  typedef FiniteElement<dim,spacedim> FEL;
  AssertThrow ((x_source_fe.get_name().find ("FE_Q<") == 0)
               ||
               (dynamic_cast<const FEQ*>(&x_source_fe) != 0),
               typename FEL::
               ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.n() == this->dofs_per_face,
	  ExcDimensionMismatch (interpolation_matrix.n(),
				this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
	  ExcDimensionMismatch (interpolation_matrix.m(),
				x_source_fe.dofs_per_face));

				   // ok, source is a Q element, so
				   // we will be able to do the work
  const FE_Q<dim,spacedim> &source_fe
    = dynamic_cast<const FE_Q<dim,spacedim>&>(x_source_fe);

				   // Make sure, that the element,
                                   // for which the DoFs should be
                                   // constrained is the one with
                                   // the higher polynomial degree.
                                   // Actually the procedure will work
                                   // also if this assertion is not
                                   // satisfied. But the matrices
                                   // produced in that case might
				   // lead to problems in the
                                   // hp procedures, which use this
				   // method.
  Assert (this->dofs_per_face <= source_fe.dofs_per_face,
	  typename FEL::
	  ExcInterpolationNotImplemented ());

                                   // generate a quadrature
                                   // with the unit support points.
                                   // This is later based as a
				   // basis for the QProjector,
				   // which returns the support
                                   // points on the face.
  Quadrature<dim-1> quad_face_support (source_fe.get_unit_face_support_points ());

				   // Rule of thumb for FP accuracy,
				   // that can be expected for a
				   // given polynomial degree.
				   // This value is used to cut
				   // off values close to zero.
  const double eps = 2e-13*this->degree*(dim-1);

				   // compute the interpolation
				   // matrix by simply taking the
                                   // value at the support points.
//TODO: Verify that all faces are the same with respect to
// these support points. Furthermore, check if something has to
// be done for the face orientation flag in 3D.
  const Quadrature<dim> face_quadrature
    = QProjector<dim>::project_to_face (quad_face_support, 0);
  for (unsigned int i=0; i<source_fe.dofs_per_face; ++i)
    {
      const Point<dim> &p = face_quadrature.point (i);

      for (unsigned int j=0; j<this->dofs_per_face; ++j)
	{
	  double matrix_entry = this->shape_value (this->face_to_cell_index(j, 0), p);

				   // Correct the interpolated
				   // value. I.e. if it is close
			     	   // to 1 or 0, make it exactly
			       	   // 1 or 0. Unfortunately, this
			       	   // is required to avoid problems
			       	   // with higher order elements.
	  if (std::fabs (matrix_entry - 1.0) < eps)
	    matrix_entry = 1.0;
	  if (std::fabs (matrix_entry) < eps)
	    matrix_entry = 0.0;

	  interpolation_matrix(i,j) = matrix_entry;
	}
    }

				   // make sure that the row sum of
				   // each of the matrices is 1 at
				   // this point. this must be so
				   // since the shape functions sum up
				   // to 1
  for (unsigned int j=0; j<source_fe.dofs_per_face; ++j)
    {
      double sum = 0.;

      for (unsigned int i=0; i<this->dofs_per_face; ++i)
        sum += interpolation_matrix(j,i);

      Assert (std::fabs(sum-1) < 2e-13*this->degree*(dim-1),
              ExcInternalError());
    }
}



template <int dim, int spacedim>
void
FE_Q<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
				  const unsigned int        subface,
				  FullMatrix<double>       &interpolation_matrix) const
{
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
	  ExcDimensionMismatch (interpolation_matrix.m(),
				x_source_fe.dofs_per_face));

				   // see if source is a Q element
  if (const FE_Q<dim,spacedim> *source_fe
      = dynamic_cast<const FE_Q<dim,spacedim> *>(&x_source_fe))
    {
				       // have this test in here since
				       // a table of size 2x0 reports
				       // its size as 0x0
      Assert (interpolation_matrix.n() == this->dofs_per_face,
	      ExcDimensionMismatch (interpolation_matrix.n(),
				    this->dofs_per_face));

				       // Make sure, that the element,
				       // for which the DoFs should be
				       // constrained is the one with
				       // the higher polynomial degree.
				       // Actually the procedure will work
				       // also if this assertion is not
				       // satisfied. But the matrices
				       // produced in that case might
				       // lead to problems in the
				       // hp procedures, which use this
				       // method.
      Assert (this->dofs_per_face <= source_fe->dofs_per_face,
	      (typename FiniteElement<dim,spacedim>::
	       ExcInterpolationNotImplemented ()));

				       // generate a point on this
				       // cell and evaluate the
				       // shape functions there
      const Quadrature<dim-1>
	quad_face_support (source_fe->get_unit_face_support_points ());

				       // Rule of thumb for FP accuracy,
				       // that can be expected for a
				       // given polynomial degree.
				       // This value is used to cut
				       // off values close to zero.
      double eps = 2e-13*this->degree*(dim-1);

				       // compute the interpolation
				       // matrix by simply taking the
				       // value at the support points.
//TODO: Verify that all faces are the same with respect to
// these support points. Furthermore, check if something has to
// be done for the face orientation flag in 3D.
      const Quadrature<dim> subface_quadrature
	= QProjector<dim>::project_to_subface (quad_face_support, 0, subface);
      for (unsigned int i=0; i<source_fe->dofs_per_face; ++i)
	{
	  const Point<dim> &p = subface_quadrature.point (i);

	  for (unsigned int j=0; j<this->dofs_per_face; ++j)
	    {
	      double matrix_entry = this->shape_value (this->face_to_cell_index(j, 0), p);

					       // Correct the interpolated
					       // value. I.e. if it is close
					       // to 1 or 0, make it exactly
					       // 1 or 0. Unfortunately, this
					       // is required to avoid problems
					       // with higher order elements.
	      if (std::fabs (matrix_entry - 1.0) < eps)
		matrix_entry = 1.0;
	      if (std::fabs (matrix_entry) < eps)
		matrix_entry = 0.0;

	      interpolation_matrix(i,j) = matrix_entry;
	    }
	}

				       // make sure that the row sum of
				       // each of the matrices is 1 at
				       // this point. this must be so
				       // since the shape functions sum up
				       // to 1
      for (unsigned int j=0; j<source_fe->dofs_per_face; ++j)
	{
	  double sum = 0.;

	  for (unsigned int i=0; i<this->dofs_per_face; ++i)
	    sum += interpolation_matrix(j,i);

	  Assert (std::fabs(sum-1) < 2e-13*this->degree*dim,
		  ExcInternalError());
	}
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&x_source_fe) != 0)
    {
				       // nothing to do here, the
				       // FE_Nothing has no degrees of
				       // freedom anyway
    }
  else
    AssertThrow (false,
		 (typename FiniteElement<dim,spacedim>::
		  ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
bool
FE_Q<dim,spacedim>::hp_constraints_are_implemented () const
{
  return true;
}




template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Q<dim,spacedim>::
hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
				   // we can presently only compute
				   // these identities if both FEs are
				   // FE_Qs or if the other one is an
				   // FE_Nothing. in the first case,
				   // there should be exactly one
				   // single DoF of each FE at a
				   // vertex, and they should have
				   // identical value
  if (dynamic_cast<const FE_Q<dim,spacedim>*>(&fe_other) != 0)
    {
      return
	std::vector<std::pair<unsigned int, unsigned int> >
	(1, std::make_pair (0U, 0U));
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
				       // the FE_Nothing has no
				       // degrees of freedom, so there
				       // are no equivalencies to be
				       // recorded
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Q<dim,spacedim>::
hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
				   // we can presently only compute
				   // these identities if both FEs are
				   // FE_Qs or if the other one is an
				   // FE_Nothing
  if (const FE_Q<dim,spacedim> *fe_q_other = dynamic_cast<const FE_Q<dim,spacedim>*>(&fe_other))
    {
				       // dofs are located along lines, so two
				       // dofs are identical if they are
				       // located at identical positions. if
				       // we had only equidistant points, we
				       // could simple check for similarity
				       // like (i+1)*q == (j+1)*p, but we
				       // might have other support points
				       // (e.g. Gauss-Lobatto
				       // points). Therefore, read the points
				       // in unit_support_points for the first
				       // coordinate direction. We take the
				       // lexicographic ordering of the points
				       // in the first direction (i.e.,
				       // x-direction), which we access
				       // between index 1 and p-1 (index 0 and
				       // p are vertex dofs).
      const unsigned int p = this->degree;
      const unsigned int q = fe_q_other->degree;

      std::vector<std::pair<unsigned int, unsigned int> > identities;

      const std::vector<unsigned int> &index_map_inverse=
	this->poly_space.get_numbering_inverse();
      const std::vector<unsigned int> &index_map_inverse_other=
	fe_q_other->poly_space.get_numbering_inverse();

      for (unsigned int i=0; i<p-1; ++i)
	for (unsigned int j=0; j<q-1; ++j)
	  if (std::fabs(this->unit_support_points[index_map_inverse[i+1]][0]-
			fe_q_other->unit_support_points[index_map_inverse_other[j+1]][0])
	      < 1e-14)
	    identities.push_back (std::make_pair(i,j));

      return identities;
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
				       // the FE_Nothing has no
				       // degrees of freedom, so there
				       // are no equivalencies to be
				       // recorded
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Q<dim,spacedim>::
hp_quad_dof_identities (const FiniteElement<dim,spacedim>        &fe_other) const
{
				   // we can presently only compute
				   // these identities if both FEs are
				   // FE_Qs or if the other one is an
				   // FE_Nothing
  if (const FE_Q<dim,spacedim> *fe_q_other = dynamic_cast<const FE_Q<dim,spacedim>*>(&fe_other))
    {
				       // this works exactly like the line
				       // case above, except that now we have
				       // to have two indices i1, i2 and j1,
				       // j2 to characterize the dofs on the
				       // face of each of the finite
				       // elements. since they are ordered
				       // lexicographically along the first
				       // line and we have a tensor product,
				       // the rest is rather straightforward
      const unsigned int p = this->degree;
      const unsigned int q = fe_q_other->degree;

      std::vector<std::pair<unsigned int, unsigned int> > identities;

      const std::vector<unsigned int> &index_map_inverse=
	this->poly_space.get_numbering_inverse();
      const std::vector<unsigned int> &index_map_inverse_other=
	fe_q_other->poly_space.get_numbering_inverse();

      for (unsigned int i1=0; i1<p-1; ++i1)
	for (unsigned int i2=0; i2<p-1; ++i2)
	  for (unsigned int j1=0; j1<q-1; ++j1)
	    for (unsigned int j2=0; j2<q-1; ++j2)
	      if ((std::fabs(this->unit_support_points[index_map_inverse[i1+1]][0]-
			     fe_q_other->unit_support_points[index_map_inverse_other[j1+1]][0])
		   < 1e-14)
		  &&
		  (std::fabs(this->unit_support_points[index_map_inverse[i2+1]][0]-
			     fe_q_other->unit_support_points[index_map_inverse_other[j2+1]][0])
		   < 1e-14))
		identities.push_back (std::make_pair(i1*(p-1)+i2,
						     j1*(q-1)+j2));

      return identities;
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
				       // the FE_Nothing has no
				       // degrees of freedom, so there
				       // are no equivalencies to be
				       // recorded
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Q<dim,spacedim>::
compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_Q<dim,spacedim> *fe_q_other
      = dynamic_cast<const FE_Q<dim,spacedim>*>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
	return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
	return FiniteElementDomination::either_element_can_dominate;
      else
	return FiniteElementDomination::other_element_dominates;
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
				       // the FE_Nothing has no
				       // degrees of freedom and it is
				       // typically used in a context
				       // where we don't require any
				       // continuity along the
				       // interface
      return FiniteElementDomination::no_requirements;
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------



template <int dim, int spacedim>
void FE_Q<dim,spacedim>::initialize_unit_support_points ()
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



template <int dim, int spacedim>
void FE_Q<dim,spacedim>::initialize_unit_support_points (const Quadrature<1> &points)
{
				   // number of points: (degree+1)^dim
  unsigned int n = this->degree+1;
  for (unsigned int i=1; i<dim; ++i)
    n *= this->degree+1;

  this->unit_support_points.resize(n);

  const std::vector<unsigned int> &index_map_inverse=
    this->poly_space.get_numbering_inverse();

  Quadrature<dim> support_quadrature(points);

  Point<dim> p;

  for (unsigned int k=0;k<n ;k++)
    {
      this->unit_support_points[index_map_inverse[k]] = support_quadrature.point(k);
    }
}



template <>
void FE_Q<1>::initialize_unit_face_support_points ()
{
				   // no faces in 1d, so nothing to do
}

template <>
void FE_Q<1>::initialize_unit_face_support_points (const Quadrature<1> &/*points*/)
{
				   // no faces in 1d, so nothing to do
}

template <>
void FE_Q<1,2>::initialize_unit_face_support_points ()
{
				   // no faces in 1d, so nothing to do
}

template <>
void FE_Q<1,2>::initialize_unit_face_support_points (const Quadrature<1> &/*points*/)
{
				   // no faces in 1d, so nothing to do
}



template <int dim, int spacedim>
void FE_Q<dim,spacedim>::initialize_unit_face_support_points ()
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



template <int dim, int spacedim>
void FE_Q<dim,spacedim>::initialize_unit_face_support_points (const Quadrature<1> &points)
{
  const unsigned int codim = dim-1;

				   // number of points: (degree+1)^codim
  unsigned int n = this->degree+1;
  for (unsigned int i=1; i<codim; ++i)
    n *= this->degree+1;

  this->unit_face_support_points.resize(n);

  const std::vector< Point<1> > edge = points.get_points();

  const std::vector<unsigned int> &face_index_map_inverse=
    FE_Q_Helper::invert_numbering(face_index_map);

  Point<codim> p;

  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((codim>2) ? this->degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((codim>1) ? this->degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=this->degree; ++ix)
	{
	  p(0) = edge[ix](0);
	  if (codim>1)
	    p(1) = edge[iy](0);
	  if (codim>2)
	    p(2) = edge[iz](0);

	  this->unit_face_support_points[face_index_map_inverse[k++]] = p;
	}
}



template <int dim, int spacedim>
void
FE_Q<dim,spacedim>::initialize_quad_dof_index_permutation ()
{
				   // general template for 1D and 2D, do nothing
}



template <>
void
FE_Q<3>::initialize_quad_dof_index_permutation ()
{
  Assert (adjust_quad_dof_index_for_face_orientation_table.n_elements()==8*this->dofs_per_quad,
	  ExcInternalError());

  const unsigned int n=this->degree-1;
  Assert(n*n==this->dofs_per_quad, ExcInternalError());

				   // alias for the table to fill
  Table<2,int> &data=this->adjust_quad_dof_index_for_face_orientation_table;

				   // the dofs on a face are connected to a n x
				   // n matrix. for example, for degree==4 we
				   // have the following dofs on a quad

				   //  ___________
				   // |           |
				   // |  6  7  8  |
				   // |           |
				   // |  3  4  5  |
				   // |           |
				   // |  0  1  2  |
				   // |___________|
				   //
				   // we have dof_no=i+n*j with index i in
				   // x-direction and index j in y-direction
				   // running from 0 to n-1.  to extract i and j
				   // we can use i=dof_no%n and j=dof_no/n. i
				   // and j can then be used to construct the
				   // rotated and mirrored numbers.


  for (unsigned int local=0; local<this->dofs_per_quad; ++local)
				     // face support points are in lexicographic
				     // ordering with x running fastest. invert
				     // that (y running fastest)
    {
      unsigned int i=local%n,
		   j=local/n;

				       // face_orientation=false, face_flip=false, face_rotation=false
      data(local,0)=j       + i      *n - local;
				       // face_orientation=false, face_flip=false, face_rotation=true
      data(local,1)=i       + (n-1-j)*n - local;
				       // face_orientation=false, face_flip=true,  face_rotation=false
      data(local,2)=(n-1-j) + (n-1-i)*n - local;
				       // face_orientation=false, face_flip=true,  face_rotation=true
      data(local,3)=(n-1-i) + j      *n - local;
				       // face_orientation=true,  face_flip=false, face_rotation=false
      data(local,4)=0;
				       // face_orientation=true,  face_flip=false, face_rotation=true
      data(local,5)=j       + (n-1-i)*n - local;
				       // face_orientation=true,  face_flip=true,  face_rotation=false
      data(local,6)=(n-1-i) + (n-1-j)*n - local;
				       // face_orientation=true,  face_flip=true,  face_rotation=true
      data(local,7)=(n-1-j) + i      *n - local;
    }

				   // aditionally initialize reordering of line
				   // dofs
  for (unsigned int i=0; i<this->dofs_per_line; ++i)
    this->adjust_line_dof_index_for_line_orientation_table[i]=this->dofs_per_line-1-i - i;
}




template <int dim, int spacedim>
std::vector<unsigned int>
FE_Q<dim,spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}


template <int dim, int spacedim>
std::vector<unsigned int>
FE_Q<dim,spacedim>::face_lexicographic_to_hierarchic_numbering (const unsigned int degree)
{
  const FiniteElementData<dim-1> face_data(FE_Q<dim-1>::get_dpo_vector(degree),1,degree);
  std::vector<unsigned int> face_renumber (face_data.dofs_per_cell);
  FETools::lexicographic_to_hierarchic_numbering (face_data, face_renumber);
  return face_renumber;
}



template <>
std::vector<unsigned int>
FE_Q<1>::face_lexicographic_to_hierarchic_numbering (const unsigned int)
{
  return std::vector<unsigned int>();
}

template <>
std::vector<unsigned int>
FE_Q<1,2>::face_lexicographic_to_hierarchic_numbering (const unsigned int)
{
  return std::vector<unsigned int>();
}



template <int dim, int spacedim>
void
FE_Q<dim,spacedim>::initialize_constraints ()
{
  QTrapez<1> trapez;
  QIterated<1> points (trapez,this->degree);
  Implementation::initialize_constraints (points, *this);
}



template <int dim, int spacedim>
void
FE_Q<dim,spacedim>::initialize_embedding ()
{
				   // compute the interpolation matrices in
				   // much the same way as we do for the
				   // constraints. it's actually simpler
				   // here, since we don't have this weird
				   // renumbering stuff going on. The trick
				   // is again that we the interpolation
				   // matrix is formed by a permutation of
				   // the indices of the cell matrix. The
				   // value eps is used a threshold to
				   // decide when certain evaluations of the
				   // Lagrange polynomials are zero or one.
  const double eps = 1e-15*this->degree*dim;

#ifdef DEBUG
				// in DEBUG mode, check that the evaluation of
				// support points in the current numbering
				// gives the identity operation
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    {
      Assert (std::fabs (1.-this->poly_space.compute_value
			 (i, this->unit_support_points[i])) < eps,
	      ExcInternalError());
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
	if (j!=i)
	  Assert (std::fabs (this->poly_space.compute_value
			     (i, this->unit_support_points[j])) < eps,
		  ExcInternalError());
    }
#endif

				// to efficiently evaluate the polynomial at
				// the subcell, make use of the tensor product
				// structure of this element and only evaluate
				// 1D information from the polynomial. This
				// makes the cost of this function almost
				// negligible also for high order elements
  const unsigned int dofs1d = this->degree+1;
  std::vector<Table<2,double> >
    subcell_evaluations (dim, Table<2,double>(dofs1d, dofs1d));
  const std::vector<unsigned int> &index_map_inverse =
    this->poly_space.get_numbering_inverse();

				//  recreate 1D polynomials
  std::vector<Polynomials::Polynomial<double> > poly_space1d =
    FE_Q_Helper::generate_poly_space1d (this->unit_support_points,
                                        index_map_inverse, dofs1d);

				// helper value: step size how to walk through
				// diagonal and how many points we have left
				// apart from the first dimension
  unsigned int step_size_diag = 0;
  {
    unsigned int factor = 1;
    for (unsigned int d=0; d<dim; ++d)
      {
	step_size_diag += factor;
	factor *= dofs1d;
      }
  }

					     // next evaluate the functions
					     // for the different refinement
					     // cases.
  for (unsigned int ref=0; ref<RefinementCase<dim>::isotropic_refinement; ++ref)
    for (unsigned int child=0; child<GeometryInfo<dim>::n_children(RefinementCase<dim>(ref+1)); ++child)
      {
				// go through the points in diagonal to
				// capture variation in all directions
				// simultaneously
	for (unsigned int j=0; j<dofs1d; ++j)
	  {
	    const unsigned int diag_comp = index_map_inverse[j*step_size_diag];
	    const Point<dim> p_subcell = this->unit_support_points[diag_comp];
	    const Point<dim> p_cell =
	      GeometryInfo<dim>::child_to_cell_coordinates (p_subcell, child,
							    RefinementCase<dim>(ref+1));
	    for (unsigned int i=0; i<dofs1d; ++i)
	      for (unsigned int d=0; d<dim; ++d)
		{
		  const double cell_value = poly_space1d[i].value (p_cell[d]);

				   // cut off values that are too
				   // small. note that we have here Lagrange
				   // interpolation functions, so they
				   // should be zero at almost all points,
				   // and one at the others, at least on the
				   // subcells. so set them to their exact
				   // values
				   //
				   // the actual cut-off value is somewhat
				   // fuzzy, but it works for
				   // 2e-13*degree*dim (see above), which
				   // is kind of reasonable given that we
				   // compute the values of the polynomials
				   // via an degree-step recursion and then
				   // multiply the 1d-values. this gives us
				   // a linear growth in degree*dim, times a
				   // small constant.
				   //
				   // the embedding matrix is given by
				   // applying the inverse of the subcell
				   // matrix on the cell_interpolation
				   // matrix. since the subcell matrix is
				   // actually only a permutation vector,
				   // all we need to do is to switch the
				   // rows we write the data into. moreover,
				   // cut off very small values here
		  if (std::fabs(cell_value) < eps)
		    subcell_evaluations[d](j,i) = 0;
		  else
		    subcell_evaluations[d](j,i) = cell_value;
		}
	  }

				// now expand from 1D info. block innermost
				// dimension (x_0) in order to avoid difficult
				// checks at innermost loop
	unsigned int j_indices[dim];
	FE_Q_Helper::zero_indices<dim> (j_indices);
	for (unsigned int j=0; j<this->dofs_per_cell; j+=dofs1d)
	  {
	    unsigned int i_indices[dim];
	    FE_Q_Helper::zero_indices<dim> (i_indices);
	    for (unsigned int i=0; i<this->dofs_per_cell; i+=dofs1d)
	      {
		double val_extra_dim = 1.;
		for (unsigned int d=1; d<dim; ++d)
		  val_extra_dim *= subcell_evaluations[d](j_indices[d-1],
							  i_indices[d-1]);

				// innermost sum where we actually
				// compute. the same as
				// this->prolongation[ref][child](j,i) =
				// this->poly_space.compute_value (i, p_cell);
		for (unsigned int jj=0; jj<dofs1d; ++jj)
		  {
		    const unsigned int j_ind = index_map_inverse[j+jj];
		    for (unsigned int ii=0; ii<dofs1d; ++ii)
		      this->prolongation[ref][child](j_ind,index_map_inverse[i+ii])
			= val_extra_dim * subcell_evaluations[0](jj,ii);
		  }

				// update indices that denote the tensor
				// product position. a bit fuzzy and therefore
				// not done for innermost x_0 direction
		FE_Q_Helper::increment_indices<dim> (i_indices, dofs1d);
	      }
	    Assert (i_indices[dim-1] == 1, ExcInternalError());
	    FE_Q_Helper::increment_indices<dim> (j_indices, dofs1d);
	  }

				   // and make sure that the row sum is
				   // 1. this must be so since for this
				   // element, the shape functions add up to
				   // one
#ifdef DEBUG
	for (unsigned int row=0; row<this->dofs_per_cell; ++row)
	  {
	    double sum = 0;
	    for (unsigned int col=0; col<this->dofs_per_cell; ++col)
	      sum += this->prolongation[ref][child](row,col);
	    Assert (std::fabs(sum-1.) < eps, ExcInternalError());
	  }
#endif
      }
}



template <int dim, int spacedim>
void
FE_Q<dim,spacedim>::initialize_restriction ()
{
                                   // for these Lagrange interpolation
                                   // polynomials, construction of the
                                   // restriction matrices is relatively
                                   // simple. the reason is that the
                                   // interpolation points on the mother cell
                                   // are (except for the case with arbitrary
                                   // nonequidistant nodes) always also
                                   // interpolation points for some shape
                                   // function on one or the other child,
                                   // because we have chosen equidistant
                                   // Lagrange interpolation points for the
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
                                   // interpolation points there
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

  const double eps = 1e-15*this->degree*dim;
  const std::vector<unsigned int> &index_map_inverse =
    this->poly_space.get_numbering_inverse();

				//  recreate 1D polynomials for faster
				//  evaluation of polynomial
  const unsigned int dofs1d = this->degree+1;
  std::vector<Polynomials::Polynomial<double> > poly_space1d =
    FE_Q_Helper::generate_poly_space1d (this->unit_support_points,
                                        index_map_inverse, dofs1d);
  std::vector<Tensor<1,dim> > evaluations1d (dofs1d);

  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    {
      unsigned int mother_dof = index_map_inverse[i];
      const Point<dim> p_cell = this->unit_support_points[mother_dof];

                                       // then find the children on
                                       // which the interpolation
                                       // point is located
      for (unsigned int ref=RefinementCase<dim>::cut_x; ref<=RefinementCase<dim>::isotropic_refinement; ++ref)
	for (unsigned int child=0; child<GeometryInfo<dim>::n_children(RefinementCase<dim>(ref)); ++child)
	  {
					     // check whether this
					     // interpolation point is
					     // inside this child cell
	    const Point<dim> p_subcell
	      = GeometryInfo<dim>::cell_to_child_coordinates (p_cell, child, RefinementCase<dim>(ref));
	    if (GeometryInfo<dim>::is_inside_unit_cell (p_subcell))
	      {
				// same logic as in initialize_embedding to
				// evaluate the polynomial faster than from
				// the tensor product: since we evaluate all
				// polynomials, it is much faster to just
				// compute the 1D values for all polynomials
				// before and then get the dim-data.
		for (unsigned int j=0; j<dofs1d; ++j)
		  for (unsigned int d=0; d<dim; ++d)
		    evaluations1d[j][d] = poly_space1d[j].value (p_subcell[d]);
		unsigned int j_indices[dim];
		FE_Q_Helper::zero_indices<dim> (j_indices);
		double sum_check = 0;
		for (unsigned int j = 0; j<this->dofs_per_cell; j += dofs1d)
		  {
		    double val_extra_dim = 1.;
		    for (unsigned int d=1; d<dim; ++d)
		      val_extra_dim *= evaluations1d[j_indices[d-1]][d];
		    for (unsigned int jj=0; jj<dofs1d; ++jj)
		      {

						 // find the child shape
						 // function(s) corresponding
						 // to this point. Usually
						 // this is just one function;
						 // however, when we use FE_Q
						 // on arbitrary nodes a
						 // parent support point might
						 // not be a child support
						 // point, and then we will
						 // get more than one nonzero
						 // value per row. Still, the
						 // values should sum up to 1
			const double val
			  = val_extra_dim * evaluations1d[jj][0];
			const unsigned int child_dof =
			  index_map_inverse[j+jj];
			if (std::fabs (val-1.) < eps)
			  this->restriction[ref-1][child](mother_dof,child_dof)=1.;
			else if(std::fabs(val) > eps)
			  this->restriction[ref-1][child](mother_dof,child_dof)=val;
			sum_check += val;
		      }
		    FE_Q_Helper::increment_indices<dim> (j_indices, dofs1d);
		  }
		Assert (std::fabs(sum_check-1) < eps,
			ExcInternalError());
	      }
	  }
    }
}



//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------


template <>
bool
FE_Q<1>::has_support_on_face (const unsigned int shape_index,
                              const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<1>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<1>::faces_per_cell));


				   // in 1d, things are simple. since
				   // there is only one degree of
				   // freedom per vertex in this
				   // class, the first is on vertex 0
				   // (==face 0 in some sense), the
				   // second on face 1:
  return (((shape_index == 0) && (face_index == 0)) ||
          ((shape_index == 1) && (face_index == 1)));
}


template <>
bool
FE_Q<1,2>::has_support_on_face (const unsigned int shape_index,
                              const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<1>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<1>::faces_per_cell));


				   // in 1d, things are simple. since
				   // there is only one degree of
				   // freedom per vertex in this
				   // class, the first is on vertex 0
				   // (==face 0 in some sense), the
				   // second on face 1:
  return (((shape_index == 0) && (face_index == 0)) ||
          ((shape_index == 1) && (face_index == 1)));
}



template <int dim, int spacedim>
bool
FE_Q<dim,spacedim>::has_support_on_face (const unsigned int shape_index,
				const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));


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

                                   // we should not have gotten here
  Assert (false, ExcInternalError());
  return false;
}



template <int dim, int spacedim>
std::size_t
FE_Q<dim,spacedim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



// explicit instantiations
#include "fe_q.inst"

DEAL_II_NAMESPACE_CLOSE
