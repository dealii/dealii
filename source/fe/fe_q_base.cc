// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/tensor_product_polynomials_const.h>
#include <deal.II/base/tensor_product_polynomials_bubbles.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/fe/fe_q_base.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/base/quadrature_lib.h>

#include <vector>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


namespace FE_Q_Helper
{
  namespace
  {
    // get the renumbering for faces
    template <int dim>
    inline
    std::vector<unsigned int>
    face_lexicographic_to_hierarchic_numbering (const unsigned int degree)
    {
      std::vector<unsigned int> dpo(dim, 1U);
      for (unsigned int i=1; i<dpo.size(); ++i)
        dpo[i]=dpo[i-1]*(degree-1);
      const dealii::FiniteElementData<dim-1> face_data(dpo,1,degree);
      std::vector<unsigned int> face_renumber (face_data.dofs_per_cell);
      FETools::lexicographic_to_hierarchic_numbering (face_data, face_renumber);
      return face_renumber;
    }

    // dummy specialization for dim == 1 to avoid linker errors
    template <>
    inline
    std::vector<unsigned int>
    face_lexicographic_to_hierarchic_numbering<1> (const unsigned int)
    {
      return std::vector<unsigned int>();
    }



    // in get_restriction_matrix() and get_prolongation_matrix(), want to undo
    // tensorization on inner loops for performance reasons. this clears a
    // dim-array
    template <int dim>
    inline
    void
    zero_indices (unsigned int (&indices)[dim])
    {
      for (unsigned int d=0; d<dim; ++d)
        indices[d] = 0;
    }



    // in get_restriction_matrix() and get_prolongation_matrix(), want to undo
    // tensorization on inner loops for performance reasons. this increments
    // tensor product indices
    template <int dim>
    inline
    void
    increment_indices (unsigned int       (&indices)[dim],
                       const unsigned int   dofs1d)
    {
      ++indices[0];
      for (int d=0; d<dim-1; ++d)
        if (indices[d]==dofs1d)
          {
            indices[d] = 0;
            indices[d+1]++;
          }
    }
  }
}



/**
 * A class with the same purpose as the similarly named class of the
 * Triangulation class. See there for more information.
 */
template <class PolynomialType, int xdim, int xspacedim>
struct FE_Q_Base<PolynomialType,xdim,xspacedim>::Implementation
{
  /**
   * Initialize the hanging node constraints matrices. Called from the
   * constructor in case the finite element is based on quadrature points.
   */
  template <int spacedim>
  static
  void initialize_constraints (const std::vector<Point<1> > &,
                               FE_Q_Base<PolynomialType,1,spacedim> &)
  {
    // no constraints in 1d
  }


  template <int spacedim>
  static
  void initialize_constraints (const std::vector<Point<1> > &/*points*/,
                               FE_Q_Base<PolynomialType,2,spacedim> &fe)
  {
    const unsigned int dim = 2;

    unsigned int q_deg = fe.degree;
    if (types_are_equal<PolynomialType, TensorProductPolynomialsBubbles<dim> >::value)
      q_deg = fe.degree-1;

    // restricted to each face, the traces of the shape functions is an
    // element of P_{k} (in 2d), or Q_{k} (in 3d), where k is the degree of
    // the element.  from this, we interpolate between mother and cell face.

    // the interpolation process works as follows: on each subface, we want
    // that finite element solutions from both sides coincide. i.e. if a and b
    // are expansion coefficients for the shape functions from both sides, we
    // seek a relation between a and b such that
    //   sum_j a_j phi^c_j(x) == sum_j b_j phi_j(x)
    // for all points x on the interface. here, phi^c_j are the shape
    // functions on the small cell on one side of the face, and phi_j those on
    // the big cell on the other side. To get this relation, it suffices to
    // look at a sufficient number of points for which this has to hold. if
    // there are n functions, then we need n evaluation points, and we choose
    // them equidistantly.
    //
    // we obtain the matrix system
    //    A a  ==  B b
    // where
    //    A_ij = phi^c_j(x_i)
    //    B_ij = phi_j(x_i)
    // and the relation we are looking for is
    //    a = A^-1 B b
    //
    // for the special case of Lagrange interpolation polynomials, A_ij
    // reduces to delta_ij, and
    //    a_i = B_ij b_j
    // Hence, interface_constraints(i,j)=B_ij.
    //
    // for the general case, where we don't have Lagrange interpolation
    // polynomials, this is a little more complicated. Then we would evaluate
    // at a number of points and invert the interpolation matrix A.
    //
    // Note, that we build up these matrices for all subfaces at once, rather
    // than considering them separately. the reason is that we finally will
    // want to have them in this order anyway, as this is the format we need
    // inside deal.II

    // In the following the points x_i are constructed in following order
    // (n=degree-1)
    // *----------*---------*
    //     1..n   0  n+1..2n
    // i.e. first the midpoint of the line, then the support points on subface
    // 0 and on subface 1
    std::vector<Point<dim-1> > constraint_points;
    // Add midpoint
    constraint_points.push_back (Point<dim-1> (0.5));

    if (q_deg>1)
      {
        const unsigned int n=q_deg-1;
        const double step=1./q_deg;
        // subface 0
        for (unsigned int i=1; i<=n; ++i)
          constraint_points.push_back (
            GeometryInfo<dim-1>::child_to_cell_coordinates(Point<dim-1>(i*step),0));
        // subface 1
        for (unsigned int i=1; i<=n; ++i)
          constraint_points.push_back (
            GeometryInfo<dim-1>::child_to_cell_coordinates(Point<dim-1>(i*step),1));
      }

    // Now construct relation between destination (child) and source (mother)
    // dofs.

    fe.interface_constraints
    .TableBase<2,double>::reinit (fe.interface_constraints_size());

    // use that the element evaluates to 1 at index 0 and along the line at
    // zero
    const std::vector<unsigned int> &index_map_inverse =
      fe.poly_space.get_numbering_inverse();
    const std::vector<unsigned int> face_index_map =
      FE_Q_Helper::face_lexicographic_to_hierarchic_numbering<dim>(q_deg);
    Assert(std::abs(fe.poly_space.compute_value(index_map_inverse[0],Point<dim>())
                    - 1.) < 1e-14,
           ExcInternalError());

    for (unsigned int i=0; i<constraint_points.size(); ++i)
      for (unsigned int j=0; j<q_deg+1; ++j)
        {
          Point<dim> p;
          p[0] = constraint_points[i](0);
          fe.interface_constraints(i,face_index_map[j]) =
            fe.poly_space.compute_value(index_map_inverse[j], p);

          // if the value is small up to round-off, then simply set it to zero
          // to avoid unwanted fill-in of the constraint matrices (which would
          // then increase the number of other DoFs a constrained DoF would
          // couple to)
          if (std::fabs(fe.interface_constraints(i,face_index_map[j])) < 1e-13)
            fe.interface_constraints(i,face_index_map[j]) = 0;
        }
  }


  template <int spacedim>
  static
  void initialize_constraints (const std::vector<Point<1> > &/*points*/,
                               FE_Q_Base<PolynomialType,3,spacedim> &fe)
  {
    const unsigned int dim = 3;

    unsigned int q_deg = fe.degree;
    if (types_are_equal<PolynomialType,TensorProductPolynomialsBubbles<dim> >::value)
      q_deg = fe.degree-1;

    // For a detailed documentation of the interpolation see the
    // FE_Q_Base<2>::initialize_constraints function.

    // In the following the points x_i are constructed in the order as
    // described in the documentation of the FiniteElement class (fe_base.h),
    // i.e.
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

    // Add midpoints of lines of "mother-face"
    constraint_points.push_back (Point<dim-1> (0, 0.5));
    constraint_points.push_back (Point<dim-1> (1, 0.5));
    constraint_points.push_back (Point<dim-1> (0.5, 0));
    constraint_points.push_back (Point<dim-1> (0.5, 1));

    if (q_deg>1)
      {
        const unsigned int n=q_deg-1;
        const double step=1./q_deg;
        std::vector<Point<dim-2> > line_support_points(n);
        for (unsigned int i=0; i<n; ++i)
          line_support_points[i](0)=(i+1)*step;
        Quadrature<dim-2> qline(line_support_points);

        // auxiliary points in 2d
        std::vector<Point<dim-1> > p_line(n);

        // Add nodes of lines interior in the "mother-face"

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

        // DoFs on bordering lines lines 9-16
        for (unsigned int face=0; face<GeometryInfo<dim-1>::faces_per_cell; ++face)
          for (unsigned int subface=0;
               subface<GeometryInfo<dim-1>::max_children_per_face; ++subface)
            {
              QProjector<dim-1>::project_to_subface(qline, face, subface, p_line);
              constraint_points.insert(constraint_points.end(),
                                       p_line.begin(), p_line.end());
            }

        // Create constraints for interior nodes
        std::vector<Point<dim-1> > inner_points(n*n);
        for (unsigned int i=0, iy=1; iy<=n; ++iy)
          for (unsigned int ix=1; ix<=n; ++ix)
            inner_points[i++] = Point<dim-1> (ix*step, iy*step);

        // at the moment do this for isotropic face refinement only
        for (unsigned int child=0;
             child<GeometryInfo<dim-1>::max_children_per_cell; ++child)
          for (unsigned int i=0; i<inner_points.size(); ++i)
            constraint_points.push_back (
              GeometryInfo<dim-1>::child_to_cell_coordinates(inner_points[i], child));
      }

    // Now construct relation between destination (child) and source (mother)
    // dofs.
    const unsigned int pnts=(q_deg+1)*(q_deg+1);

    // use that the element evaluates to 1 at index 0 and along the line at
    // zero
    const std::vector<unsigned int> &index_map_inverse =
      fe.poly_space.get_numbering_inverse();
    const std::vector<unsigned int> face_index_map =
      FE_Q_Helper::face_lexicographic_to_hierarchic_numbering<dim>(q_deg);
    Assert(std::abs(fe.poly_space.compute_value(index_map_inverse[0],Point<dim>())
                    - 1.) < 1e-14,
           ExcInternalError());

    fe.interface_constraints
    .TableBase<2,double>::reinit (fe.interface_constraints_size());

    for (unsigned int i=0; i<constraint_points.size(); ++i)
      {
        const double interval = (double) (q_deg * 2);
        bool mirror[dim - 1];
        Point<dim> constraint_point;

        // Eliminate FP errors in constraint points. Due to their origin, they
        // must all be fractions of the unit interval. If we have polynomial
        // degree 4, the refined element has 8 intervals.  Hence the
        // coordinates must be 0, 0.125, 0.25, 0.375 etc.  Now the coordinates
        // of the constraint points will be multiplied by the inverse of the
        // interval size (in the example by 8).  After that the coordinates
        // must be integral numbers. Hence a normal truncation is performed
        // and the coordinates will be scaled back. The equal treatment of all
        // coordinates should eliminate any FP errors.
        for (unsigned int k=0; k<dim-1; ++k)
          {
            const int coord_int =
              static_cast<int> (constraint_points[i](k) * interval + 0.25);
            constraint_point(k) = 1.*coord_int / interval;

            // The following lines of code should eliminate the problems with
            // the Constraint-Matrix, which appeared for P>=4. The
            // ConstraintMatrix class complained about different constraints
            // for the same entry of the Constraint-Matrix.  Actually this
            // difference could be attributed to FP errors, as it was in the
            // range of 1.0e-16. These errors originate in the loss of
            // symmetry in the FP approximation of the shape-functions.
            // Considering a 3rd order shape function in 1D, we have
            // N0(x)=N3(1-x) and N1(x)=N2(1-x).  For higher order polynomials
            // the FP approximations of the shape functions do not satisfy
            // these equations any more!  Thus in the following code
            // everything is computed in the interval x \in [0..0.5], which is
            // sufficient to express all values that could come out from a
            // computation of any shape function in the full interval
            // [0..1]. If x > 0.5 the computation is done for 1-x with the
            // shape function N_{p-n} instead of N_n.  Hence symmetry is
            // preserved and everything works fine...
            //
            // For a different explanation of the problem, see the discussion
            // in the FiniteElement class for constraint matrices in 3d.
            mirror[k] = (constraint_point(k) > 0.5);
            if (mirror[k])
              constraint_point(k) = 1.0 - constraint_point(k);
          }

        for (unsigned int j=0; j<pnts; ++j)
          {
            unsigned int indices[2] = { j % (q_deg+1), j / (q_deg+1) };

            for (unsigned int k = 0; k<2; ++k)
              if (mirror[k])
                indices[k] = q_deg - indices[k];

            const unsigned int
            new_index = indices[1] * (q_deg + 1) + indices[0];

            fe.interface_constraints(i,face_index_map[j]) =
              fe.poly_space.compute_value (index_map_inverse[new_index],
                                           constraint_point);

            // if the value is small up to round-off, then simply set it to
            // zero to avoid unwanted fill-in of the constraint matrices
            // (which would then increase the number of other DoFs a
            // constrained DoF would couple to)
            if (std::fabs(fe.interface_constraints(i,face_index_map[j])) < 1e-13)
              fe.interface_constraints(i,face_index_map[j]) = 0;
          }
      }
  }
};



template <class PolynomialType, int dim, int spacedim>
FE_Q_Base<PolynomialType,dim,spacedim>::FE_Q_Base
(const PolynomialType         &poly_space,
 const FiniteElementData<dim> &fe_data,
 const std::vector<bool>      &restriction_is_additive_flags)
  :
  FE_Poly<PolynomialType,dim,spacedim>(poly_space, fe_data, restriction_is_additive_flags,
                                       std::vector<ComponentMask>(1, std::vector<bool>(1,true))),
  q_degree (types_are_equal<PolynomialType, TensorProductPolynomialsBubbles<dim> >::value
            ?this->degree-1
            :this->degree)
{}



template <class PolynomialType, int dim, int spacedim>
void
FE_Q_Base<PolynomialType,dim,spacedim>::initialize (const std::vector<Point<1> > &points)
{
  Assert (points[0][0] == 0,
          ExcMessage ("The first support point has to be zero."));
  Assert (points.back()[0] == 1,
          ExcMessage ("The last support point has to be one."));

  // distinguish q/q_dg0 case: need to be flexible enough to allow more
  // degrees of freedom than there are FE_Q degrees of freedom for derived
  // class FE_Q_DG0 that otherwise shares 95% of the code.
  const unsigned int q_dofs_per_cell = Utilities::fixed_power<dim>(q_degree+1);
  Assert(q_dofs_per_cell == this->dofs_per_cell ||
         q_dofs_per_cell+1 == this->dofs_per_cell ||
         q_dofs_per_cell+dim == this->dofs_per_cell , ExcInternalError());

  {
    std::vector<unsigned int> renumber(q_dofs_per_cell);
    const FiniteElementData<dim> fe(get_dpo_vector(q_degree),1,
                                    q_degree);
    FETools::hierarchic_to_lexicographic_numbering (fe, renumber);
    for (unsigned int i= q_dofs_per_cell; i<this->dofs_per_cell; ++i)
      renumber.push_back(i);
    this->poly_space.set_numbering(renumber);
  }

  // finally fill in support points on cell and face
  initialize_unit_support_points (points);
  initialize_unit_face_support_points (points);

  // reinit constraints
  initialize_constraints (points);

  // do not initialize embedding and restriction here. these matrices are
  // initialized on demand in get_restriction_matrix and
  // get_prolongation_matrix

  this->initialize_quad_dof_index_permutation();
}



template <class PolynomialType, int dim, int spacedim>
void
FE_Q_Base<PolynomialType,dim,spacedim>::
get_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
                          FullMatrix<double>                &interpolation_matrix) const
{
  // go through the list of elements we can interpolate from
  if (const FE_Q_Base<PolynomialType,dim,spacedim> *source_fe
      = dynamic_cast<const FE_Q_Base<PolynomialType,dim,spacedim>*>(&x_source_fe))
    {
      // ok, source is a Q element, so we will be able to do the work
      Assert (interpolation_matrix.m() == this->dofs_per_cell,
              ExcDimensionMismatch (interpolation_matrix.m(),
                                    this->dofs_per_cell));
      Assert (interpolation_matrix.n() == x_source_fe.dofs_per_cell,
              ExcDimensionMismatch (interpolation_matrix.m(),
                                    x_source_fe.dofs_per_cell));

      // only evaluate Q dofs
      const unsigned int q_dofs_per_cell = Utilities::fixed_power<dim>(q_degree+1);
      const unsigned int source_q_dofs_per_cell = Utilities::fixed_power<dim>(source_fe->degree+1);

      // evaluation is simply done by evaluating the other FE's basis functions on
      // the unit support points (FE_Q has the property that the cell
      // interpolation matrix is a unit matrix, so no need to evaluate it and
      // invert it)
      for (unsigned int j=0; j<q_dofs_per_cell; ++j)
        {
          // read in a point on this cell and evaluate the shape functions there
          const Point<dim> p = this->unit_support_points[j];

          // FE_Q element evaluates to 1 in unit support point and to zero in all
          // other points by construction
          Assert(std::abs(this->poly_space.compute_value (j, p)-1.)<1e-13,
                 ExcInternalError());

          for (unsigned int i=0; i<source_q_dofs_per_cell; ++i)
            interpolation_matrix(j,i) = source_fe->poly_space.compute_value (i, p);
        }

      // for FE_Q_DG0, add one last row of identity
      if (q_dofs_per_cell < this->dofs_per_cell)
        {
          AssertDimension(source_q_dofs_per_cell+1, source_fe->dofs_per_cell);
          for (unsigned int i=0; i<source_q_dofs_per_cell; ++i)
            interpolation_matrix(q_dofs_per_cell, i) = 0.;
          for (unsigned int j=0; j<q_dofs_per_cell; ++j)
            interpolation_matrix(j, source_q_dofs_per_cell) = 0.;
          interpolation_matrix(q_dofs_per_cell, source_q_dofs_per_cell) = 1.;
        }

      // cut off very small values
      const double eps = 2e-13*q_degree*dim;
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        for (unsigned int j=0; j<source_fe->dofs_per_cell; ++j)
          if (std::fabs(interpolation_matrix(i,j)) < eps)
            interpolation_matrix(i,j) = 0.;

      // make sure that the row sum of each of the matrices is 1 at this
      // point. this must be so since the shape functions sum up to 1
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        {
          double sum = 0.;
          for (unsigned int j=0; j<source_fe->dofs_per_cell; ++j)
            sum += interpolation_matrix(i,j);

          Assert (std::fabs(sum-1) < eps, ExcInternalError());
        }
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&x_source_fe))
    {
      // the element we want to interpolate from is an FE_Nothing. this
      // element represents a function that is constant zero and has no
      // degrees of freedom, so the interpolation is simply a multiplication
      // with a n_dofs x 0 matrix. there is nothing to do here

      // we would like to verify that the number of rows and columns of
      // the matrix equals this->dofs_per_cell and zero. unfortunately,
      // whenever we do FullMatrix::reinit(m,0), it sets both rows and
      // columns to zero, instead of m and zero. thus, only test the
      // number of columns
      Assert (interpolation_matrix.n() == x_source_fe.dofs_per_cell,
              ExcDimensionMismatch (interpolation_matrix.m(),
                                    x_source_fe.dofs_per_cell));

    }
  else
    AssertThrow (false,
                 (typename FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented()));

}



template <class PolynomialType, int dim, int spacedim>
void
FE_Q_Base<PolynomialType,dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source_fe,
                               FullMatrix<double>                &interpolation_matrix) const
{
  Assert (dim > 1, ExcImpossibleInDim(1));
  get_subface_interpolation_matrix (source_fe, numbers::invalid_unsigned_int,
                                    interpolation_matrix);
}



template <class PolynomialType, int dim, int spacedim>
void
FE_Q_Base<PolynomialType,dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
                                  const unsigned int                 subface,
                                  FullMatrix<double>                &interpolation_matrix) const
{
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // see if source is a Q element
  if (const FE_Q_Base<PolynomialType,dim,spacedim> *source_fe
      = dynamic_cast<const FE_Q_Base<PolynomialType,dim,spacedim> *>(&x_source_fe))
    {
      // have this test in here since a table of size 2x0 reports its size as
      // 0x0
      Assert (interpolation_matrix.n() == this->dofs_per_face,
              ExcDimensionMismatch (interpolation_matrix.n(),
                                    this->dofs_per_face));

      // Make sure that the element for which the DoFs should be constrained
      // is the one with the higher polynomial degree.  Actually the procedure
      // will work also if this assertion is not satisfied. But the matrices
      // produced in that case might lead to problems in the hp procedures,
      // which use this method.
      Assert (this->dofs_per_face <= source_fe->dofs_per_face,
              (typename FiniteElement<dim,spacedim>::
               ExcInterpolationNotImplemented ()));

      // generate a point on this cell and evaluate the shape functions there
      const Quadrature<dim-1>
      quad_face_support (source_fe->get_unit_face_support_points ());

      // Rule of thumb for FP accuracy, that can be expected for a given
      // polynomial degree.  This value is used to cut off values close to
      // zero.
      double eps = 2e-13*q_degree*(dim-1);

      // compute the interpolation matrix by simply taking the value at the
      // support points.
//TODO: Verify that all faces are the same with respect to
// these support points. Furthermore, check if something has to
// be done for the face orientation flag in 3D.
      const Quadrature<dim> subface_quadrature
        = subface == numbers::invalid_unsigned_int
          ?
          QProjector<dim>::project_to_face (quad_face_support, 0)
          :
          QProjector<dim>::project_to_subface (quad_face_support, 0, subface);
      for (unsigned int i=0; i<source_fe->dofs_per_face; ++i)
        {
          const Point<dim> &p = subface_quadrature.point (i);

          for (unsigned int j=0; j<this->dofs_per_face; ++j)
            {
              double matrix_entry = this->shape_value (this->face_to_cell_index(j, 0), p);

              // Correct the interpolated value. I.e. if it is close to 1 or
              // 0, make it exactly 1 or 0. Unfortunately, this is required to
              // avoid problems with higher order elements.
              if (std::fabs (matrix_entry - 1.0) < eps)
                matrix_entry = 1.0;
              if (std::fabs (matrix_entry) < eps)
                matrix_entry = 0.0;

              interpolation_matrix(i,j) = matrix_entry;
            }
        }

      // make sure that the row sum of each of the matrices is 1 at this
      // point. this must be so since the shape functions sum up to 1
      for (unsigned int j=0; j<source_fe->dofs_per_face; ++j)
        {
          double sum = 0.;

          for (unsigned int i=0; i<this->dofs_per_face; ++i)
            sum += interpolation_matrix(j,i);

          Assert (std::fabs(sum-1) < eps, ExcInternalError());
        }
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&x_source_fe) != 0)
    {
      // nothing to do here, the FE_Nothing has no degrees of freedom anyway
    }
  else
    AssertThrow (false,(typename FiniteElement<dim,spacedim>::
                        ExcInterpolationNotImplemented()));
}



template <class PolynomialType, int dim, int spacedim>
bool
FE_Q_Base<PolynomialType,dim,spacedim>::hp_constraints_are_implemented () const
{
  return true;
}




template <class PolynomialType, int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Q_Base<PolynomialType,dim,spacedim>::
hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
  // we can presently only compute these identities if both FEs are FE_Qs or
  // if the other one is an FE_Nothing. in the first case, there should be
  // exactly one single DoF of each FE at a vertex, and they should have
  // identical value
  if (dynamic_cast<const FE_Q_Base<PolynomialType,dim,spacedim>*>(&fe_other) != 0)
    {
      return
        std::vector<std::pair<unsigned int, unsigned int> >
        (1, std::make_pair (0U, 0U));
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
      // the FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else if (fe_other.dofs_per_face == 0)
    {
      // if the other element has no elements on faces at all,
      // then it would be impossible to enforce any kind of
      // continuity even if we knew exactly what kind of element
      // we have -- simply because the other element declares
      // that it is discontinuous because it has no DoFs on
      // its faces. in that case, just state that we have no
      // constraints to declare
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <class PolynomialType, int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Q_Base<PolynomialType,dim,spacedim>::
hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
  // we can presently only compute these identities if both FEs are FE_Qs or
  // if the other one is an FE_Nothing
  if (const FE_Q_Base<PolynomialType,dim,spacedim> *fe_q_other
      = dynamic_cast<const FE_Q_Base<PolynomialType,dim,spacedim>*>(&fe_other))
    {
      // dofs are located along lines, so two dofs are identical if they are
      // located at identical positions. if we had only equidistant points, we
      // could simply check for similarity like (i+1)*q == (j+1)*p, but we
      // might have other support points (e.g. Gauss-Lobatto
      // points). Therefore, read the points in unit_support_points for the
      // first coordinate direction. We take the lexicographic ordering of the
      // points in the first direction (i.e., x-direction), which we access
      // between index 1 and p-1 (index 0 and p are vertex dofs).
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
      // the FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else if (fe_other.dofs_per_face == 0)
    {
      // if the other element has no elements on faces at all,
      // then it would be impossible to enforce any kind of
      // continuity even if we knew exactly what kind of element
      // we have -- simply because the other element declares
      // that it is discontinuous because it has no DoFs on
      // its faces. in that case, just state that we have no
      // constraints to declare
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <class PolynomialType, int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Q_Base<PolynomialType,dim,spacedim>::
hp_quad_dof_identities (const FiniteElement<dim,spacedim>        &fe_other) const
{
  // we can presently only compute these identities if both FEs are FE_Qs or
  // if the other one is an FE_Nothing
  if (const FE_Q_Base<PolynomialType,dim,spacedim> *fe_q_other
      = dynamic_cast<const FE_Q_Base<PolynomialType,dim,spacedim>*>(&fe_other))
    {
      // this works exactly like the line case above, except that now we have
      // to have two indices i1, i2 and j1, j2 to characterize the dofs on the
      // face of each of the finite elements. since they are ordered
      // lexicographically along the first line and we have a tensor product,
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
      // the FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else if (fe_other.dofs_per_face == 0)
    {
      // if the other element has no elements on faces at all,
      // then it would be impossible to enforce any kind of
      // continuity even if we knew exactly what kind of element
      // we have -- simply because the other element declares
      // that it is discontinuous because it has no DoFs on
      // its faces. in that case, just state that we have no
      // constraints to declare
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}



template <class PolynomialType, int dim, int spacedim>
FiniteElementDomination::Domination
FE_Q_Base<PolynomialType,dim,spacedim>::
compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_Q_Base<PolynomialType,dim,spacedim> *fe_q_other
      = dynamic_cast<const FE_Q_Base<PolynomialType,dim,spacedim>*>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim> *fe_nothing = dynamic_cast<const FE_Nothing<dim>*>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        {
          return FiniteElementDomination::other_element_dominates;
        }
      else
        {
          // the FE_Nothing has no degrees of freedom and it is typically used in
          // a context where we don't require any continuity along the interface
          return FiniteElementDomination::no_requirements;
        }
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------



template <class PolynomialType, int dim, int spacedim>
void FE_Q_Base<PolynomialType,dim,spacedim>::initialize_unit_support_points
(const std::vector<Point<1> > &points)
{
  const std::vector<unsigned int> &index_map_inverse=
    this->poly_space.get_numbering_inverse();

  Quadrature<1> support_1d(points);
  Quadrature<dim> support_quadrature(support_1d);
  this->unit_support_points.resize(support_quadrature.size());

  for (unsigned int k=0; k<support_quadrature.size(); k++)
    this->unit_support_points[index_map_inverse[k]] = support_quadrature.point(k);
}



template <class PolynomialType, int dim, int spacedim>
void FE_Q_Base<PolynomialType,dim,spacedim>::initialize_unit_face_support_points
(const std::vector<Point<1> > &points)
{
  // no faces in 1d, so nothing to do
  if (dim == 1)
    return;

  const unsigned int codim = dim-1;
  this->unit_face_support_points.resize(Utilities::fixed_power<codim>(q_degree+1));

  // find renumbering of faces and assign from values of quadrature
  std::vector<unsigned int> face_index_map =
    FE_Q_Helper::face_lexicographic_to_hierarchic_numbering<dim>(q_degree);
  Quadrature<1> support_1d(points);
  Quadrature<codim> support_quadrature(support_1d);
  this->unit_face_support_points.resize(support_quadrature.size());

  for (unsigned int k=0; k<support_quadrature.size(); k++)
    this->unit_face_support_points[face_index_map[k]] = support_quadrature.point(k);
}



template <class PolynomialType, int dim, int spacedim>
void
FE_Q_Base<PolynomialType,dim,spacedim>::initialize_quad_dof_index_permutation ()
{
  // for 1D and 2D, do nothing
  if (dim < 3)
    return;

  Assert (this->adjust_quad_dof_index_for_face_orientation_table.n_elements()==8*this->dofs_per_quad,
          ExcInternalError());

  const unsigned int n=q_degree-1;
  Assert(n*n==this->dofs_per_quad, ExcInternalError());

  // alias for the table to fill
  Table<2,int> &data=this->adjust_quad_dof_index_for_face_orientation_table;

  // the dofs on a face are connected to a n x n matrix. for example, for
  // degree==4 we have the following dofs on a quad

  //  ___________
  // |           |
  // |  6  7  8  |
  // |           |
  // |  3  4  5  |
  // |           |
  // |  0  1  2  |
  // |___________|
  //
  // we have dof_no=i+n*j with index i in x-direction and index j in
  // y-direction running from 0 to n-1.  to extract i and j we can use
  // i=dof_no%n and j=dof_no/n. i and j can then be used to construct the
  // rotated and mirrored numbers.


  for (unsigned int local=0; local<this->dofs_per_quad; ++local)
    // face support points are in lexicographic ordering with x running
    // fastest. invert that (y running fastest)
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

  // additionally initialize reordering of line dofs
  for (unsigned int i=0; i<this->dofs_per_line; ++i)
    this->adjust_line_dof_index_for_line_orientation_table[i]=this->dofs_per_line-1-i - i;
}



template <class PolynomialType, int dim, int spacedim>
unsigned int
FE_Q_Base<PolynomialType,dim,spacedim>::
face_to_cell_index (const unsigned int face_index,
                    const unsigned int face,
                    const bool face_orientation,
                    const bool face_flip,
                    const bool face_rotation) const
{
  Assert (face_index < this->dofs_per_face,
          ExcIndexRange(face_index, 0, this->dofs_per_face));
  Assert (face < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange(face, 0, GeometryInfo<dim>::faces_per_cell));

//TODO: we could presumably solve the 3d case below using the
// adjust_quad_dof_index_for_face_orientation_table field. for the
// 2d case, we can't use adjust_line_dof_index_for_line_orientation_table
// since that array is empty (presumably because we thought that
// there are no flipped edges in 2d, but these can happen in
// DoFTools::make_periodicity_constraints, for example). so we
// would need to either fill this field, or rely on derived classes
// implementing this function, as we currently do

  // we need to distinguish between DoFs on vertices, lines and in 3d quads.
  // do so in a sequence of if-else statements
  if (face_index < this->first_face_line_index)
    // DoF is on a vertex
    {
      // get the number of the vertex on the face that corresponds to this DoF,
      // along with the number of the DoF on this vertex
      const unsigned int face_vertex         = face_index / this->dofs_per_vertex;
      const unsigned int dof_index_on_vertex = face_index % this->dofs_per_vertex;

      // then get the number of this vertex on the cell and translate
      // this to a DoF number on the cell
      return (GeometryInfo<dim>::face_to_cell_vertices(face, face_vertex,
                                                       face_orientation,
                                                       face_flip,
                                                       face_rotation)
              * this->dofs_per_vertex
              +
              dof_index_on_vertex);
    }
  else if (face_index < this->first_face_quad_index)
    // DoF is on a face
    {
      // do the same kind of translation as before. we need to only consider
      // DoFs on the lines, i.e., ignoring those on the vertices
      const unsigned int index = face_index - this->first_face_line_index;

      const unsigned int face_line         = index / this->dofs_per_line;
      const unsigned int dof_index_on_line = index % this->dofs_per_line;

      // we now also need to adjust the line index for the case of
      // face orientation, face flips, etc
      unsigned int adjusted_dof_index_on_line;
      switch (dim)
        {
        case 1:
          Assert (false, ExcInternalError());

        case 2:
          // in 2d, only face_flip has a meaning. if it is set, consider
          // dofs in reverse order
          if (face_flip == false)
            adjusted_dof_index_on_line = dof_index_on_line;
          else
            adjusted_dof_index_on_line = this->dofs_per_line - dof_index_on_line - 1;
          break;

        case 3:
          // in 3d, things are difficult. someone will have to think
          // about how this code here should look like, by drawing a bunch
          // of pictures of how all the faces can look like with the various
          // flips and rotations.
          //
          // that said, the Q2 case is easy enough to implement, as is the case
          // where everything is in standard orientation
          Assert ((this->dofs_per_line <= 1) ||
                  ((face_orientation == true) &&
                   (face_flip == false) &&
                   (face_rotation == false)),
                  ExcNotImplemented());
          adjusted_dof_index_on_line = dof_index_on_line;
          break;
        }

      return (this->first_line_index
              + GeometryInfo<dim>::face_to_cell_lines(face, face_line,
                                                      face_orientation,
                                                      face_flip,
                                                      face_rotation)
              * this->dofs_per_line
              +
              adjusted_dof_index_on_line);
    }
  else
    // DoF is on a quad
    {
      Assert (dim >= 3, ExcInternalError());

      // ignore vertex and line dofs
      const unsigned int index = face_index - this->first_face_quad_index;

      // the same is true here as above for the 3d case -- someone will
      // just have to draw a bunch of pictures. in the meantime,
      // we can implement the Q2 case in which it is simple
      Assert ((this->dofs_per_quad <= 1) ||
              ((face_orientation == true) &&
               (face_flip == false) &&
               (face_rotation == false)),
              ExcNotImplemented());
      return (this->first_quad_index
              + face * this->dofs_per_quad
              + index);
    }
}




template <class PolynomialType, int dim, int spacedim>
std::vector<unsigned int>
FE_Q_Base<PolynomialType,dim,spacedim>::get_dpo_vector(const unsigned int deg)
{
  AssertThrow(deg>0,ExcMessage("FE_Q needs to be of degree > 0."));
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}



template <class PolynomialType, int dim, int spacedim>
void
FE_Q_Base<PolynomialType,dim,spacedim>::initialize_constraints
(const std::vector<Point<1> > &points)
{
  Implementation::initialize_constraints (points, *this);
}



template <class PolynomialType, int dim, int spacedim>
const FullMatrix<double> &
FE_Q_Base<PolynomialType,dim,spacedim>
::get_prolongation_matrix (const unsigned int child,
                           const RefinementCase<dim> &refinement_case) const
{
  Assert (refinement_case<RefinementCase<dim>::isotropic_refinement+1,
          ExcIndexRange(refinement_case,0,RefinementCase<dim>::isotropic_refinement+1));
  Assert (refinement_case!=RefinementCase<dim>::no_refinement,
          ExcMessage("Prolongation matrices are only available for refined cells!"));
  Assert (child<GeometryInfo<dim>::n_children(refinement_case),
          ExcIndexRange(child,0,GeometryInfo<dim>::n_children(refinement_case)));

  // initialization upon first request
  if (this->prolongation[refinement_case-1][child].n() == 0)
    {
      Threads::Mutex::ScopedLock lock(this->mutex);

      // if matrix got updated while waiting for the lock
      if (this->prolongation[refinement_case-1][child].n() ==
          this->dofs_per_cell)
        return this->prolongation[refinement_case-1][child];

      // distinguish q/q_dg0 case: only treat Q dofs first
      const unsigned int q_dofs_per_cell = Utilities::fixed_power<dim>(q_degree+1);

      // compute the interpolation matrices in much the same way as we do for
      // the constraints. it's actually simpler here, since we don't have this
      // weird renumbering stuff going on. The trick is again that we the
      // interpolation matrix is formed by a permutation of the indices of the
      // cell matrix. The value eps is used a threshold to decide when certain
      // evaluations of the Lagrange polynomials are zero or one.
      const double eps = 1e-15*q_degree*dim;

#ifdef DEBUG
      // in DEBUG mode, check that the evaluation of support points in the
      // current numbering gives the identity operation
      for (unsigned int i=0; i<q_dofs_per_cell; ++i)
        {
          Assert (std::fabs (1.-this->poly_space.compute_value
                             (i, this->unit_support_points[i])) < eps,
                  ExcInternalError("The Lagrange polynomial does not evaluate "
                                   "to one or zero in a nodal point. "
                                   "This typically indicates that the "
                                   "polynomial interpolation is "
                                   "ill-conditioned such that round-off "
                                   "prevents the sum to be one."));
          for (unsigned int j=0; j<q_dofs_per_cell; ++j)
            if (j!=i)
              Assert (std::fabs (this->poly_space.compute_value
                                 (i, this->unit_support_points[j])) < eps,
                      ExcInternalError("The Lagrange polynomial does not evaluate "
                                       "to one or zero in a nodal point. "
                                       "This typically indicates that the "
                                       "polynomial interpolation is "
                                       "ill-conditioned such that round-off "
                                       "prevents the sum to be one."));
        }
#endif

      // to efficiently evaluate the polynomial at the subcell, make use of
      // the tensor product structure of this element and only evaluate 1D
      // information from the polynomial. This makes the cost of this function
      // almost negligible also for high order elements
      const unsigned int dofs1d = q_degree+1;
      std::vector<Table<2,double> >
      subcell_evaluations (dim, Table<2,double>(dofs1d, dofs1d));
      const std::vector<unsigned int> &index_map_inverse =
        this->poly_space.get_numbering_inverse();

      // helper value: step size how to walk through diagonal and how many
      // points we have left apart from the first dimension
      unsigned int step_size_diag = 0;
      {
        unsigned int factor = 1;
        for (unsigned int d=0; d<dim; ++d)
          {
            step_size_diag += factor;
            factor *= dofs1d;
          }
      }

      FullMatrix<double> prolongate (this->dofs_per_cell, this->dofs_per_cell);

      // go through the points in diagonal to capture variation in all
      // directions simultaneously
      for (unsigned int j=0; j<dofs1d; ++j)
        {
          const unsigned int diag_comp = index_map_inverse[j*step_size_diag];
          const Point<dim> p_subcell = this->unit_support_points[diag_comp];
          const Point<dim> p_cell =
            GeometryInfo<dim>::child_to_cell_coordinates (p_subcell, child,
                                                          refinement_case);
          for (unsigned int i=0; i<dofs1d; ++i)
            for (unsigned int d=0; d<dim; ++d)
              {
                // evaluate along line where only x is different from zero
                Point<dim> point;
                point[0] = p_cell[d];
                const double cell_value =
                  this->poly_space.compute_value(index_map_inverse[i], point);

                // cut off values that are too small. note that we have here
                // Lagrange interpolation functions, so they should be zero at
                // almost all points, and one at the others, at least on the
                // subcells. so set them to their exact values
                //
                // the actual cut-off value is somewhat fuzzy, but it works
                // for 2e-13*degree*dim (see above), which is kind of
                // reasonable given that we compute the values of the
                // polynomials via an degree-step recursion and then multiply
                // the 1d-values. this gives us a linear growth in degree*dim,
                // times a small constant.
                //
                // the embedding matrix is given by applying the inverse of
                // the subcell matrix on the cell_interpolation matrix. since
                // the subcell matrix is actually only a permutation vector,
                // all we need to do is to switch the rows we write the data
                // into. moreover, cut off very small values here
                if (std::fabs(cell_value) < eps)
                  subcell_evaluations[d](j,i) = 0;
                else
                  subcell_evaluations[d](j,i) = cell_value;
              }
        }

      // now expand from 1D info. block innermost dimension (x_0) in order to
      // avoid difficult checks at innermost loop
      unsigned int j_indices[dim];
      FE_Q_Helper::zero_indices<dim> (j_indices);
      for (unsigned int j=0; j<q_dofs_per_cell; j+=dofs1d)
        {
          unsigned int i_indices[dim];
          FE_Q_Helper::zero_indices<dim> (i_indices);
          for (unsigned int i=0; i<q_dofs_per_cell; i+=dofs1d)
            {
              double val_extra_dim = 1.;
              for (unsigned int d=1; d<dim; ++d)
                val_extra_dim *= subcell_evaluations[d](j_indices[d-1],
                                                        i_indices[d-1]);

              // innermost sum where we actually compute. the same as
              // prolongate(j,i) = this->poly_space.compute_value (i, p_cell)
              for (unsigned int jj=0; jj<dofs1d; ++jj)
                {
                  const unsigned int j_ind = index_map_inverse[j+jj];
                  for (unsigned int ii=0; ii<dofs1d; ++ii)
                    prolongate(j_ind,index_map_inverse[i+ii])
                      = val_extra_dim * subcell_evaluations[0](jj,ii);
                }

              // update indices that denote the tensor product position. a bit
              // fuzzy and therefore not done for innermost x_0 direction
              FE_Q_Helper::increment_indices<dim> (i_indices, dofs1d);
            }
          Assert (i_indices[dim-1] == 1, ExcInternalError());
          FE_Q_Helper::increment_indices<dim> (j_indices, dofs1d);
        }

      // the discontinuous node is simply mapped on the discontinuous node on
      // the child element
      if (q_dofs_per_cell < this->dofs_per_cell)
        prolongate(q_dofs_per_cell,q_dofs_per_cell) = 1.;

      // and make sure that the row sum is 1. this must be so since for this
      // element, the shape functions add up to one
#ifdef DEBUG
      for (unsigned int row=0; row<this->dofs_per_cell; ++row)
        {
          double sum = 0;
          for (unsigned int col=0; col<this->dofs_per_cell; ++col)
            sum += prolongate(row,col);
          Assert (std::fabs(sum-1.) <
                  std::max(eps, 5e-16*std::sqrt(this->dofs_per_cell)),
                  ExcInternalError("The entries in a row of the local "
                                   "prolongation matrix do not add to one. "
                                   "This typically indicates that the "
                                   "polynomial interpolation is "
                                   "ill-conditioned such that round-off "
                                   "prevents the sum to be one."));
        }
#endif

      // swap matrices
      prolongate.swap(const_cast<FullMatrix<double> &>
                      (this->prolongation[refinement_case-1][child]));
    }

  // finally return the matrix
  return this->prolongation[refinement_case-1][child];
}



template <class PolynomialType, int dim, int spacedim>
const FullMatrix<double> &
FE_Q_Base<PolynomialType,dim,spacedim>
::get_restriction_matrix (const unsigned int child,
                          const RefinementCase<dim> &refinement_case) const
{
  Assert (refinement_case<RefinementCase<dim>::isotropic_refinement+1,
          ExcIndexRange(refinement_case,0,RefinementCase<dim>::isotropic_refinement+1));
  Assert (refinement_case!=RefinementCase<dim>::no_refinement,
          ExcMessage("Restriction matrices are only available for refined cells!"));
  Assert (child<GeometryInfo<dim>::n_children(refinement_case),
          ExcIndexRange(child,0,GeometryInfo<dim>::n_children(refinement_case)));

  // initialization upon first request
  if (this->restriction[refinement_case-1][child].n() == 0)
    {
      Threads::Mutex::ScopedLock lock(this->mutex);

      // if matrix got updated while waiting for the lock...
      if (this->restriction[refinement_case-1][child].n() ==
          this->dofs_per_cell)
        return this->restriction[refinement_case-1][child];

      FullMatrix<double> restriction(this->dofs_per_cell, this->dofs_per_cell);
      // distinguish q/q_dg0 case
      const unsigned int q_dofs_per_cell = Utilities::fixed_power<dim>(q_degree+1);

      // for Lagrange interpolation polynomials based on equidistant points,
      // construction of the restriction matrices is relatively simple. the
      // reason is that in this case the interpolation points on the mother
      // cell are always also interpolation points for some shape function on
      // one or the other child.
      //
      // in the general case with non-equidistant points, we need to actually
      // do an interpolation. thus, we take the interpolation points on the
      // mother cell and evaluate the shape functions of the child cell on
      // those points. it does not hurt in the equidistant case because then
      // simple one shape function evaluates to one and the others to zero.
      //
      // this element is non-additive in all its degrees of freedom by
      // default, which requires care in downstream use. fortunately, even the
      // interpolation on non-equidistant points is invariant under the
      // assumption that whenever a row makes a non-zero contribution to the
      // mother's residual, the correct value is interpolated.

      const double eps = 1e-15*q_degree*dim;
      const std::vector<unsigned int> &index_map_inverse =
        this->poly_space.get_numbering_inverse();

      const unsigned int dofs1d = q_degree+1;
      std::vector<Tensor<1,dim> > evaluations1d (dofs1d);

      restriction.reinit(this->dofs_per_cell, this->dofs_per_cell);

      for (unsigned int i=0; i<q_dofs_per_cell; ++i)
        {
          unsigned int mother_dof = index_map_inverse[i];
          const Point<dim> p_cell = this->unit_support_points[mother_dof];

          // check whether this interpolation point is inside this child cell
          const Point<dim> p_subcell
            = GeometryInfo<dim>::cell_to_child_coordinates (p_cell, child,
                                                            refinement_case);
          if (GeometryInfo<dim>::is_inside_unit_cell (p_subcell))
            {
              // same logic as in initialize_embedding to evaluate the
              // polynomial faster than from the tensor product: since we
              // evaluate all polynomials, it is much faster to just compute
              // the 1D values for all polynomials before and then get the
              // dim-data.
              for (unsigned int j=0; j<dofs1d; ++j)
                for (unsigned int d=0; d<dim; ++d)
                  {
                    Point<dim> point;
                    point[0] = p_subcell[d];
                    evaluations1d[j][d] =
                      this->poly_space.compute_value(index_map_inverse[j], point);
                  }
              unsigned int j_indices[dim];
              FE_Q_Helper::zero_indices<dim> (j_indices);
              double sum_check = 0;
              for (unsigned int j = 0; j<q_dofs_per_cell; j += dofs1d)
                {
                  double val_extra_dim = 1.;
                  for (unsigned int d=1; d<dim; ++d)
                    val_extra_dim *= evaluations1d[j_indices[d-1]][d];
                  for (unsigned int jj=0; jj<dofs1d; ++jj)
                    {

                      // find the child shape function(s) corresponding to
                      // this point. Usually this is just one function;
                      // however, when we use FE_Q on arbitrary nodes a parent
                      // support point might not be a child support point, and
                      // then we will get more than one nonzero value per
                      // row. Still, the values should sum up to 1
                      const double val
                        = val_extra_dim * evaluations1d[jj][0];
                      const unsigned int child_dof =
                        index_map_inverse[j+jj];
                      if (std::fabs (val-1.) < eps)
                        restriction(mother_dof,child_dof)=1.;
                      else if (std::fabs(val) > eps)
                        restriction(mother_dof,child_dof)=val;
                      sum_check += val;
                    }
                  FE_Q_Helper::increment_indices<dim> (j_indices, dofs1d);
                }
              Assert (std::fabs(sum_check-1) <
                      std::max(eps, 5e-16*std::sqrt(this->dofs_per_cell)),
                      ExcInternalError("The entries in a row of the local "
                                       "restriction matrix do not add to one. "
                                       "This typically indicates that the "
                                       "polynomial interpolation is "
                                       "ill-conditioned such that round-off "
                                       "prevents the sum to be one."));
            }

          // part for FE_Q_DG0
          if (q_dofs_per_cell < this->dofs_per_cell)
            restriction(this->dofs_per_cell-1,this->dofs_per_cell-1) =
              1./GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case));
        }

      // swap matrices
      restriction.swap(const_cast<FullMatrix<double> &>
                       (this->restriction[refinement_case-1][child]));
    }

  return this->restriction[refinement_case-1][child];
}



//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------


template <class PolynomialType, int dim, int spacedim>
bool
FE_Q_Base<PolynomialType,dim,spacedim>::has_support_on_face
(const unsigned int shape_index,
 const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
          ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

  // in 1d, things are simple. since there is only one degree of freedom per
  // vertex in this class, the first is on vertex 0 (==face 0 in some sense),
  // the second on face 1:
  if (dim == 1)
    return (((shape_index == 0) && (face_index == 0)) ||
            ((shape_index == 1) && (face_index == 1)));

  // first, special-case interior shape functions, since they have no support
  // no-where on the boundary
  if (((dim==2) && (shape_index>=this->first_quad_index))
      ||
      ((dim==3) && (shape_index>=this->first_hex_index)))
    return false;

  // let's see whether this is a vertex
  if (shape_index < this->first_line_index)
    {
      // for Q elements, there is one dof per vertex, so
      // shape_index==vertex_number. check whether this vertex is on the given
      // face. thus, for each face, give a list of vertices
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

      // in 2d, the line is the face, so get the line index
      if (dim == 2)
        return (line_index == face_index);
      else if (dim == 3)
        {
          // silence compiler warning
          const unsigned int lines_per_face =
            dim == 3 ? GeometryInfo<dim>::lines_per_face : 1;
          // see whether the given line is on the given face.
          for (unsigned int l=0; l<lines_per_face; ++l)
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

      // in 2d, cell bubble are zero on all faces. but we have treated this
      // case above already
      Assert (dim != 2, ExcInternalError());

      // in 3d, quad_index=face_index
      if (dim == 3)
        return (quad_index == face_index);
      else
        Assert (false, ExcNotImplemented());
    }
  else
    // dof on hex
    {
      // can only happen in 3d, but this case has already been covered above
      Assert (false, ExcNotImplemented());
      return false;
    }

  // we should not have gotten here
  Assert (false, ExcInternalError());
  return false;
}



template <typename PolynomialType, int dim, int spacedim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FE_Q_Base<PolynomialType,dim,spacedim>::get_constant_modes () const
{
  Table<2,bool> constant_modes(1, this->dofs_per_cell);
  // We here just care for the constant mode due to the polynomial space
  // without any enrichments
  // As there may be more constant modes derived classes may to implement this
  // themselves. An example for this is FE_Q_DG0.
  for (unsigned int i=0; i<Utilities::fixed_power<dim>(q_degree+1); ++i)
    constant_modes(0, i) = true;
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (constant_modes, std::vector<unsigned int>(1, 0));
}


// explicit instantiations
#include "fe_q_base.inst"

DEAL_II_NAMESPACE_CLOSE
