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


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <cmath>
#include <algorithm>


DEAL_II_NAMESPACE_OPEN




template<int dim, int spacedim>
MappingQ1<dim,spacedim>::MappingQ1 ()
  :
  MappingQGeneric<dim,spacedim> (1)
{}



template<int dim, int spacedim>
MappingQ1<dim,spacedim>::MappingQ1 (const unsigned int p)
  :
  MappingQGeneric<dim,spacedim> (p)
{}



namespace internal
{
  namespace MappingQ1
  {
    // These are left as templates on the spatial dimension (even though dim
    // == spacedim must be true for them to make sense) because templates are
    // expanded before the compiler eliminates code due to the 'if (dim ==
    // spacedim)' statement (see the body of the general
    // transform_real_to_unit_cell).
    template<int spacedim>
    Point<1>
    transform_real_to_unit_cell
    (const std_cxx11::array<Point<spacedim>, GeometryInfo<1>::vertices_per_cell> &vertices,
     const Point<spacedim> &p)
    {
      Assert(spacedim == 1, ExcInternalError());
      return Point<1>((p[0] - vertices[0](0))/(vertices[1](0) - vertices[0](0)));
    }



    template<int spacedim>
    Point<2>
    transform_real_to_unit_cell
    (const std_cxx11::array<Point<spacedim>, GeometryInfo<2>::vertices_per_cell> &vertices,
     const Point<spacedim> &p)
    {
      Assert(spacedim == 2, ExcInternalError());
      const double x = p(0);
      const double y = p(1);

      const double x0 = vertices[0](0);
      const double x1 = vertices[1](0);
      const double x2 = vertices[2](0);
      const double x3 = vertices[3](0);

      const double y0 = vertices[0](1);
      const double y1 = vertices[1](1);
      const double y2 = vertices[2](1);
      const double y3 = vertices[3](1);

      const double a = (x1 - x3)*(y0 - y2) - (x0 - x2)*(y1 - y3);
      const double b = -(x0 - x1 - x2 + x3)*y + (x - 2*x1 + x3)*y0 - (x - 2*x0 + x2)*y1
                       - (x - x1)*y2 + (x - x0)*y3;
      const double c = (x0 - x1)*y - (x - x1)*y0 + (x - x0)*y1;

      const double discriminant = b*b - 4*a*c;
      // exit if the point is not in the cell (this is the only case where the
      // discriminant is negative)
      if (discriminant < 0.0)
        {
          AssertThrow (false,
                       (typename Mapping<spacedim,spacedim>::ExcTransformationFailed()));
        }

      double eta1;
      double eta2;
      // special case #1: if a is zero, then use the linear formula
      if (a == 0.0 && b != 0.0)
        {
          eta1 = -c/b;
          eta2 = -c/b;
        }
      // special case #2: if c is very small:
      else if (std::abs(c/b) < 1e-12)
        {
          eta1 = (-b - std::sqrt(discriminant)) / (2*a);
          eta2 = (-b + std::sqrt(discriminant)) / (2*a);
        }
      // finally, use the numerically stable version of the quadratic formula:
      else
        {
          eta1 = 2*c / (-b - std::sqrt(discriminant));
          eta2 = 2*c / (-b + std::sqrt(discriminant));
        }
      // pick the one closer to the center of the cell.
      const double eta = (std::abs(eta1 - 0.5) < std::abs(eta2 - 0.5)) ? eta1 : eta2;

      /*
       * There are two ways to compute xi from eta, but either one may have a
       * zero denominator.
       */
      const double subexpr0 = -eta*x2 + x0*(eta - 1);
      const double xi_denominator0 = eta*x3 - x1*(eta - 1) + subexpr0;
      const double max_x = std::max(std::max(std::abs(x0), std::abs(x1)),
                                    std::max(std::abs(x2), std::abs(x3)));

      if (std::abs(xi_denominator0) > 1e-10*max_x)
        {
          const double xi = (x + subexpr0)/xi_denominator0;
          return Point<2>(xi, eta);
        }
      else
        {
          const double max_y = std::max(std::max(std::abs(y0), std::abs(y1)),
                                        std::max(std::abs(y2), std::abs(y3)));
          const double subexpr1 = -eta*y2 + y0*(eta - 1);
          const double xi_denominator1 = eta*y3 - y1*(eta - 1) + subexpr1;
          if (std::abs(xi_denominator1) > 1e-10*max_y)
            {
              const double xi = (subexpr1 + y)/xi_denominator1;
              return Point<2>(xi, eta);
            }
          else // give up and try Newton iteration
            {
              AssertThrow (false,
                           (typename Mapping<spacedim,spacedim>::ExcTransformationFailed()));
            }
        }
      // bogus return to placate compiler. It should not be possible to get
      // here.
      Assert(false, ExcInternalError());
      return Point<2>(std::numeric_limits<double>::quiet_NaN(),
                      std::numeric_limits<double>::quiet_NaN());
    }



    template<int spacedim>
    Point<3>
    transform_real_to_unit_cell
    (const std_cxx11::array<Point<spacedim>, GeometryInfo<3>::vertices_per_cell> &/*vertices*/,
     const Point<spacedim> &/*p*/)
    {
      // It should not be possible to get here
      Assert(false, ExcInternalError());
      return Point<3>();
    }
  }
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  std::vector<Point<spacedim> > &a) const
{
  std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
  vertices = this->get_vertices(cell);

  a.resize(GeometryInfo<dim>::vertices_per_cell);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a[i] = vertices[i];
}



/* For an explanation of the  KA and Kb
   arrays see the comments in the declaration of
   transform_real_to_unit_cell_initial_guess */
namespace
{
  template <int dim>
  struct TransformR2UInitialGuess
  {
    static const double KA[GeometryInfo<dim>::vertices_per_cell][dim];
    static const double Kb[GeometryInfo<dim>::vertices_per_cell];
  };


  /*
    Octave code:
    M=[0 1; 1 1];
    K1 = transpose(M) * inverse (M*transpose(M));
    printf ("{%f, %f},\n", K1' );
  */
  template <>
  const double
  TransformR2UInitialGuess<1>::
  KA[GeometryInfo<1>::vertices_per_cell][1] =
  {
    {-1.000000},
    {1.000000}
  };

  template <>
  const double
  TransformR2UInitialGuess<1>::
  Kb[GeometryInfo<1>::vertices_per_cell] = {1.000000, 0.000000};


  /*
    Octave code:
    M=[0 1 0 1;0 0 1 1;1 1 1 1];
    K2 = transpose(M) * inverse (M*transpose(M));
    printf ("{%f, %f, %f},\n", K2' );
  */
  template <>
  const double
  TransformR2UInitialGuess<2>::
  KA[GeometryInfo<2>::vertices_per_cell][2] =
  {
    {-0.500000, -0.500000},
    { 0.500000, -0.500000},
    {-0.500000,  0.500000},
    { 0.500000,  0.500000}
  };

  /*
    Octave code:
    M=[0 1 0 1 0 1 0 1;0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1; 1 1 1 1 1 1 1 1];
    K3 = transpose(M) * inverse (M*transpose(M))
    printf ("{%f, %f, %f, %f},\n", K3' );
  */
  template <>
  const double
  TransformR2UInitialGuess<2>::
  Kb[GeometryInfo<2>::vertices_per_cell] =
  {0.750000,0.250000,0.250000,-0.250000 };


  template <>
  const double
  TransformR2UInitialGuess<3>::
  KA[GeometryInfo<3>::vertices_per_cell][3] =
  {
    {-0.250000, -0.250000, -0.250000},
    { 0.250000, -0.250000, -0.250000},
    {-0.250000,  0.250000, -0.250000},
    { 0.250000,  0.250000, -0.250000},
    {-0.250000, -0.250000,  0.250000},
    { 0.250000, -0.250000,  0.250000},
    {-0.250000,  0.250000,  0.250000},
    { 0.250000,  0.250000,  0.250000}

  };


  template <>
  const double
  TransformR2UInitialGuess<3>::
  Kb[GeometryInfo<3>::vertices_per_cell] =
  {0.500000,0.250000,0.250000,0.000000,0.250000,0.000000,0.000000,-0.250000};

}


template<int dim, int spacedim>
Point<dim>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell_initial_guess (const std::vector<Point<spacedim> > &vertex,
                                           const Point<spacedim>               &p) const
{
  Point<dim> p_unit;

  FullMatrix<double>  KA(GeometryInfo<dim>::vertices_per_cell, dim);
  Vector <double>  Kb(GeometryInfo<dim>::vertices_per_cell);

  KA.fill( (double *)(TransformR2UInitialGuess<dim>::KA) );
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    Kb(i)=(TransformR2UInitialGuess<dim>::Kb)[i];

  FullMatrix<double> Y(spacedim, GeometryInfo<dim>::vertices_per_cell);
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
    for (unsigned int i=0; i<spacedim; ++i)
      Y(i,v) = vertex[v][i];

  FullMatrix<double> A(spacedim,dim);
  Y.mmult(A,KA); // A = Y*KA
  Vector< double > b(spacedim);
  Y.vmult(b,Kb); // b = Y*Kb

  for (unsigned int i=0; i<spacedim; ++i)
    b(i) -= p[i];
  b*=-1;

  Vector< double > dest(dim);

  FullMatrix<double> A_1(dim,spacedim);
  if (dim<spacedim)
    A_1.left_invert(A);
  else
    A_1.invert(A);

  A_1.vmult(dest,b); //A^{-1}*b

  for (unsigned int i=0; i<dim; ++i)
    p_unit[i]=dest(i);

  return p_unit;
}



template<int dim, int spacedim>
Point<dim>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<spacedim>                            &p) const
{
  // Use the exact formula if available
  if (dim == spacedim && (dim == 1 || dim == 2))
    {
      // The dimension-dependent algorithms are much faster (about 25-45x in
      // 2D) but fail most of the time when the given point (p) is not in the
      // cell. The dimension-independent Newton algorithm given below is
      // slower, but more robust (though it still sometimes fails). Therefore
      // this function implements the following strategy based on the
      // p's dimension:
      //
      // * In 1D this mapping is linear, so the mapping is always invertible
      //   (and the exact formula is known) as long as the cell has non-zero
      //   length.
      // * In 2D the exact (quadratic) formula is called first. If either the
      //   exact formula does not succeed (negative discriminant in the
      //   quadratic formula) or succeeds but finds a solution outside of the
      //   unit cell, then the Newton solver is called. The rationale for the
      //   second choice is that the exact formula may provide two different
      //   answers when mapping a point outside of the real cell, but the
      //   Newton solver (if it converges) will only return one answer.
      //   Otherwise the exact formula successfully found a point in the unit
      //   cell and that value is returned.
      // * In 3D there is no (known to the authors) exact formula, so the Newton
      //   algorithm is used.
      const std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
      vertices = this->get_vertices(cell);
      try
        {
          Point<dim> point = internal::MappingQ1::transform_real_to_unit_cell(vertices, p);

          if (dim == 1)
            {
              // formula not subject to any issues
              return point;
            }
          else if (dim == 2)
            {
              // formula not guaranteed to work for points outside of the cell
              const double eps = 1e-15;
              if (-eps <= point(1) && point(1) <= 1 + eps &&
                  -eps <= point(0) && point(0) <= 1 + eps)
                {
                  return point;
                }
            }
          else
            {
              Assert(false, ExcInternalError());
            }
        }
      catch (const typename Mapping<spacedim,spacedim>::ExcTransformationFailed &)
        {
          // continue on to the standard Newton code
        }
    }
  // Find the initial value for the
  // Newton iteration by a normal projection
  // to the least square plane determined by
  // the vertices of the cell
  std::vector<Point<spacedim> > a;
  compute_mapping_support_points (cell,a);
  Point<dim> initial_p_unit =
    transform_real_to_unit_cell_initial_guess(a,p);

  // if dim==1 there is nothing
  // else to do to the initial
  // value, and it is the answer
  if (dim == 1)
    return initial_p_unit;
  else
    {
      // use the full mapping. in case the
      // function above should have given us
      // something back that lies outside the
      // unit cell (that might happen because
      // either the function computing an
      // initial guess gave us a poor initial
      // guess or for the following reason:
      // we call this function here in the Q1
      // mapping to produce an initial guess
      // for a higher order mapping, but
      // we may have given a point 'p' that
      // lies inside the cell with the higher
      // order mapping, but outside the
      // Q1-mapped reference cell), then
      // project it back into the reference
      // cell in hopes that this gives a
      // better starting point to the
      // following iteration
//TODO: the following line was added in r25581 but it leads to
// changes in the test results. investigate why this is so --
// it shouldn't really make any difference...
//      initial_p_unit = GeometryInfo<dim>::project_to_unit_cell(initial_p_unit);

      const Quadrature<dim> point_quadrature(initial_p_unit);

      UpdateFlags update_flags = update_quadrature_points | update_jacobians;
      if (spacedim>dim)
        update_flags |= update_jacobian_grads;

      //TODO: Use get_data() here once MappingQ is no longer derived from
      //MappingQ1. this doesn't currently work because we here really need
      //a Q1 InternalData, but MappingQGeneric produces one with the
      //polynomial degree of the MappingQ
      std_cxx11::unique_ptr<InternalData> mdata (new InternalData(1));
      mdata->initialize (this->requires_update_flags (update_flags), point_quadrature, 1);

      compute_mapping_support_points (cell, mdata->mapping_support_points);
      // The support points have to be at
      // least as many as there are
      // vertices.
      Assert(mdata->mapping_support_points.size() >=
             GeometryInfo<dim>::vertices_per_cell,
             ExcInternalError());
      // Ignore non vertex support points.
      mdata->mapping_support_points.resize(GeometryInfo<dim>::vertices_per_cell);

      // perform the Newton iteration and
      // return the result. note that this
      // statement may throw an exception, which
      // we simply pass up to the caller
      return transform_real_to_unit_cell_internal(cell, p, initial_p_unit,
                                                  *mdata);
    }
}



template<int dim, int spacedim>
Point<dim>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell_internal
(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
 const Point<spacedim>                            &p,
 const Point<dim>                                 &initial_p_unit,
 InternalData                                     &mdata) const
{
  const unsigned int n_shapes=mdata.shape_values.size();
  (void)n_shapes;
  Assert(n_shapes!=0, ExcInternalError());
  AssertDimension (mdata.shape_derivatives.size(), n_shapes);

  std::vector<Point<spacedim> > &points=mdata.mapping_support_points;
  AssertDimension (points.size(), n_shapes);


  // Newton iteration to solve
  //    f(x)=p(x)-p=0
  // where we are looking for 'x' and p(x) is the forward transformation
  // from unit to real cell. We solve this using a Newton iteration
  //    x_{n+1}=x_n-[f'(x)]^{-1}f(x)
  // The start value is set to be the linear approximation to the cell

  // The shape values and derivatives of the mapping at this point are
  // previously computed.

  Point<dim> p_unit = initial_p_unit;

  mdata.compute_shape_function_values(std::vector<Point<dim> > (1, p_unit));

  Point<spacedim> p_real = transform_unit_to_real_cell_internal(mdata);
  Tensor<1,spacedim> f = p_real-p;

  // early out if we already have our point
  if (f.norm_square() < 1e-24 * cell->diameter() * cell->diameter())
    return p_unit;

  // we need to compare the position of the computed p(x) against the given
  // point 'p'. We will terminate the iteration and return 'x' if they are
  // less than eps apart. The question is how to choose eps -- or, put maybe
  // more generally: in which norm we want these 'p' and 'p(x)' to be eps
  // apart.
  //
  // the question is difficult since we may have to deal with very elongated
  // cells where we may achieve 1e-12*h for the distance of these two points
  // in the 'long' direction, but achieving this tolerance in the 'short'
  // direction of the cell may not be possible
  //
  // what we do instead is then to terminate iterations if
  //    \| p(x) - p \|_A < eps
  // where the A-norm is somehow induced by the transformation of the cell.
  // in particular, we want to measure distances relative to the sizes of
  // the cell in its principal directions.
  //
  // to define what exactly A should be, note that to first order we have
  // the following (assuming that x* is the solution of the problem, i.e.,
  // p(x*)=p):
  //    p(x) - p = p(x) - p(x*)
  //             = -grad p(x) * (x*-x) + higher order terms
  // This suggest to measure with a norm that corresponds to
  //    A = {[grad p(x]^T [grad p(x)]}^{-1}
  // because then
  //    \| p(x) - p \|_A  \approx  \| x - x* \|
  // Consequently, we will try to enforce that
  //    \| p(x) - p \|_A  =  \| f \|  <=  eps
  //
  // Note that using this norm is a bit dangerous since the norm changes
  // in every iteration (A isn't fixed by depends on xk). However, if the
  // cell is not too deformed (it may be stretched, but not twisted) then
  // the mapping is almost linear and A is indeed constant or nearly so.
  const double eps = 1.e-11;
  const unsigned int newton_iteration_limit = 20;

  unsigned int newton_iteration = 0;
  double last_f_weighted_norm;
  do
    {
#ifdef DEBUG_TRANSFORM_REAL_TO_UNIT_CELL
      std::cout << "Newton iteration " << newton_iteration << std::endl;
#endif

      // f'(x)
      Tensor<2,spacedim> df;
      for (unsigned int k=0; k<mdata.n_shape_functions; ++k)
        {
          const Tensor<1,dim> &grad_transform=mdata.derivative(0,k);
          const Point<spacedim> &point=points[k];

          for (unsigned int i=0; i<spacedim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              df[i][j]+=point[i]*grad_transform[j];
        }

      // Solve  [f'(x)]d=f(x)
      Tensor<1,spacedim> delta;
      Tensor<2,spacedim> df_inverse = invert(df);
      contract (delta, df_inverse, static_cast<const Tensor<1,spacedim>&>(f));

#ifdef DEBUG_TRANSFORM_REAL_TO_UNIT_CELL
      std::cout << "   delta=" << delta  << std::endl;
#endif

      // do a line search
      double step_length = 1;
      do
        {
          // update of p_unit. The spacedim-th component of transformed point
          // is simply ignored in codimension one case. When this component is
          // not zero, then we are projecting the point to the surface or
          // curve identified by the cell.
          Point<dim> p_unit_trial = p_unit;
          for (unsigned int i=0; i<dim; ++i)
            p_unit_trial[i] -= step_length * delta[i];

          // shape values and derivatives
          // at new p_unit point
          mdata.compute_shape_function_values(std::vector<Point<dim> > (1, p_unit_trial));

          // f(x)
          Point<spacedim> p_real_trial = transform_unit_to_real_cell_internal(mdata);
          const Tensor<1,spacedim> f_trial = p_real_trial-p;

#ifdef DEBUG_TRANSFORM_REAL_TO_UNIT_CELL
          std::cout << "     step_length=" << step_length << std::endl
                    << "       ||f ||   =" << f.norm() << std::endl
                    << "       ||f*||   =" << f_trial.norm() << std::endl
                    << "       ||f*||_A =" << (df_inverse * f_trial).norm() << std::endl;
#endif

          // see if we are making progress with the current step length
          // and if not, reduce it by a factor of two and try again
          //
          // strictly speaking, we should probably use the same norm as we use
          // for the outer algorithm. in practice, line search is just a
          // crutch to find a "reasonable" step length, and so using the l2
          // norm is probably just fine
          if (f_trial.norm() < f.norm())
            {
              p_real = p_real_trial;
              p_unit = p_unit_trial;
              f = f_trial;
              break;
            }
          else if (step_length > 0.05)
            step_length /= 2;
          else
            AssertThrow (false,
                         (typename Mapping<dim,spacedim>::ExcTransformationFailed()));
        }
      while (true);

      ++newton_iteration;
      if (newton_iteration > newton_iteration_limit)
        AssertThrow (false,
                     (typename Mapping<dim,spacedim>::ExcTransformationFailed()));
      last_f_weighted_norm = (df_inverse * f).norm();
    }
  while (last_f_weighted_norm > eps);

  return p_unit;
}



/*
  This function becomes a little tricky in dimension <2,3>.
  There is a surface embedded in R^3 and we pass a point p in R^3, that
  is most likely not lying on the surface.
  We then ask,
  what point in R^2 (hopefully in the unit cell) satisfies that
  map(x) = p.

  An appropriate modification of this question is:
  Find x in R^2 and alpha in R such that

  map(x) + alpha * normal(x) = p


 */

template<>
Point<2>
MappingQ1<2,3>::
transform_real_to_unit_cell_internal (const Triangulation<2,3>::cell_iterator &cell,
                                      const Point<3> &p,
                                      const Point<2> &initial_p_unit,
                                      InternalData   &mdata) const
{
  return
    transform_real_to_unit_cell_internal_codim1(cell, p, initial_p_unit,
                                                mdata);
}




template<>
Point<1>
MappingQ1<1,2>::
transform_real_to_unit_cell_internal (const Triangulation<1,2>::cell_iterator &cell,
                                      const Point<2> &p,
                                      const Point<1> &initial_p_unit,
                                      InternalData   &mdata) const
{
  return
    transform_real_to_unit_cell_internal_codim1(cell, p, initial_p_unit,
                                                mdata);
}


template<>
Point<1>
MappingQ1<1,3>::
transform_real_to_unit_cell_internal (const Triangulation<1,3>::cell_iterator &/*cell*/,
                                      const Point<3> &/*p*/,
                                      const Point<1> &/*initial_p_unit*/,
                                      InternalData   &/*mdata*/) const
{
  Assert(false, ExcNotImplemented());
  return Point<1>();
}




template<int dim, int spacedim>
template<int dim_>
Point<dim_>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell_internal_codim1
(const typename Triangulation<dim_,dim_ + 1>::cell_iterator &cell,
 const Point<dim_ + 1>                                      &p,
 const Point<dim_ >                                         &initial_p_unit,
 typename MappingQ1<dim,spacedim>::InternalData                      &mdata) const
{
  const unsigned int spacedim1 = dim_+1;
  const unsigned int dim1 = dim_;


  const unsigned int n_shapes=mdata.shape_values.size();
  (void)n_shapes;
  Assert(n_shapes!=0, ExcInternalError());
  Assert(mdata.shape_derivatives.size()==n_shapes, ExcInternalError());
  Assert(mdata.shape_second_derivatives.size()==n_shapes, ExcInternalError());

  std::vector<Point<spacedim1> > &points=mdata.mapping_support_points;
  Assert(points.size()==n_shapes, ExcInternalError());

  Point<spacedim1> p_minus_F;

  Tensor<1,spacedim1>  DF[dim1];
  Tensor<1,spacedim1>  D2F[dim1][dim1];

  Point<dim1> p_unit = initial_p_unit;
  Point<dim1> f;
  Tensor<2,dim1>  df;

  //Evaluate first and second derivatives
  mdata.compute_shape_function_values(std::vector<Point<dim> > (1, p_unit));

  for (unsigned int k=0; k<mdata.n_shape_functions; ++k)
    {
      const Tensor<1,dim1>   &grad_phi_k = mdata.derivative(0,k);
      const Tensor<2,dim1>   &hessian_k  = mdata.second_derivative(0,k);
      const Point<spacedim1> &point_k = points[k];

      for (unsigned int j=0; j<dim1; ++j)
        {
          DF[j] += grad_phi_k[j] * point_k;
          for (unsigned int l=0; l<dim1; ++l)
            D2F[j][l] += hessian_k[j][l] * point_k;
        }
    }

  p_minus_F = p;
  p_minus_F -= transform_unit_to_real_cell_internal(mdata);


  for (unsigned int j=0; j<dim1; ++j)
    f[j] = DF[j] * p_minus_F;

  for (unsigned int j=0; j<dim1; ++j)
    {
      f[j] = DF[j] * p_minus_F;
      for (unsigned int l=0; l<dim1; ++l)
        df[j][l] = -DF[j]*DF[l] + D2F[j][l] * p_minus_F;
    }


  const double eps = 1.e-12*cell->diameter();
  const unsigned int loop_limit = 10;

  unsigned int loop=0;

  while (f.norm()>eps && loop++<loop_limit)
    {
      // Solve  [df(x)]d=f(x)
      Tensor<1,dim1> d;
      Tensor<2,dim1> df_1;

      df_1 = invert(df);
      contract (d, df_1, static_cast<const Tensor<1,dim1>&>(f));
      p_unit -= d;

      for (unsigned int j=0; j<dim1; ++j)
        {
          DF[j].clear();
          for (unsigned int l=0; l<dim1; ++l)
            D2F[j][l].clear();
        }

      mdata.compute_shape_function_values(std::vector<Point<dim> > (1, p_unit));

      for (unsigned int k=0; k<mdata.n_shape_functions; ++k)
        {
          const Tensor<1,dim1>   &grad_phi_k = mdata.derivative(0,k);
          const Tensor<2,dim1>   &hessian_k  = mdata.second_derivative(0,k);
          const Point<spacedim1> &point_k = points[k];

          for (unsigned int j=0; j<dim1; ++j)
            {
              DF[j] += grad_phi_k[j] * point_k;
              for (unsigned int l=0; l<dim1; ++l)
                D2F[j][l] += hessian_k[j][l] * point_k;
            }
        }

//TODO: implement a line search here in much the same way as for
// the corresponding function above that does so for codim==0.
      p_minus_F = p;
      p_minus_F -= transform_unit_to_real_cell_internal(mdata);

      for (unsigned int j=0; j<dim1; ++j)
        {
          f[j] = DF[j] * p_minus_F;
          for (unsigned int l=0; l<dim1; ++l)
            df[j][l] = -DF[j]*DF[l] + D2F[j][l] * p_minus_F;
        }

    }


  // Here we check that in the last
  // execution of while the first
  // condition was already wrong,
  // meaning the residual was below
  // eps. Only if the first condition
  // failed, loop will have been
  // increased and tested, and thus
  // have reached the limit.
  AssertThrow (loop<loop_limit, (typename Mapping<dim,spacedim>::ExcTransformationFailed()));

  return p_unit;
}





template<int dim, int spacedim>
Mapping<dim,spacedim> *
MappingQ1<dim,spacedim>::clone () const
{
  return new MappingQ1<dim,spacedim>(*this);
}

//---------------------------------------------------------------------------


template<int dim, int spacedim>
MappingQ1<dim,spacedim> StaticMappingQ1<dim,spacedim>::mapping;



//--------------------------- Explicit instantiations -----------------------
#include "mapping_q1.inst"


DEAL_II_NAMESPACE_CLOSE
