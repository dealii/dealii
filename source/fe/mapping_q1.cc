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



/**
 * Compute an initial guess to pass to the Newton method in
 * transform_real_to_unit_cell.  For the initial guess we proceed in the
 * following way:
 * <ul>
 * <li> find the least square dim-dimensional plane approximating the cell
 * vertices, i.e. we find an affine map A x_hat + b from the reference cell
 * to the real space.
 * <li> Solve the equation A x_hat + b = p for x_hat
 * <li> This x_hat is the initial solution used for the Newton Method.
 * </ul>
 *
 * @note if dim<spacedim we first project p onto the plane.
 *
 * @note if dim==1 (for any spacedim) the initial guess is the exact
 * solution and no Newton iteration is needed.
 *
 * Some details about how we compute the least square plane. We look
 * for a spacedim x (dim + 1) matrix X such that X * M = Y where M is
 * a (dim+1) x n_vertices matrix and Y a spacedim x n_vertices.  And:
 * The i-th column of M is unit_vertex[i] and the last row all
 * 1's. The i-th column of Y is real_vertex[i].  If we split X=[A|b],
 * the least square approx is A x_hat+b Classically X = Y * (M^t (M
 * M^t)^{-1}) Let K = M^t * (M M^t)^{-1} = [KA Kb] this can be
 * precomputed, and that is exactly what we do.  Finally A = Y*KA and
 * b = Y*Kb.
 */
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

  template<int dim, int spacedim>
  Point<dim>
  transform_real_to_unit_cell_initial_guess (const std::vector<Point<spacedim> > &vertex,
                                             const Point<spacedim>               &p)
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

  // Find the initial value for the Newton iteration by a normal
  // projection to the least square plane determined by the vertices
  // of the cell
  std::vector<Point<spacedim> > a;
  compute_mapping_support_points (cell,a);
  Point<dim> initial_p_unit =
    transform_real_to_unit_cell_initial_guess<dim,spacedim>(a,p);

  // if dim==1 there is nothing else to do to the initial value, and
  // it is the answer
  if (dim == 1)
    return initial_p_unit;
  else
    {
      // use the full mapping. in case the function above should have
      // given us something back that lies outside the unit cell (that
      // might happen because either the function computing an initial
      // guess gave us a poor initial guess or for the following
      // reason: we call this function here in the Q1 mapping to
      // produce an initial guess for a higher order mapping, but we
      // may have given a point 'p' that lies inside the cell with the
      // higher order mapping, but outside the Q1-mapped reference
      // cell), then project it back into the reference cell in hopes
      // that this gives a better starting point to the following
      // iteration
      const Quadrature<dim> point_quadrature(initial_p_unit);

      UpdateFlags update_flags = update_quadrature_points | update_jacobians;
      if (spacedim>dim)
        update_flags |= update_jacobian_grads;

      //TODO: Use get_data() here once MappingQ is no longer derived from
      //MappingQ1. this doesn't currently work because we here really need
      //a Q1 InternalData, but MappingQGeneric produces one with the
      //polynomial degree of the MappingQ
      std_cxx11::unique_ptr<typename MappingQGeneric<dim,spacedim>::InternalData>
      mdata (new typename MappingQGeneric<dim,spacedim>::InternalData(1));
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
      return this->transform_real_to_unit_cell_internal(cell, p, initial_p_unit,
                                                        *mdata);
    }
}





template<int dim, int spacedim>
MappingQ1<dim,spacedim> *
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
