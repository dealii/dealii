// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_mapping_q_internal_h
#define dealii_mapping_q_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/fe_evaluation_data.h>
#include <deal.II/matrix_free/mapping_info_storage.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <array>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * Internal namespace to implement methods specific to MappingQ1, in
   * particular an explicit formula for the transformation from the real to
   * the unit cell in 2D.
   */
  namespace MappingQ1
  {
    // These are left as templates on the spatial dimension (even though dim
    // == spacedim must be true for them to make sense) because templates are
    // expanded before the compiler eliminates code due to the 'if (dim ==
    // spacedim)' statement (see the body of the general
    // transform_real_to_unit_cell).
    template <int spacedim>
    inline Point<1>
    transform_real_to_unit_cell(
      const std::array<Point<spacedim>, GeometryInfo<1>::vertices_per_cell>
        &                    vertices,
      const Point<spacedim> &p)
    {
      Assert(spacedim == 1, ExcInternalError());
      return Point<1>((p[0] - vertices[0](0)) /
                      (vertices[1](0) - vertices[0](0)));
    }



    template <int spacedim>
    inline Point<2>
    transform_real_to_unit_cell(
      const std::array<Point<spacedim>, GeometryInfo<2>::vertices_per_cell>
        &                    vertices,
      const Point<spacedim> &p)
    {
      Assert(spacedim == 2, ExcInternalError());

      // For accuracy reasons, we do all arithmetic in extended precision
      // (long double). This has a noticeable effect on the hit rate for
      // borderline cases and thus makes the algorithm more robust.
      const long double x = p(0);
      const long double y = p(1);

      const long double x0 = vertices[0](0);
      const long double x1 = vertices[1](0);
      const long double x2 = vertices[2](0);
      const long double x3 = vertices[3](0);

      const long double y0 = vertices[0](1);
      const long double y1 = vertices[1](1);
      const long double y2 = vertices[2](1);
      const long double y3 = vertices[3](1);

      const long double a = (x1 - x3) * (y0 - y2) - (x0 - x2) * (y1 - y3);
      const long double b = -(x0 - x1 - x2 + x3) * y + (x - 2 * x1 + x3) * y0 -
                            (x - 2 * x0 + x2) * y1 - (x - x1) * y2 +
                            (x - x0) * y3;
      const long double c = (x0 - x1) * y - (x - x1) * y0 + (x - x0) * y1;

      const long double discriminant = b * b - 4 * a * c;
      // exit if the point is not in the cell (this is the only case where the
      // discriminant is negative)
      AssertThrow(
        discriminant > 0.0,
        (typename Mapping<spacedim, spacedim>::ExcTransformationFailed()));

      long double       eta1;
      long double       eta2;
      const long double sqrt_discriminant = std::sqrt(discriminant);
      // special case #1: if a is near-zero to make the discriminant exactly
      // equal b, then use the linear formula
      if (b != 0.0 && std::abs(b) == sqrt_discriminant)
        {
          eta1 = -c / b;
          eta2 = -c / b;
        }
      // special case #2: a is zero for parallelograms and very small for
      // near-parallelograms:
      else if (std::abs(a) < 1e-8 * std::abs(b))
        {
          // if both a and c are very small then the root should be near
          // zero: this first case will capture that
          eta1 = 2 * c / (-b - sqrt_discriminant);
          eta2 = 2 * c / (-b + sqrt_discriminant);
        }
      // finally, use the plain version:
      else
        {
          eta1 = (-b - sqrt_discriminant) / (2 * a);
          eta2 = (-b + sqrt_discriminant) / (2 * a);
        }
      // pick the one closer to the center of the cell.
      const long double eta =
        (std::abs(eta1 - 0.5) < std::abs(eta2 - 0.5)) ? eta1 : eta2;

      /*
       * There are two ways to compute xi from eta, but either one may have a
       * zero denominator.
       */
      const long double subexpr0        = -eta * x2 + x0 * (eta - 1);
      const long double xi_denominator0 = eta * x3 - x1 * (eta - 1) + subexpr0;
      const long double max_x = std::max(std::max(std::abs(x0), std::abs(x1)),
                                         std::max(std::abs(x2), std::abs(x3)));

      if (std::abs(xi_denominator0) > 1e-10 * max_x)
        {
          const double xi = (x + subexpr0) / xi_denominator0;
          return {xi, static_cast<double>(eta)};
        }
      else
        {
          const long double max_y =
            std::max(std::max(std::abs(y0), std::abs(y1)),
                     std::max(std::abs(y2), std::abs(y3)));
          const long double subexpr1 = -eta * y2 + y0 * (eta - 1);
          const long double xi_denominator1 =
            eta * y3 - y1 * (eta - 1) + subexpr1;
          if (std::abs(xi_denominator1) > 1e-10 * max_y)
            {
              const double xi = (subexpr1 + y) / xi_denominator1;
              return {xi, static_cast<double>(eta)};
            }
          else // give up and try Newton iteration
            {
              AssertThrow(
                false,
                (typename Mapping<spacedim,
                                  spacedim>::ExcTransformationFailed()));
            }
        }
      // bogus return to placate compiler. It should not be possible to get
      // here.
      Assert(false, ExcInternalError());
      return {std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN()};
    }



    template <int spacedim>
    inline Point<3>
    transform_real_to_unit_cell(
      const std::array<Point<spacedim>, GeometryInfo<3>::vertices_per_cell>
        & /*vertices*/,
      const Point<spacedim> & /*p*/)
    {
      // It should not be possible to get here
      Assert(false, ExcInternalError());
      return {std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN()};
    }
  } // namespace MappingQ1



  /**
   * Internal namespace to implement methods of MappingQ, such as the
   * evaluation of the mapping and the transformation between real and unit
   * cell.
   */
  namespace MappingQImplementation
  {
    /**
     * This function generates the reference cell support points from the 1d
     * support points by expanding the tensor product.
     */
    template <int dim>
    std::vector<Point<dim>>
    unit_support_points(const std::vector<Point<1>> &    line_support_points,
                        const std::vector<unsigned int> &renumbering)
    {
      AssertDimension(Utilities::pow(line_support_points.size(), dim),
                      renumbering.size());
      std::vector<Point<dim>> points(renumbering.size());
      const unsigned int      n1 = line_support_points.size();
      for (unsigned int q2 = 0, q = 0; q2 < (dim > 2 ? n1 : 1); ++q2)
        for (unsigned int q1 = 0; q1 < (dim > 1 ? n1 : 1); ++q1)
          for (unsigned int q0 = 0; q0 < n1; ++q0, ++q)
            {
              points[renumbering[q]][0] = line_support_points[q0][0];
              if (dim > 1)
                points[renumbering[q]][1] = line_support_points[q1][0];
              if (dim > 2)
                points[renumbering[q]][2] = line_support_points[q2][0];
            }
      return points;
    }



    /**
     * This function is needed by the constructor of
     * <tt>MappingQ<dim,spacedim></tt> for <tt>dim=</tt> 2 and 3.
     *
     * For the definition of the @p support_point_weights_on_quad please
     * refer to the description of TransfiniteInterpolationManifold.
     */
    inline dealii::Table<2, double>
    compute_support_point_weights_on_quad(const unsigned int polynomial_degree)
    {
      dealii::Table<2, double> loqvs;

      // we are asked to compute weights for interior support points, but
      // there are no interior points if degree==1
      if (polynomial_degree == 1)
        return loqvs;

      const unsigned int M          = polynomial_degree - 1;
      const unsigned int n_inner_2d = M * M;
      const unsigned int n_outer_2d = 4 + 4 * M;

      // set the weights of transfinite interpolation
      loqvs.reinit(n_inner_2d, n_outer_2d);
      QGaussLobatto<2> gl(polynomial_degree + 1);
      for (unsigned int i = 0; i < M; ++i)
        for (unsigned int j = 0; j < M; ++j)
          {
            const Point<2> &p =
              gl.point((i + 1) * (polynomial_degree + 1) + (j + 1));
            const unsigned int index_table = i * M + j;
            for (unsigned int v = 0; v < 4; ++v)
              loqvs(index_table, v) =
                -GeometryInfo<2>::d_linear_shape_function(p, v);
            loqvs(index_table, 4 + i)         = 1. - p[0];
            loqvs(index_table, 4 + i + M)     = p[0];
            loqvs(index_table, 4 + j + 2 * M) = 1. - p[1];
            loqvs(index_table, 4 + j + 3 * M) = p[1];
          }

      // the sum of weights of the points at the outer rim should be one.
      // check this
      for (unsigned int unit_point = 0; unit_point < n_inner_2d; ++unit_point)
        Assert(std::fabs(std::accumulate(loqvs[unit_point].begin(),
                                         loqvs[unit_point].end(),
                                         0.) -
                         1) < 1e-13 * polynomial_degree,
               ExcInternalError());

      return loqvs;
    }



    /**
     * This function is needed by the constructor of <tt>MappingQ<3></tt>.
     *
     * For the definition of the @p support_point_weights_on_quad please
     * refer to the description of TransfiniteInterpolationManifold.
     */
    inline dealii::Table<2, double>
    compute_support_point_weights_on_hex(const unsigned int polynomial_degree)
    {
      dealii::Table<2, double> lohvs;

      // we are asked to compute weights for interior support points, but
      // there are no interior points if degree==1
      if (polynomial_degree == 1)
        return lohvs;

      const unsigned int M = polynomial_degree - 1;

      const unsigned int n_inner = Utilities::fixed_power<3>(M);
      const unsigned int n_outer = 8 + 12 * M + 6 * M * M;

      // set the weights of transfinite interpolation
      lohvs.reinit(n_inner, n_outer);
      QGaussLobatto<3> gl(polynomial_degree + 1);
      for (unsigned int i = 0; i < M; ++i)
        for (unsigned int j = 0; j < M; ++j)
          for (unsigned int k = 0; k < M; ++k)
            {
              const Point<3> &   p = gl.point((i + 1) * (M + 2) * (M + 2) +
                                           (j + 1) * (M + 2) + (k + 1));
              const unsigned int index_table = i * M * M + j * M + k;

              // vertices
              for (unsigned int v = 0; v < 8; ++v)
                lohvs(index_table, v) =
                  GeometryInfo<3>::d_linear_shape_function(p, v);

              // lines
              {
                constexpr std::array<unsigned int, 4> line_coordinates_y(
                  {{0, 1, 4, 5}});
                const Point<2> py(p[0], p[2]);
                for (unsigned int l = 0; l < 4; ++l)
                  lohvs(index_table, 8 + line_coordinates_y[l] * M + j) =
                    -GeometryInfo<2>::d_linear_shape_function(py, l);
              }

              {
                constexpr std::array<unsigned int, 4> line_coordinates_x(
                  {{2, 3, 6, 7}});
                const Point<2> px(p[1], p[2]);
                for (unsigned int l = 0; l < 4; ++l)
                  lohvs(index_table, 8 + line_coordinates_x[l] * M + k) =
                    -GeometryInfo<2>::d_linear_shape_function(px, l);
              }

              {
                constexpr std::array<unsigned int, 4> line_coordinates_z(
                  {{8, 9, 10, 11}});
                const Point<2> pz(p[0], p[1]);
                for (unsigned int l = 0; l < 4; ++l)
                  lohvs(index_table, 8 + line_coordinates_z[l] * M + i) =
                    -GeometryInfo<2>::d_linear_shape_function(pz, l);
              }

              // quads
              lohvs(index_table, 8 + 12 * M + 0 * M * M + i * M + j) =
                1. - p[0];
              lohvs(index_table, 8 + 12 * M + 1 * M * M + i * M + j) = p[0];
              lohvs(index_table, 8 + 12 * M + 2 * M * M + k * M + i) =
                1. - p[1];
              lohvs(index_table, 8 + 12 * M + 3 * M * M + k * M + i) = p[1];
              lohvs(index_table, 8 + 12 * M + 4 * M * M + j * M + k) =
                1. - p[2];
              lohvs(index_table, 8 + 12 * M + 5 * M * M + j * M + k) = p[2];
            }

      // the sum of weights of the points at the outer rim should be one.
      // check this
      for (unsigned int unit_point = 0; unit_point < n_inner; ++unit_point)
        Assert(std::fabs(std::accumulate(lohvs[unit_point].begin(),
                                         lohvs[unit_point].end(),
                                         0.) -
                         1) < 1e-13 * polynomial_degree,
               ExcInternalError());

      return lohvs;
    }



    /**
     * This function collects the output of
     * compute_support_point_weights_on_{quad,hex} in a single data structure.
     */
    inline std::vector<dealii::Table<2, double>>
    compute_support_point_weights_perimeter_to_interior(
      const unsigned int polynomial_degree,
      const unsigned int dim)
    {
      Assert(dim > 0 && dim <= 3, ExcImpossibleInDim(dim));
      std::vector<dealii::Table<2, double>> output(dim);
      if (polynomial_degree <= 1)
        return output;

      // fill the 1D interior weights
      QGaussLobatto<1> quadrature(polynomial_degree + 1);
      output[0].reinit(polynomial_degree - 1,
                       GeometryInfo<1>::vertices_per_cell);
      for (unsigned int q = 0; q < polynomial_degree - 1; ++q)
        for (const unsigned int i : GeometryInfo<1>::vertex_indices())
          output[0](q, i) =
            GeometryInfo<1>::d_linear_shape_function(quadrature.point(q + 1),
                                                     i);

      if (dim > 1)
        output[1] = compute_support_point_weights_on_quad(polynomial_degree);

      if (dim > 2)
        output[2] = compute_support_point_weights_on_hex(polynomial_degree);

      return output;
    }



    /**
     * Collects all interior points for the various dimensions.
     */
    template <int dim>
    inline dealii::Table<2, double>
    compute_support_point_weights_cell(const unsigned int polynomial_degree)
    {
      Assert(dim > 0 && dim <= 3, ExcImpossibleInDim(dim));
      if (polynomial_degree <= 1)
        return dealii::Table<2, double>();

      QGaussLobatto<dim>              quadrature(polynomial_degree + 1);
      const std::vector<unsigned int> h2l =
        FETools::hierarchic_to_lexicographic_numbering<dim>(polynomial_degree);

      dealii::Table<2, double> output(quadrature.size() -
                                        GeometryInfo<dim>::vertices_per_cell,
                                      GeometryInfo<dim>::vertices_per_cell);
      for (unsigned int q = 0; q < output.size(0); ++q)
        for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
          output(q, i) = GeometryInfo<dim>::d_linear_shape_function(
            quadrature.point(h2l[q + GeometryInfo<dim>::vertices_per_cell]), i);

      return output;
    }



    /**
     * Using the relative weights of the shape functions evaluated at
     * one point on the reference cell (and stored in data.shape_values
     * and accessed via data.shape(0,i)) and the locations of mapping
     * support points (stored in data.mapping_support_points), compute
     * the mapped location of that point in real space.
     */
    template <int dim, int spacedim>
    inline Point<spacedim>
    compute_mapped_location_of_point(
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data)
    {
      AssertDimension(data.shape_values.size(),
                      data.mapping_support_points.size());

      // use now the InternalData to compute the point in real space.
      Point<spacedim> p_real;
      for (unsigned int i = 0; i < data.mapping_support_points.size(); ++i)
        p_real += data.mapping_support_points[i] * data.shape(0, i);

      return p_real;
    }



    /**
     * Implementation of transform_real_to_unit_cell for either type double
     * or VectorizedArray<double>
     */
    template <int dim, int spacedim, typename Number>
    inline Point<dim, Number>
    do_transform_real_to_unit_cell_internal(
      const Point<spacedim, Number> &                     p,
      const Point<dim, Number> &                          initial_p_unit,
      const std::vector<Point<spacedim>> &                points,
      const std::vector<Polynomials::Polynomial<double>> &polynomials_1d,
      const std::vector<unsigned int> &                   renumber,
      const bool print_iterations_to_deallog = false)
    {
      if (print_iterations_to_deallog)
        deallog << "Start MappingQ::do_transform_real_to_unit_cell for real "
                << "point [ " << p << " ] " << std::endl;

      AssertDimension(points.size(),
                      Utilities::pow(polynomials_1d.size(), dim));

      // Newton iteration to solve
      //    f(x)=p(x)-p=0
      // where we are looking for 'x' and p(x) is the forward transformation
      // from unit to real cell. We solve this using a Newton iteration
      //    x_{n+1}=x_n-[f'(x)]^{-1}f(x)
      // The start value is set to be the linear approximation to the cell

      // The shape values and derivatives of the mapping at this point are
      // previously computed.

      Point<dim, Number> p_unit = initial_p_unit;
      auto p_real = internal::evaluate_tensor_product_value_and_gradient(
        polynomials_1d, points, p_unit, polynomials_1d.size() == 2, renumber);

      Tensor<1, spacedim, Number> f = p_real.first - p;

      // early out if we already have our point in all SIMD lanes, i.e.,
      // f.norm_square() < 1e-24 * p_real.second[0].norm_square(). To enable
      // this step also for VectorizedArray where we do not have operator <,
      // compare the result to zero which is defined for SIMD types
      if (std::max(Number(0.),
                   f.norm_square() - 1e-24 * p_real.second[0].norm_square()) ==
          Number(0.))
        return p_unit;

      // we need to compare the position of the computed p(x) against the
      // given point 'p'. We will terminate the iteration and return 'x' if
      // they are less than eps apart. The question is how to choose eps --
      // or, put maybe more generally: in which norm we want these 'p' and
      // 'p(x)' to be eps apart.
      //
      // the question is difficult since we may have to deal with very
      // elongated cells where we may achieve 1e-12*h for the distance of
      // these two points in the 'long' direction, but achieving this
      // tolerance in the 'short' direction of the cell may not be possible
      //
      // what we do instead is then to terminate iterations if
      //    \| p(x) - p \|_A < eps
      // where the A-norm is somehow induced by the transformation of the
      // cell. in particular, we want to measure distances relative to the
      // sizes of the cell in its principal directions.
      //
      // to define what exactly A should be, note that to first order we have
      // the following (assuming that x* is the solution of the problem, i.e.,
      // p(x*)=p):
      //    p(x) - p = p(x) - p(x*)
      //             = -grad p(x) * (x*-x) + higher order terms
      // This suggest to measure with a norm that corresponds to
      //    A = {[grad p(x)]^T [grad p(x)]}^{-1}
      // because then
      //    \| p(x) - p \|_A  \approx  \| x - x* \|
      // Consequently, we will try to enforce that
      //    \| p(x) - p \|_A  =  \| f \|  <=  eps
      //
      // Note that using this norm is a bit dangerous since the norm changes
      // in every iteration (A isn't fixed by depending on xk). However, if
      // the cell is not too deformed (it may be stretched, but not twisted)
      // then the mapping is almost linear and A is indeed constant or
      // nearly so.
      const double       eps                    = 1.e-11;
      const unsigned int newton_iteration_limit = 20;

      Point<dim, Number> invalid_point;
      invalid_point[0]                = std::numeric_limits<double>::infinity();
      bool tried_project_to_unit_cell = false;

      unsigned int newton_iteration            = 0;
      Number       f_weighted_norm_square      = 1.;
      Number       last_f_weighted_norm_square = 1.;

      do
        {
          if (print_iterations_to_deallog)
            deallog << "Newton iteration " << newton_iteration
                    << " for unit point guess " << p_unit << std::endl;

          // f'(x)
          Tensor<2, spacedim, Number> df;
          for (unsigned int d = 0; d < spacedim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
              df[d][e] = p_real.second[e][d];

          // check determinand(df) > 0 on all SIMD lanes
          if (!(std::min(determinant(df),
                         Number(std::numeric_limits<double>::min())) ==
                Number(std::numeric_limits<double>::min())))
            {
              // We allow to enter this function with an initial guess
              // outside the unit cell. We might have run into invalid
              // Jacobians due to that, so we should at least try once to go
              // back to the unit cell and go on with the Newton iteration
              // from there. Since the outside case is unlikely, we can
              // afford spending the extra effort at this place.
              if (tried_project_to_unit_cell == false)
                {
                  p_unit = GeometryInfo<dim>::project_to_unit_cell(p_unit);
                  p_real = internal::evaluate_tensor_product_value_and_gradient(
                    polynomials_1d,
                    points,
                    p_unit,
                    polynomials_1d.size() == 2,
                    renumber);
                  f                           = p_real.first - p;
                  f_weighted_norm_square      = 1.;
                  last_f_weighted_norm_square = 1;
                  tried_project_to_unit_cell  = true;
                  continue;
                }
              else
                return invalid_point;
            }

          // Solve  [f'(x)]d=f(x)
          const Tensor<2, spacedim, Number> df_inverse = invert(df);
          const Tensor<1, spacedim, Number> delta      = df_inverse * f;
          last_f_weighted_norm_square                  = delta.norm_square();

          if (print_iterations_to_deallog)
            deallog << "   delta=" << delta << std::endl;

          // do a line search
          double step_length = 1.0;
          do
            {
              // update of p_unit. The spacedim-th component of transformed
              // point is simply ignored in codimension one case. When this
              // component is not zero, then we are projecting the point to
              // the surface or curve identified by the cell.
              Point<dim, Number> p_unit_trial = p_unit;
              for (unsigned int i = 0; i < dim; ++i)
                p_unit_trial[i] -= step_length * delta[i];

              // shape values and derivatives at new p_unit point
              const auto p_real_trial =
                internal::evaluate_tensor_product_value_and_gradient(
                  polynomials_1d,
                  points,
                  p_unit_trial,
                  polynomials_1d.size() == 2,
                  renumber);
              const Tensor<1, spacedim, Number> f_trial =
                p_real_trial.first - p;
              f_weighted_norm_square = (df_inverse * f_trial).norm_square();

              if (print_iterations_to_deallog)
                {
                  deallog << "     step_length=" << step_length << std::endl;
                  if (step_length == 1.0)
                    deallog << "       ||f ||   =" << f.norm() << std::endl;
                  deallog << "       ||f*||   =" << f_trial.norm() << std::endl
                          << "       ||f*||_A ="
                          << std::sqrt(f_weighted_norm_square) << std::endl;
                }

              // See if we are making progress with the current step length
              // and if not, reduce it by a factor of two and try again.
              //
              // Strictly speaking, we should probably use the same norm as we
              // use for the outer algorithm. In practice, line search is just
              // a crutch to find a "reasonable" step length, and so using the
              // l2 norm is probably just fine.
              //
              // check f_trial.norm() < f.norm() in SIMD form. This is a bit
              // more complicated because some SIMD lanes might not be doing
              // any progress any more as they have already reached roundoff
              // accuracy. We define that as the case of updates less than
              // 1e-6. The tolerance might seem coarse but since we are
              // dealing with a Newton iteration of a polynomial function we
              // either converge quadratically or not any more. Thus, our
              // selection is to terminate if either the norm of f is
              // decreasing or that threshold of 1e-6 is reached.
              if (std::max(f_weighted_norm_square - 1e-6 * 1e-6, Number(0.)) *
                    std::max(f_trial.norm_square() - f.norm_square(),
                             Number(0.)) ==
                  Number(0.))
                {
                  p_real = p_real_trial;
                  p_unit = p_unit_trial;
                  f      = f_trial;
                  break;
                }
              else if (step_length > 0.05)
                step_length *= 0.5;
              else
                break;
            }
          while (true);

          // In case we terminated the line search due to the step becoming
          // too small, we give the iteration another try with the
          // projection of the initial guess to the unit cell before we give
          // up, just like for the negative determinant case.
          if (step_length <= 0.05 && tried_project_to_unit_cell == false)
            {
              p_unit = GeometryInfo<dim>::project_to_unit_cell(p_unit);
              p_real = internal::evaluate_tensor_product_value_and_gradient(
                polynomials_1d,
                points,
                p_unit,
                polynomials_1d.size() == 2,
                renumber);
              f                           = p_real.first - p;
              f_weighted_norm_square      = 1.;
              last_f_weighted_norm_square = 1;
              tried_project_to_unit_cell  = true;
              continue;
            }
          else if (step_length <= 0.05)
            return invalid_point;

          ++newton_iteration;
          if (newton_iteration > newton_iteration_limit)
            return invalid_point;
        }
      // Stop if f_weighted_norm_square <= eps^2 on all SIMD lanes or if the
      // weighted norm is less than 1e-6 and the improvement against the
      // previous step was less than a factor of 10 (in that regime, we
      // either have quadratic convergence or roundoff errors due to a bad
      // mapping)
      while (
        !(std::max(f_weighted_norm_square - eps * eps, Number(0.)) *
            std::max(last_f_weighted_norm_square -
                       std::min(f_weighted_norm_square, Number(1e-6 * 1e-6)) *
                         100.,
                     Number(0.)) ==
          Number(0.)));

      if (print_iterations_to_deallog)
        deallog << "Iteration converged for p_unit = [ " << p_unit
                << " ] and iteration error "
                << std::sqrt(f_weighted_norm_square) << std::endl;

      return p_unit;
    }



    /**
     * Implementation of transform_real_to_unit_cell for dim==spacedim-1
     */
    template <int dim>
    inline Point<dim>
    do_transform_real_to_unit_cell_internal_codim1(
      const typename dealii::Triangulation<dim, dim + 1>::cell_iterator &cell,
      const Point<dim + 1> &                                             p,
      const Point<dim> &                                     initial_p_unit,
      typename dealii::MappingQ<dim, dim + 1>::InternalData &mdata)
    {
      const unsigned int spacedim = dim + 1;

      const unsigned int n_shapes = mdata.shape_values.size();
      (void)n_shapes;
      Assert(n_shapes != 0, ExcInternalError());
      Assert(mdata.shape_derivatives.size() == n_shapes, ExcInternalError());
      Assert(mdata.shape_second_derivatives.size() == n_shapes,
             ExcInternalError());

      std::vector<Point<spacedim>> &points = mdata.mapping_support_points;
      Assert(points.size() == n_shapes, ExcInternalError());

      Point<spacedim> p_minus_F;

      Tensor<1, spacedim> DF[dim];
      Tensor<1, spacedim> D2F[dim][dim];

      Point<dim>     p_unit = initial_p_unit;
      Point<dim>     f;
      Tensor<2, dim> df;

      // Evaluate first and second derivatives
      mdata.compute_shape_function_values(std::vector<Point<dim>>(1, p_unit));

      for (unsigned int k = 0; k < mdata.n_shape_functions; ++k)
        {
          const Tensor<1, dim> & grad_phi_k = mdata.derivative(0, k);
          const Tensor<2, dim> & hessian_k  = mdata.second_derivative(0, k);
          const Point<spacedim> &point_k    = points[k];

          for (unsigned int j = 0; j < dim; ++j)
            {
              DF[j] += grad_phi_k[j] * point_k;
              for (unsigned int l = 0; l < dim; ++l)
                D2F[j][l] += hessian_k[j][l] * point_k;
            }
        }

      p_minus_F = p;
      p_minus_F -= compute_mapped_location_of_point<dim, spacedim>(mdata);


      for (unsigned int j = 0; j < dim; ++j)
        f[j] = DF[j] * p_minus_F;

      for (unsigned int j = 0; j < dim; ++j)
        {
          f[j] = DF[j] * p_minus_F;
          for (unsigned int l = 0; l < dim; ++l)
            df[j][l] = -DF[j] * DF[l] + D2F[j][l] * p_minus_F;
        }


      const double       eps        = 1.e-12 * cell->diameter();
      const unsigned int loop_limit = 10;

      unsigned int loop = 0;

      while (f.norm() > eps && loop++ < loop_limit)
        {
          // Solve  [df(x)]d=f(x)
          const Tensor<1, dim> d =
            invert(df) * static_cast<const Tensor<1, dim> &>(f);
          p_unit -= d;

          for (unsigned int j = 0; j < dim; ++j)
            {
              DF[j].clear();
              for (unsigned int l = 0; l < dim; ++l)
                D2F[j][l].clear();
            }

          mdata.compute_shape_function_values(
            std::vector<Point<dim>>(1, p_unit));

          for (unsigned int k = 0; k < mdata.n_shape_functions; ++k)
            {
              const Tensor<1, dim> & grad_phi_k = mdata.derivative(0, k);
              const Tensor<2, dim> & hessian_k  = mdata.second_derivative(0, k);
              const Point<spacedim> &point_k    = points[k];

              for (unsigned int j = 0; j < dim; ++j)
                {
                  DF[j] += grad_phi_k[j] * point_k;
                  for (unsigned int l = 0; l < dim; ++l)
                    D2F[j][l] += hessian_k[j][l] * point_k;
                }
            }

          // TODO: implement a line search here in much the same way as for
          // the corresponding function above that does so for dim==spacedim
          p_minus_F = p;
          p_minus_F -= compute_mapped_location_of_point<dim, spacedim>(mdata);

          for (unsigned int j = 0; j < dim; ++j)
            {
              f[j] = DF[j] * p_minus_F;
              for (unsigned int l = 0; l < dim; ++l)
                df[j][l] = -DF[j] * DF[l] + D2F[j][l] * p_minus_F;
            }
        }


      // Here we check that in the last execution of while the first
      // condition was already wrong, meaning the residual was below
      // eps. Only if the first condition failed, loop will have been
      // increased and tested, and thus have reached the limit.
      AssertThrow(loop < loop_limit,
                  (typename Mapping<dim, spacedim>::ExcTransformationFailed()));

      return p_unit;
    }



    /**
     * A class to compute a quadratic approximation to the inverse map from
     * real to unit points by a least-squares fit along the mapping support
     * points. The least squares fit is special in the sense that the
     * approximation is constructed for the inverse function of a
     * MappingQ, which is generally a rational function. This allows
     * for a very cheap evaluation of the inverse map by a simple polynomial
     * interpolation, which can be used as a better initial guess for
     * transforming points from real to unit coordinates than an affine
     * approximation.
     *
     * Far away outside the unit cell, this approximation can become
     * inaccurate for non-affine cell shapes. This must be expected from a
     * fit of a polynomial to a rational function, and due to the fact that
     * the region of the least squares fit, the unit cell, is left. Hence,
     * use this function with care in those situations.
     */
    template <int dim, int spacedim>
    class InverseQuadraticApproximation
    {
    public:
      /**
       * Number of basis functions in the quadratic approximation.
       */
      static constexpr unsigned int n_functions =
        (spacedim == 1 ? 3 : (spacedim == 2 ? 6 : 10));

      /**
       * Constructor.
       *
       * @param real_support_points The position of the mapping support points
       * in real space, queried by
       * MappingQ::compute_mapping_support_points().
       *
       * @param unit_support_points The location of the support points in
       * reference coordinates $[0, 1]^d$ that map to the mapping support
       * points in real space by a polynomial map.
       */
      InverseQuadraticApproximation(
        const std::vector<Point<spacedim>> &real_support_points,
        const std::vector<Point<dim>> &     unit_support_points)
        : normalization_shift(real_support_points[0])
        , normalization_length(
            1. / real_support_points[0].distance(real_support_points[1]))
        , is_affine(true)
      {
        AssertDimension(real_support_points.size(), unit_support_points.size());

        // For the bi-/trilinear approximation, we cannot build a quadratic
        // polynomial due to a lack of points (interpolation matrix would
        // get singular). Similarly, it is not entirely clear how to gather
        // enough information for the case dim < spacedim.
        //
        // In both cases we require the vector real_support_points to
        // contain the vertex positions and fall back to an affine
        // approximation:
        Assert(dim == spacedim || real_support_points.size() ==
                                    GeometryInfo<dim>::vertices_per_cell,
               ExcInternalError());
        if (real_support_points.size() == GeometryInfo<dim>::vertices_per_cell)
          {
            const auto affine = GridTools::affine_cell_approximation<dim>(
              make_array_view(real_support_points));
            DerivativeForm<1, spacedim, dim> A_inv =
              affine.first.covariant_form().transpose();

            // The code for evaluation assumes an additional transformation of
            // the form (x - normalization_shift) * normalization_length --
            // account for this in the definition of the coefficients.
            coefficients[0] =
              apply_transformation(A_inv, normalization_shift - affine.second);
            for (unsigned int d = 0; d < spacedim; ++d)
              for (unsigned int e = 0; e < dim; ++e)
                coefficients[1 + d][e] =
                  A_inv[e][d] * (1.0 / normalization_length);
            is_affine = true;
            return;
          }

        Tensor<2, n_functions>          matrix;
        std::array<double, n_functions> shape_values;
        for (unsigned int q = 0; q < unit_support_points.size(); ++q)
          {
            // Evaluate quadratic shape functions in point, with the
            // normalization applied in order to avoid roundoff issues with
            // scaling far away from 1.
            shape_values[0] = 1.;
            const Tensor<1, spacedim> p_scaled =
              normalization_length *
              (real_support_points[q] - normalization_shift);
            for (unsigned int d = 0; d < spacedim; ++d)
              shape_values[1 + d] = p_scaled[d];
            for (unsigned int d = 0, c = 0; d < spacedim; ++d)
              for (unsigned int e = 0; e <= d; ++e, ++c)
                shape_values[1 + spacedim + c] = p_scaled[d] * p_scaled[e];

            // Build lower diagonal of least squares matrix and rhs, the
            // essential part being that we construct the matrix with the
            // real points and the right hand side by comparing to the
            // reference point positions which sets up an inverse
            // interpolation.
            for (unsigned int i = 0; i < n_functions; ++i)
              for (unsigned int j = 0; j < n_functions; ++j)
                matrix[i][j] += shape_values[i] * shape_values[j];
            for (unsigned int i = 0; i < n_functions; ++i)
              coefficients[i] += shape_values[i] * unit_support_points[q];
          }

        // Factorize the matrix A = L * L^T in-place with the
        // Cholesky-Banachiewicz algorithm. The implementation is similar to
        // FullMatrix::cholesky() but re-implemented to avoid memory
        // allocations and some unnecessary divisions which we can do here as
        // we only need to solve with dim right hand sides.
        for (unsigned int i = 0; i < n_functions; ++i)
          {
            double Lij_sum = 0;
            for (unsigned int j = 0; j < i; ++j)
              {
                double Lik_Ljk_sum = 0;
                for (unsigned int k = 0; k < j; ++k)
                  Lik_Ljk_sum += matrix[i][k] * matrix[j][k];
                matrix[i][j] = matrix[j][j] * (matrix[i][j] - Lik_Ljk_sum);
                Lij_sum += matrix[i][j] * matrix[i][j];
              }
            AssertThrow(matrix[i][i] - Lij_sum >= 0,
                        ExcMessage("Matrix of normal equations not positive "
                                   "definite"));

            // Store the inverse in the diagonal since that is the quantity
            // needed later in the factorization as well as the forward and
            // backward substitution, minimizing the number of divisions.
            matrix[i][i] = 1. / std::sqrt(matrix[i][i] - Lij_sum);
          }

        // Solve lower triangular part, L * y = rhs.
        for (unsigned int i = 0; i < n_functions; ++i)
          {
            Point<dim> sum = coefficients[i];
            for (unsigned int j = 0; j < i; ++j)
              sum -= matrix[i][j] * coefficients[j];
            coefficients[i] = sum * matrix[i][i];
          }

        // Solve upper triangular part, L^T * x = y (i.e., x = A^{-1} * rhs)
        for (unsigned int i = n_functions; i > 0;)
          {
            --i;
            Point<dim> sum = coefficients[i];
            for (unsigned int j = i + 1; j < n_functions; ++j)
              sum -= matrix[j][i] * coefficients[j];
            coefficients[i] = sum * matrix[i][i];
          }

        // Check whether the approximation is indeed affine, allowing to
        // skip the quadratic terms.
        is_affine = true;
        for (unsigned int i = dim + 1; i < n_functions; ++i)
          if (coefficients[i].norm_square() > 1e-20)
            {
              is_affine = false;
              break;
            }
      }

      /**
       * Copy constructor.
       */
      InverseQuadraticApproximation(const InverseQuadraticApproximation &) =
        default;

      /**
       * Evaluate the quadratic approximation.
       */
      template <typename Number>
      Point<dim, Number>
      compute(const Point<spacedim, Number> &p) const
      {
        Point<dim, Number> result;
        for (unsigned int d = 0; d < dim; ++d)
          result[d] = coefficients[0][d];

        // Apply the normalization to ensure a good conditioning. Since Number
        // might be a vectorized array whereas the normalization is a point of
        // doubles, we cannot use the overload of operator- and must instead
        // loop over the components of the point.
        Point<spacedim, Number> p_scaled;
        for (unsigned int d = 0; d < spacedim; ++d)
          p_scaled[d] = (p[d] - normalization_shift[d]) * normalization_length;

        for (unsigned int d = 0; d < spacedim; ++d)
          result += coefficients[1 + d] * p_scaled[d];

        if (!is_affine)
          {
            Point<dim, Number> result_affine = result;
            for (unsigned int d = 0, c = 0; d < spacedim; ++d)
              for (unsigned int e = 0; e <= d; ++e, ++c)
                result +=
                  coefficients[1 + spacedim + c] * (p_scaled[d] * p_scaled[e]);

            // Check if the quadratic approximation ends up considerably
            // farther outside the unit cell on some or all SIMD lanes than
            // the affine approximation - in that case, we switch those
            // components back to the affine approximation. Note that the
            // quadratic approximation will grow more quickly away from the
            // unit cell. We make the selection for each SIMD lane with a
            // ternary operation.
            const Number distance_to_unit_cell = result.distance_square(
              GeometryInfo<dim>::project_to_unit_cell(result));
            const Number affine_distance_to_unit_cell =
              result_affine.distance_square(
                GeometryInfo<dim>::project_to_unit_cell(result_affine));
            for (unsigned int d = 0; d < dim; ++d)
              result[d] = compare_and_apply_mask<SIMDComparison::greater_than>(
                distance_to_unit_cell,
                affine_distance_to_unit_cell + 0.5,
                result_affine[d],
                result[d]);
          }
        return result;
      }

    private:
      /**
       * In order to guarantee a good conditioning, we need to apply a
       * transformation to the points in real space that is computed by a
       * shift vector normalization_shift (first point of the mapping support
       * points in real space) and an inverse length scale called
       * `length_normalization` as the distance between the first two points.
       */
      const Point<spacedim> normalization_shift;

      /**
       * See the documentation of `normalization_shift` above.
       */
      const double normalization_length;

      /**
       * The vector of coefficients in the quadratic approximation.
       */
      std::array<Point<dim>, n_functions> coefficients;

      /**
       * In case the quadratic approximation is not possible due to an
       * insufficient number of support points, we switch to an affine
       * approximation that always works but is less accurate.
       */
      bool is_affine;
    };



    /**
     * In case the quadrature formula is a tensor product, this is a
     * replacement for maybe_compute_q_points(), maybe_update_Jacobians() and
     * maybe_update_jacobian_grads()
     */
    template <int dim, int spacedim>
    inline void
    maybe_update_q_points_Jacobians_and_grads_tensor(
      const CellSimilarity::Similarity cell_similarity,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      std::vector<Point<spacedim>> &                 quadrature_points,
      std::vector<DerivativeForm<2, dim, spacedim>> &jacobian_grads)
    {
      const UpdateFlags update_flags = data.update_each;

      using VectorizedArrayType =
        typename dealii::MappingQ<dim,
                                  spacedim>::InternalData::VectorizedArrayType;
      const unsigned int     n_shape_values = data.n_shape_functions;
      const unsigned int     n_q_points     = data.shape_info.n_q_points;
      constexpr unsigned int n_lanes        = VectorizedArrayType::size();
      constexpr unsigned int n_comp         = 1 + (spacedim - 1) / n_lanes;
      constexpr unsigned int n_hessians     = (dim * (dim + 1)) / 2;

      EvaluationFlags::EvaluationFlags evaluation_flag =
        (update_flags & update_quadrature_points ? EvaluationFlags::values :
                                                   EvaluationFlags::nothing) |
        ((cell_similarity != CellSimilarity::translation) &&
             (update_flags & update_contravariant_transformation) ?
           EvaluationFlags::gradients :
           EvaluationFlags::nothing) |
        ((cell_similarity != CellSimilarity::translation) &&
             (update_flags & update_jacobian_grads) ?
           EvaluationFlags::hessians :
           EvaluationFlags::nothing);

      Assert(!(evaluation_flag & EvaluationFlags::values) || n_q_points > 0,
             ExcInternalError());
      Assert(!(evaluation_flag & EvaluationFlags::values) ||
               n_q_points == quadrature_points.size(),
             ExcDimensionMismatch(n_q_points, quadrature_points.size()));
      Assert(!(evaluation_flag & EvaluationFlags::gradients) ||
               data.n_shape_functions > 0,
             ExcInternalError());
      Assert(!(evaluation_flag & EvaluationFlags::gradients) ||
               n_q_points == data.contravariant.size(),
             ExcDimensionMismatch(n_q_points, data.contravariant.size()));
      Assert(!(evaluation_flag & EvaluationFlags::hessians) ||
               n_q_points == jacobian_grads.size(),
             ExcDimensionMismatch(n_q_points, jacobian_grads.size()));

      // shortcut in case we have an identity interpolation and only request
      // the quadrature points
      if (evaluation_flag == EvaluationFlags::values &&
          data.shape_info.element_type ==
            internal::MatrixFreeFunctions::tensor_symmetric_collocation)
        {
          for (unsigned int q = 0; q < n_q_points; ++q)
            quadrature_points[q] =
              data.mapping_support_points[data.shape_info
                                            .lexicographic_numbering[q]];
          return;
        }

      FEEvaluationData<dim, VectorizedArrayType, false> eval(data.shape_info);

      // prepare arrays
      if (evaluation_flag != EvaluationFlags::nothing)
        {
          eval.set_data_pointers(&data.scratch, n_comp);

          // make sure to initialize on all lanes also when some are unused in
          // the code below
          for (unsigned int i = 0; i < n_shape_values * n_comp; ++i)
            eval.begin_dof_values()[i] = VectorizedArrayType();

          const std::vector<unsigned int> &renumber_to_lexicographic =
            data.shape_info.lexicographic_numbering;
          for (unsigned int i = 0; i < n_shape_values; ++i)
            for (unsigned int d = 0; d < spacedim; ++d)
              {
                const unsigned int in_comp  = d % n_lanes;
                const unsigned int out_comp = d / n_lanes;
                eval
                  .begin_dof_values()[out_comp * n_shape_values + i][in_comp] =
                  data.mapping_support_points[renumber_to_lexicographic[i]][d];
              }

          // do the actual tensorized evaluation
          internal::FEEvaluationFactory<dim, VectorizedArrayType>::evaluate(
            n_comp, evaluation_flag, eval.begin_dof_values(), eval);
        }

      // do the postprocessing
      if (evaluation_flag & EvaluationFlags::values)
        {
          for (unsigned int out_comp = 0; out_comp < n_comp; ++out_comp)
            for (unsigned int i = 0; i < n_q_points; ++i)
              for (unsigned int in_comp = 0;
                   in_comp < n_lanes && in_comp < spacedim - out_comp * n_lanes;
                   ++in_comp)
                quadrature_points[i][out_comp * n_lanes + in_comp] =
                  eval.begin_values()[out_comp * n_q_points + i][in_comp];
        }

      if (evaluation_flag & EvaluationFlags::gradients)
        {
          std::fill(data.contravariant.begin(),
                    data.contravariant.end(),
                    DerivativeForm<1, dim, spacedim>());
          // We need to reinterpret the data after evaluate has been applied.
          for (unsigned int out_comp = 0; out_comp < n_comp; ++out_comp)
            for (unsigned int point = 0; point < n_q_points; ++point)
              for (unsigned int j = 0; j < dim; ++j)
                for (unsigned int in_comp = 0;
                     in_comp < n_lanes &&
                     in_comp < spacedim - out_comp * n_lanes;
                     ++in_comp)
                  {
                    const unsigned int total_number = point * dim + j;
                    const unsigned int new_comp     = total_number / n_q_points;
                    const unsigned int new_point    = total_number % n_q_points;
                    data.contravariant[new_point][out_comp * n_lanes + in_comp]
                                      [new_comp] =
                      eval.begin_gradients()[(out_comp * n_q_points + point) *
                                               dim +
                                             j][in_comp];
                  }
        }
      if (update_flags & update_covariant_transformation)
        if (cell_similarity != CellSimilarity::translation)
          for (unsigned int point = 0; point < n_q_points; ++point)
            data.covariant[point] =
              (data.contravariant[point]).covariant_form();

      if (update_flags & update_volume_elements)
        if (cell_similarity != CellSimilarity::translation)
          for (unsigned int point = 0; point < n_q_points; ++point)
            data.volume_elements[point] =
              data.contravariant[point].determinant();

      if (evaluation_flag & EvaluationFlags::hessians)
        {
          constexpr int desymmetrize_3d[6][2] = {
            {0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
          constexpr int desymmetrize_2d[3][2] = {{0, 0}, {1, 1}, {0, 1}};

          // We need to reinterpret the data after evaluate has been applied.
          for (unsigned int out_comp = 0; out_comp < n_comp; ++out_comp)
            for (unsigned int point = 0; point < n_q_points; ++point)
              for (unsigned int j = 0; j < n_hessians; ++j)
                for (unsigned int in_comp = 0;
                     in_comp < n_lanes &&
                     in_comp < spacedim - out_comp * n_lanes;
                     ++in_comp)
                  {
                    const unsigned int total_number = point * n_hessians + j;
                    const unsigned int new_point    = total_number % n_q_points;
                    const unsigned int new_hessian_comp =
                      total_number / n_q_points;
                    const unsigned int new_hessian_comp_i =
                      dim == 2 ? desymmetrize_2d[new_hessian_comp][0] :
                                 desymmetrize_3d[new_hessian_comp][0];
                    const unsigned int new_hessian_comp_j =
                      dim == 2 ? desymmetrize_2d[new_hessian_comp][1] :
                                 desymmetrize_3d[new_hessian_comp][1];
                    const double value =
                      eval.begin_hessians()[(out_comp * n_q_points + point) *
                                              n_hessians +
                                            j][in_comp];
                    jacobian_grads[new_point][out_comp * n_lanes + in_comp]
                                  [new_hessian_comp_i][new_hessian_comp_j] =
                                    value;
                    jacobian_grads[new_point][out_comp * n_lanes + in_comp]
                                  [new_hessian_comp_j][new_hessian_comp_i] =
                                    value;
                  }
        }
    }


    /**
     * Compute the locations of quadrature points on the object described by
     * the first argument (and the cell for which the mapping support points
     * have already been set), but only if the update_flags of the @p data
     * argument indicate so.
     */
    template <int dim, int spacedim>
    inline void
    maybe_compute_q_points(
      const typename QProjector<dim>::DataSetDescriptor             data_set,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      std::vector<Point<spacedim>> &quadrature_points)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_quadrature_points)
        for (unsigned int point = 0; point < quadrature_points.size(); ++point)
          {
            const double *  shape = &data.shape(point + data_set, 0);
            Point<spacedim> result =
              (shape[0] * data.mapping_support_points[0]);
            for (unsigned int k = 1; k < data.n_shape_functions; ++k)
              for (unsigned int i = 0; i < spacedim; ++i)
                result[i] += shape[k] * data.mapping_support_points[k][i];
            quadrature_points[point] = result;
          }
    }



    /**
     * Update the co- and contravariant matrices as well as their determinant,
     * for the cell
     * described stored in the data object, but only if the update_flags of the @p data
     * argument indicate so.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    inline void
    maybe_update_Jacobians(
      const CellSimilarity::Similarity                          cell_similarity,
      const typename dealii::QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_contravariant_transformation)
        // if the current cell is just a
        // translation of the previous one, no
        // need to recompute jacobians...
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();

            std::fill(data.contravariant.begin(),
                      data.contravariant.end(),
                      DerivativeForm<1, dim, spacedim>());

            Assert(data.n_shape_functions > 0, ExcInternalError());

            for (unsigned int point = 0; point < n_q_points; ++point)
              {
                double result[spacedim][dim];

                // peel away part of sum to avoid zeroing the
                // entries and adding for the first time
                for (unsigned int i = 0; i < spacedim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    result[i][j] = data.derivative(point + data_set, 0)[j] *
                                   data.mapping_support_points[0][i];
                for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      result[i][j] += data.derivative(point + data_set, k)[j] *
                                      data.mapping_support_points[k][i];

                // write result into contravariant data. for
                // j=dim in the case dim<spacedim, there will
                // never be any nonzero data that arrives in
                // here, so it is ok anyway because it was
                // initialized to zero at the initialization
                for (unsigned int i = 0; i < spacedim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    data.contravariant[point][i][j] = result[i][j];
              }
          }

      if (update_flags & update_covariant_transformation)
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();
            for (unsigned int point = 0; point < n_q_points; ++point)
              {
                data.covariant[point] =
                  (data.contravariant[point]).covariant_form();
              }
          }

      if (update_flags & update_volume_elements)
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();
            for (unsigned int point = 0; point < n_q_points; ++point)
              data.volume_elements[point] =
                data.contravariant[point].determinant();
          }
    }



    /**
     * Update the Hessian of the transformation from unit to real cell, the
     * Jacobian gradients.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    inline void
    maybe_update_jacobian_grads(
      const CellSimilarity::Similarity                  cell_similarity,
      const typename QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      std::vector<DerivativeForm<2, dim, spacedim>> &jacobian_grads)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_grads)
        {
          const unsigned int n_q_points = jacobian_grads.size();

          if (cell_similarity != CellSimilarity::translation)
            for (unsigned int point = 0; point < n_q_points; ++point)
              {
                const Tensor<2, dim> *second =
                  &data.second_derivative(point + data_set, 0);
                double result[spacedim][dim][dim];
                for (unsigned int i = 0; i < spacedim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    for (unsigned int l = 0; l < dim; ++l)
                      result[i][j][l] =
                        (second[0][j][l] * data.mapping_support_points[0][i]);
                for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        result[i][j][l] +=
                          (second[k][j][l] * data.mapping_support_points[k][i]);

                for (unsigned int i = 0; i < spacedim; ++i)
                  for (unsigned int j = 0; j < dim; ++j)
                    for (unsigned int l = 0; l < dim; ++l)
                      jacobian_grads[point][i][j][l] = result[i][j][l];
              }
        }
    }



    /**
     * Update the Hessian of the transformation from unit to real cell, the
     * Jacobian gradients, pushed forward to the real cell coordinates.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    inline void
    maybe_update_jacobian_pushed_forward_grads(
      const CellSimilarity::Similarity                  cell_similarity,
      const typename QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      std::vector<Tensor<3, spacedim>> &jacobian_pushed_forward_grads)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_pushed_forward_grads)
        {
          const unsigned int n_q_points = jacobian_pushed_forward_grads.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              double tmp[spacedim][spacedim][spacedim];
              for (unsigned int point = 0; point < n_q_points; ++point)
                {
                  const Tensor<2, dim> *second =
                    &data.second_derivative(point + data_set, 0);
                  double result[spacedim][dim][dim];
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        result[i][j][l] =
                          (second[0][j][l] * data.mapping_support_points[0][i]);
                  for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          result[i][j][l] +=
                            (second[k][j][l] *
                             data.mapping_support_points[k][i]);

                  // first push forward the j-components
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        {
                          tmp[i][j][l] =
                            result[i][0][l] * data.covariant[point][j][0];
                          for (unsigned int jr = 1; jr < dim; ++jr)
                            {
                              tmp[i][j][l] +=
                                result[i][jr][l] * data.covariant[point][j][jr];
                            }
                        }

                  // now, pushing forward the l-components
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < spacedim; ++l)
                        {
                          jacobian_pushed_forward_grads[point][i][j][l] =
                            tmp[i][j][0] * data.covariant[point][l][0];
                          for (unsigned int lr = 1; lr < dim; ++lr)
                            {
                              jacobian_pushed_forward_grads[point][i][j][l] +=
                                tmp[i][j][lr] * data.covariant[point][l][lr];
                            }
                        }
                }
            }
        }
    }



    /**
     * Update the third derivatives of the transformation from unit to real
     * cell, the Jacobian hessians.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    inline void
    maybe_update_jacobian_2nd_derivatives(
      const CellSimilarity::Similarity                  cell_similarity,
      const typename QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      std::vector<DerivativeForm<3, dim, spacedim>> &jacobian_2nd_derivatives)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_2nd_derivatives)
        {
          const unsigned int n_q_points = jacobian_2nd_derivatives.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              for (unsigned int point = 0; point < n_q_points; ++point)
                {
                  const Tensor<3, dim> *third =
                    &data.third_derivative(point + data_set, 0);
                  double result[spacedim][dim][dim][dim];
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          result[i][j][l][m] =
                            (third[0][j][l][m] *
                             data.mapping_support_points[0][i]);
                  for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            result[i][j][l][m] +=
                              (third[k][j][l][m] *
                               data.mapping_support_points[k][i]);

                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          jacobian_2nd_derivatives[point][i][j][l][m] =
                            result[i][j][l][m];
                }
            }
        }
    }



    /**
     * Update the Hessian of the Hessian of the transformation from unit
     * to real cell, the Jacobian Hessian gradients, pushed forward to the
     * real cell coordinates.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    inline void
    maybe_update_jacobian_pushed_forward_2nd_derivatives(
      const CellSimilarity::Similarity                  cell_similarity,
      const typename QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      std::vector<Tensor<4, spacedim>> &jacobian_pushed_forward_2nd_derivatives)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_pushed_forward_2nd_derivatives)
        {
          const unsigned int n_q_points =
            jacobian_pushed_forward_2nd_derivatives.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              double tmp[spacedim][spacedim][spacedim][spacedim];
              for (unsigned int point = 0; point < n_q_points; ++point)
                {
                  const Tensor<3, dim> *third =
                    &data.third_derivative(point + data_set, 0);
                  double result[spacedim][dim][dim][dim];
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          result[i][j][l][m] =
                            (third[0][j][l][m] *
                             data.mapping_support_points[0][i]);
                  for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            result[i][j][l][m] +=
                              (third[k][j][l][m] *
                               data.mapping_support_points[k][i]);

                  // push forward the j-coordinate
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          {
                            jacobian_pushed_forward_2nd_derivatives
                              [point][i][j][l][m] = result[i][0][l][m] *
                                                    data.covariant[point][j][0];
                            for (unsigned int jr = 1; jr < dim; ++jr)
                              jacobian_pushed_forward_2nd_derivatives[point][i]
                                                                     [j][l]
                                                                     [m] +=
                                result[i][jr][l][m] *
                                data.covariant[point][j][jr];
                          }

                  // push forward the l-coordinate
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < spacedim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          {
                            tmp[i][j][l][m] =
                              jacobian_pushed_forward_2nd_derivatives[point][i]
                                                                     [j][0][m] *
                              data.covariant[point][l][0];
                            for (unsigned int lr = 1; lr < dim; ++lr)
                              tmp[i][j][l][m] +=
                                jacobian_pushed_forward_2nd_derivatives[point]
                                                                       [i][j]
                                                                       [lr][m] *
                                data.covariant[point][l][lr];
                          }

                  // push forward the m-coordinate
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < spacedim; ++l)
                        for (unsigned int m = 0; m < spacedim; ++m)
                          {
                            jacobian_pushed_forward_2nd_derivatives
                              [point][i][j][l][m] =
                                tmp[i][j][l][0] * data.covariant[point][m][0];
                            for (unsigned int mr = 1; mr < dim; ++mr)
                              jacobian_pushed_forward_2nd_derivatives[point][i]
                                                                     [j][l]
                                                                     [m] +=
                                tmp[i][j][l][mr] * data.covariant[point][m][mr];
                          }
                }
            }
        }
    }



    /**
     * Update the fourth derivatives of the transformation from unit to real
     * cell, the Jacobian hessian gradients.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    inline void
    maybe_update_jacobian_3rd_derivatives(
      const CellSimilarity::Similarity                  cell_similarity,
      const typename QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      std::vector<DerivativeForm<4, dim, spacedim>> &jacobian_3rd_derivatives)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_3rd_derivatives)
        {
          const unsigned int n_q_points = jacobian_3rd_derivatives.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              for (unsigned int point = 0; point < n_q_points; ++point)
                {
                  const Tensor<4, dim> *fourth =
                    &data.fourth_derivative(point + data_set, 0);
                  double result[spacedim][dim][dim][dim][dim];
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          for (unsigned int n = 0; n < dim; ++n)
                            result[i][j][l][m][n] =
                              (fourth[0][j][l][m][n] *
                               data.mapping_support_points[0][i]);
                  for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            for (unsigned int n = 0; n < dim; ++n)
                              result[i][j][l][m][n] +=
                                (fourth[k][j][l][m][n] *
                                 data.mapping_support_points[k][i]);

                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          for (unsigned int n = 0; n < dim; ++n)
                            jacobian_3rd_derivatives[point][i][j][l][m][n] =
                              result[i][j][l][m][n];
                }
            }
        }
    }



    /**
     * Update the Hessian gradient of the transformation from unit to real
     * cell, the Jacobian Hessians, pushed forward to the real cell
     * coordinates.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    inline void
    maybe_update_jacobian_pushed_forward_3rd_derivatives(
      const CellSimilarity::Similarity                  cell_similarity,
      const typename QProjector<dim>::DataSetDescriptor data_set,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      std::vector<Tensor<5, spacedim>> &jacobian_pushed_forward_3rd_derivatives)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_pushed_forward_3rd_derivatives)
        {
          const unsigned int n_q_points =
            jacobian_pushed_forward_3rd_derivatives.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              double tmp[spacedim][spacedim][spacedim][spacedim][spacedim];
              for (unsigned int point = 0; point < n_q_points; ++point)
                {
                  const Tensor<4, dim> *fourth =
                    &data.fourth_derivative(point + data_set, 0);
                  double result[spacedim][dim][dim][dim][dim];
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < dim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          for (unsigned int n = 0; n < dim; ++n)
                            result[i][j][l][m][n] =
                              (fourth[0][j][l][m][n] *
                               data.mapping_support_points[0][i]);
                  for (unsigned int k = 1; k < data.n_shape_functions; ++k)
                    for (unsigned int i = 0; i < spacedim; ++i)
                      for (unsigned int j = 0; j < dim; ++j)
                        for (unsigned int l = 0; l < dim; ++l)
                          for (unsigned int m = 0; m < dim; ++m)
                            for (unsigned int n = 0; n < dim; ++n)
                              result[i][j][l][m][n] +=
                                (fourth[k][j][l][m][n] *
                                 data.mapping_support_points[k][i]);

                  // push-forward the j-coordinate
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < dim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          for (unsigned int n = 0; n < dim; ++n)
                            {
                              tmp[i][j][l][m][n] = result[i][0][l][m][n] *
                                                   data.covariant[point][j][0];
                              for (unsigned int jr = 1; jr < dim; ++jr)
                                tmp[i][j][l][m][n] +=
                                  result[i][jr][l][m][n] *
                                  data.covariant[point][j][jr];
                            }

                  // push-forward the l-coordinate
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < spacedim; ++l)
                        for (unsigned int m = 0; m < dim; ++m)
                          for (unsigned int n = 0; n < dim; ++n)
                            {
                              jacobian_pushed_forward_3rd_derivatives
                                [point][i][j][l][m][n] =
                                  tmp[i][j][0][m][n] *
                                  data.covariant[point][l][0];
                              for (unsigned int lr = 1; lr < dim; ++lr)
                                jacobian_pushed_forward_3rd_derivatives[point]
                                                                       [i][j][l]
                                                                       [m][n] +=
                                  tmp[i][j][lr][m][n] *
                                  data.covariant[point][l][lr];
                            }

                  // push-forward the m-coordinate
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < spacedim; ++l)
                        for (unsigned int m = 0; m < spacedim; ++m)
                          for (unsigned int n = 0; n < dim; ++n)
                            {
                              tmp[i][j][l][m][n] =
                                jacobian_pushed_forward_3rd_derivatives[point]
                                                                       [i][j][l]
                                                                       [0][n] *
                                data.covariant[point][m][0];
                              for (unsigned int mr = 1; mr < dim; ++mr)
                                tmp[i][j][l][m][n] +=
                                  jacobian_pushed_forward_3rd_derivatives
                                    [point][i][j][l][mr][n] *
                                  data.covariant[point][m][mr];
                            }

                  // push-forward the n-coordinate
                  for (unsigned int i = 0; i < spacedim; ++i)
                    for (unsigned int j = 0; j < spacedim; ++j)
                      for (unsigned int l = 0; l < spacedim; ++l)
                        for (unsigned int m = 0; m < spacedim; ++m)
                          for (unsigned int n = 0; n < spacedim; ++n)
                            {
                              jacobian_pushed_forward_3rd_derivatives
                                [point][i][j][l][m][n] =
                                  tmp[i][j][l][m][0] *
                                  data.covariant[point][n][0];
                              for (unsigned int nr = 1; nr < dim; ++nr)
                                jacobian_pushed_forward_3rd_derivatives[point]
                                                                       [i][j][l]
                                                                       [m][n] +=
                                  tmp[i][j][l][m][nr] *
                                  data.covariant[point][n][nr];
                            }
                }
            }
        }
    }



    /**
     * Depending on what information is called for in the update flags of the
     * @p data object, compute the various pieces of information that is
     * required by the fill_fe_face_values() and fill_fe_subface_values()
     * functions. This function simply unifies the work that would be done by
     * those two functions.
     *
     * The resulting data is put into the @p output_data argument.
     */
    template <int dim, int spacedim>
    inline void
    maybe_compute_face_data(
      const dealii::MappingQ<dim, spacedim> &mapping,
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
      const unsigned int                                            face_no,
      const unsigned int                                            subface_no,
      const unsigned int                                            n_q_points,
      const std::vector<double> &                                   weights,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
        &output_data)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags &
          (update_boundary_forms | update_normal_vectors | update_jacobians |
           update_JxW_values | update_inverse_jacobians))
        {
          if (update_flags & update_boundary_forms)
            AssertDimension(output_data.boundary_forms.size(), n_q_points);
          if (update_flags & update_normal_vectors)
            AssertDimension(output_data.normal_vectors.size(), n_q_points);
          if (update_flags & update_JxW_values)
            AssertDimension(output_data.JxW_values.size(), n_q_points);

          Assert(data.aux.size() + 1 >= dim, ExcInternalError());

          // first compute some common data that is used for evaluating
          // all of the flags below

          // map the unit tangentials to the real cell. checking for d!=dim-1
          // eliminates compiler warnings regarding unsigned int expressions <
          // 0.
          for (unsigned int d = 0; d != dim - 1; ++d)
            {
              Assert(face_no + GeometryInfo<dim>::faces_per_cell * d <
                       data.unit_tangentials.size(),
                     ExcInternalError());
              Assert(
                data.aux[d].size() <=
                  data
                    .unit_tangentials[face_no +
                                      GeometryInfo<dim>::faces_per_cell * d]
                    .size(),
                ExcInternalError());

              mapping.transform(
                make_array_view(
                  data.unit_tangentials[face_no +
                                        GeometryInfo<dim>::faces_per_cell * d]),
                mapping_contravariant,
                data,
                make_array_view(data.aux[d].begin(), data.aux[d].end()));
            }

          if (update_flags & update_boundary_forms)
            {
              // if dim==spacedim, we can use the unit tangentials to compute
              // the boundary form by simply taking the cross product
              if (dim == spacedim)
                {
                  for (unsigned int i = 0; i < n_q_points; ++i)
                    switch (dim)
                      {
                        case 1:
                          // in 1d, we don't have access to any of the
                          // data.aux fields (because it has only dim-1
                          // components), but we can still compute the
                          // boundary form by simply looking at the number of
                          // the face
                          output_data.boundary_forms[i][0] =
                            (face_no == 0 ? -1 : +1);
                          break;
                        case 2:
                          output_data.boundary_forms[i] =
                            cross_product_2d(data.aux[0][i]);
                          break;
                        case 3:
                          output_data.boundary_forms[i] =
                            cross_product_3d(data.aux[0][i], data.aux[1][i]);
                          break;
                        default:
                          Assert(false, ExcNotImplemented());
                      }
                }
              else //(dim < spacedim)
                {
                  // in the codim-one case, the boundary form results from the
                  // cross product of all the face tangential vectors and the
                  // cell normal vector
                  //
                  // to compute the cell normal, use the same method used in
                  // fill_fe_values for cells above
                  AssertDimension(data.contravariant.size(), n_q_points);

                  for (unsigned int point = 0; point < n_q_points; ++point)
                    {
                      if (dim == 1)
                        {
                          // J is a tangent vector
                          output_data.boundary_forms[point] =
                            data.contravariant[point].transpose()[0];
                          output_data.boundary_forms[point] /=
                            (face_no == 0 ? -1. : +1.) *
                            output_data.boundary_forms[point].norm();
                        }

                      if (dim == 2)
                        {
                          const DerivativeForm<1, spacedim, dim> DX_t =
                            data.contravariant[point].transpose();

                          Tensor<1, spacedim> cell_normal =
                            cross_product_3d(DX_t[0], DX_t[1]);
                          cell_normal /= cell_normal.norm();

                          // then compute the face normal from the face
                          // tangent and the cell normal:
                          output_data.boundary_forms[point] =
                            cross_product_3d(data.aux[0][point], cell_normal);
                        }
                    }
                }
            }

          if (update_flags & update_JxW_values)
            for (unsigned int i = 0; i < output_data.boundary_forms.size(); ++i)
              {
                output_data.JxW_values[i] =
                  output_data.boundary_forms[i].norm() * weights[i];

                if (subface_no != numbers::invalid_unsigned_int)
                  {
                    const double area_ratio = GeometryInfo<dim>::subface_ratio(
                      cell->subface_case(face_no), subface_no);
                    output_data.JxW_values[i] *= area_ratio;
                  }
              }

          if (update_flags & update_normal_vectors)
            for (unsigned int i = 0; i < output_data.normal_vectors.size(); ++i)
              output_data.normal_vectors[i] =
                Point<spacedim>(output_data.boundary_forms[i] /
                                output_data.boundary_forms[i].norm());

          if (update_flags & update_jacobians)
            for (unsigned int point = 0; point < n_q_points; ++point)
              output_data.jacobians[point] = data.contravariant[point];

          if (update_flags & update_inverse_jacobians)
            for (unsigned int point = 0; point < n_q_points; ++point)
              output_data.inverse_jacobians[point] =
                data.covariant[point].transpose();
        }
    }


    /**
     * Do the work of MappingQ::fill_fe_face_values() and
     * MappingQ::fill_fe_subface_values() in a generic way,
     * using the 'data_set' to differentiate whether we will
     * work on a face (and if so, which one) or subface.
     */
    template <int dim, int spacedim>
    inline void
    do_fill_fe_face_values(
      const dealii::MappingQ<dim, spacedim> &mapping,
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
      const unsigned int                                            face_no,
      const unsigned int                                            subface_no,
      const typename QProjector<dim>::DataSetDescriptor             data_set,
      const Quadrature<dim - 1> &                                   quadrature,
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data,
      internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
        &output_data)
    {
      if (dim > 1 && data.tensor_product_quadrature)
        {
          maybe_update_q_points_Jacobians_and_grads_tensor<dim, spacedim>(
            CellSimilarity::none,
            data,
            output_data.quadrature_points,
            output_data.jacobian_grads);
        }
      else
        {
          maybe_compute_q_points<dim, spacedim>(data_set,
                                                data,
                                                output_data.quadrature_points);
          maybe_update_Jacobians<dim, spacedim>(CellSimilarity::none,
                                                data_set,
                                                data);
          maybe_update_jacobian_grads<dim, spacedim>(
            CellSimilarity::none, data_set, data, output_data.jacobian_grads);
        }
      maybe_update_jacobian_pushed_forward_grads<dim, spacedim>(
        CellSimilarity::none,
        data_set,
        data,
        output_data.jacobian_pushed_forward_grads);
      maybe_update_jacobian_2nd_derivatives<dim, spacedim>(
        CellSimilarity::none,
        data_set,
        data,
        output_data.jacobian_2nd_derivatives);
      maybe_update_jacobian_pushed_forward_2nd_derivatives<dim, spacedim>(
        CellSimilarity::none,
        data_set,
        data,
        output_data.jacobian_pushed_forward_2nd_derivatives);
      maybe_update_jacobian_3rd_derivatives<dim, spacedim>(
        CellSimilarity::none,
        data_set,
        data,
        output_data.jacobian_3rd_derivatives);
      maybe_update_jacobian_pushed_forward_3rd_derivatives<dim, spacedim>(
        CellSimilarity::none,
        data_set,
        data,
        output_data.jacobian_pushed_forward_3rd_derivatives);

      maybe_compute_face_data(mapping,
                              cell,
                              face_no,
                              subface_no,
                              quadrature.size(),
                              quadrature.get_weights(),
                              data,
                              output_data);
    }



    /**
     * Implementation of MappingQ::transform() for generic tensors.
     */
    template <int dim, int spacedim, int rank>
    inline void
    transform_fields(
      const ArrayView<const Tensor<rank, dim>> &               input,
      const MappingKind                                        mapping_kind,
      const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
      const ArrayView<Tensor<rank, spacedim>> &                output)
    {
      AssertDimension(input.size(), output.size());
      Assert((dynamic_cast<
                const typename dealii::MappingQ<dim, spacedim>::InternalData *>(
                &mapping_data) != nullptr),
             ExcInternalError());
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data =
        static_cast<
          const typename dealii::MappingQ<dim, spacedim>::InternalData &>(
          mapping_data);

      switch (mapping_kind)
        {
          case mapping_contravariant:
            {
              Assert(data.update_each & update_contravariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_contravariant_transformation"));

              for (unsigned int i = 0; i < output.size(); ++i)
                output[i] =
                  apply_transformation(data.contravariant[i], input[i]);

              return;
            }

          case mapping_piola:
            {
              Assert(data.update_each & update_contravariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_contravariant_transformation"));
              Assert(data.update_each & update_volume_elements,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_volume_elements"));
              Assert(rank == 1, ExcMessage("Only for rank 1"));
              if (rank != 1)
                return;

              for (unsigned int i = 0; i < output.size(); ++i)
                {
                  output[i] =
                    apply_transformation(data.contravariant[i], input[i]);
                  output[i] /= data.volume_elements[i];
                }
              return;
            }
          // We still allow this operation as in the
          // reference cell Derivatives are Tensor
          // rather than DerivativeForm
          case mapping_covariant:
            {
              Assert(data.update_each & update_contravariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_covariant_transformation"));

              for (unsigned int i = 0; i < output.size(); ++i)
                output[i] = apply_transformation(data.covariant[i], input[i]);

              return;
            }

          default:
            Assert(false, ExcNotImplemented());
        }
    }



    /**
     * Implementation of MappingQ::transform() for gradients.
     */
    template <int dim, int spacedim, int rank>
    inline void
    transform_gradients(
      const ArrayView<const Tensor<rank, dim>> &               input,
      const MappingKind                                        mapping_kind,
      const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
      const ArrayView<Tensor<rank, spacedim>> &                output)
    {
      AssertDimension(input.size(), output.size());
      Assert((dynamic_cast<
                const typename dealii::MappingQ<dim, spacedim>::InternalData *>(
                &mapping_data) != nullptr),
             ExcInternalError());
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data =
        static_cast<
          const typename dealii::MappingQ<dim, spacedim>::InternalData &>(
          mapping_data);

      switch (mapping_kind)
        {
          case mapping_contravariant_gradient:
            {
              Assert(data.update_each & update_covariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_covariant_transformation"));
              Assert(data.update_each & update_contravariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_contravariant_transformation"));
              Assert(rank == 2, ExcMessage("Only for rank 2"));

              for (unsigned int i = 0; i < output.size(); ++i)
                {
                  const DerivativeForm<1, spacedim, dim> A =
                    apply_transformation(data.contravariant[i],
                                         transpose(input[i]));
                  output[i] =
                    apply_transformation(data.covariant[i], A.transpose());
                }

              return;
            }

          case mapping_covariant_gradient:
            {
              Assert(data.update_each & update_covariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_covariant_transformation"));
              Assert(rank == 2, ExcMessage("Only for rank 2"));

              for (unsigned int i = 0; i < output.size(); ++i)
                {
                  const DerivativeForm<1, spacedim, dim> A =
                    apply_transformation(data.covariant[i],
                                         transpose(input[i]));
                  output[i] =
                    apply_transformation(data.covariant[i], A.transpose());
                }

              return;
            }

          case mapping_piola_gradient:
            {
              Assert(data.update_each & update_covariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_covariant_transformation"));
              Assert(data.update_each & update_contravariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_contravariant_transformation"));
              Assert(data.update_each & update_volume_elements,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_volume_elements"));
              Assert(rank == 2, ExcMessage("Only for rank 2"));

              for (unsigned int i = 0; i < output.size(); ++i)
                {
                  const DerivativeForm<1, spacedim, dim> A =
                    apply_transformation(data.covariant[i], input[i]);
                  const Tensor<2, spacedim> T =
                    apply_transformation(data.contravariant[i], A.transpose());

                  output[i] = transpose(T);
                  output[i] /= data.volume_elements[i];
                }

              return;
            }

          default:
            Assert(false, ExcNotImplemented());
        }
    }



    /**
     * Implementation of MappingQ::transform() for hessians.
     */
    template <int dim, int spacedim>
    inline void
    transform_hessians(
      const ArrayView<const Tensor<3, dim>> &                  input,
      const MappingKind                                        mapping_kind,
      const typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
      const ArrayView<Tensor<3, spacedim>> &                   output)
    {
      AssertDimension(input.size(), output.size());
      Assert((dynamic_cast<
                const typename dealii::MappingQ<dim, spacedim>::InternalData *>(
                &mapping_data) != nullptr),
             ExcInternalError());
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data =
        static_cast<
          const typename dealii::MappingQ<dim, spacedim>::InternalData &>(
          mapping_data);

      switch (mapping_kind)
        {
          case mapping_contravariant_hessian:
            {
              Assert(data.update_each & update_covariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_covariant_transformation"));
              Assert(data.update_each & update_contravariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_contravariant_transformation"));

              for (unsigned int q = 0; q < output.size(); ++q)
                for (unsigned int i = 0; i < spacedim; ++i)
                  {
                    double tmp1[dim][dim];
                    for (unsigned int J = 0; J < dim; ++J)
                      for (unsigned int K = 0; K < dim; ++K)
                        {
                          tmp1[J][K] =
                            data.contravariant[q][i][0] * input[q][0][J][K];
                          for (unsigned int I = 1; I < dim; ++I)
                            tmp1[J][K] +=
                              data.contravariant[q][i][I] * input[q][I][J][K];
                        }
                    for (unsigned int j = 0; j < spacedim; ++j)
                      {
                        double tmp2[dim];
                        for (unsigned int K = 0; K < dim; ++K)
                          {
                            tmp2[K] = data.covariant[q][j][0] * tmp1[0][K];
                            for (unsigned int J = 1; J < dim; ++J)
                              tmp2[K] += data.covariant[q][j][J] * tmp1[J][K];
                          }
                        for (unsigned int k = 0; k < spacedim; ++k)
                          {
                            output[q][i][j][k] =
                              data.covariant[q][k][0] * tmp2[0];
                            for (unsigned int K = 1; K < dim; ++K)
                              output[q][i][j][k] +=
                                data.covariant[q][k][K] * tmp2[K];
                          }
                      }
                  }
              return;
            }

          case mapping_covariant_hessian:
            {
              Assert(data.update_each & update_covariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_covariant_transformation"));

              for (unsigned int q = 0; q < output.size(); ++q)
                for (unsigned int i = 0; i < spacedim; ++i)
                  {
                    double tmp1[dim][dim];
                    for (unsigned int J = 0; J < dim; ++J)
                      for (unsigned int K = 0; K < dim; ++K)
                        {
                          tmp1[J][K] =
                            data.covariant[q][i][0] * input[q][0][J][K];
                          for (unsigned int I = 1; I < dim; ++I)
                            tmp1[J][K] +=
                              data.covariant[q][i][I] * input[q][I][J][K];
                        }
                    for (unsigned int j = 0; j < spacedim; ++j)
                      {
                        double tmp2[dim];
                        for (unsigned int K = 0; K < dim; ++K)
                          {
                            tmp2[K] = data.covariant[q][j][0] * tmp1[0][K];
                            for (unsigned int J = 1; J < dim; ++J)
                              tmp2[K] += data.covariant[q][j][J] * tmp1[J][K];
                          }
                        for (unsigned int k = 0; k < spacedim; ++k)
                          {
                            output[q][i][j][k] =
                              data.covariant[q][k][0] * tmp2[0];
                            for (unsigned int K = 1; K < dim; ++K)
                              output[q][i][j][k] +=
                                data.covariant[q][k][K] * tmp2[K];
                          }
                      }
                  }

              return;
            }

          case mapping_piola_hessian:
            {
              Assert(data.update_each & update_covariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_covariant_transformation"));
              Assert(data.update_each & update_contravariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_contravariant_transformation"));
              Assert(data.update_each & update_volume_elements,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_volume_elements"));

              for (unsigned int q = 0; q < output.size(); ++q)
                for (unsigned int i = 0; i < spacedim; ++i)
                  {
                    double factor[dim];
                    for (unsigned int I = 0; I < dim; ++I)
                      factor[I] =
                        data.contravariant[q][i][I] / data.volume_elements[q];
                    double tmp1[dim][dim];
                    for (unsigned int J = 0; J < dim; ++J)
                      for (unsigned int K = 0; K < dim; ++K)
                        {
                          tmp1[J][K] = factor[0] * input[q][0][J][K];
                          for (unsigned int I = 1; I < dim; ++I)
                            tmp1[J][K] += factor[I] * input[q][I][J][K];
                        }
                    for (unsigned int j = 0; j < spacedim; ++j)
                      {
                        double tmp2[dim];
                        for (unsigned int K = 0; K < dim; ++K)
                          {
                            tmp2[K] = data.covariant[q][j][0] * tmp1[0][K];
                            for (unsigned int J = 1; J < dim; ++J)
                              tmp2[K] += data.covariant[q][j][J] * tmp1[J][K];
                          }
                        for (unsigned int k = 0; k < spacedim; ++k)
                          {
                            output[q][i][j][k] =
                              data.covariant[q][k][0] * tmp2[0];
                            for (unsigned int K = 1; K < dim; ++K)
                              output[q][i][j][k] +=
                                data.covariant[q][k][K] * tmp2[K];
                          }
                      }
                  }

              return;
            }

          default:
            Assert(false, ExcNotImplemented());
        }
    }



    /**
     * Implementation of MappingQ::transform() for DerivativeForm
     * arguments.
     */
    template <int dim, int spacedim, int rank>
    inline void
    transform_differential_forms(
      const ArrayView<const DerivativeForm<rank, dim, spacedim>> &input,
      const MappingKind                                           mapping_kind,
      const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_data,
      const ArrayView<Tensor<rank + 1, spacedim>> &               output)
    {
      AssertDimension(input.size(), output.size());
      Assert((dynamic_cast<
                const typename dealii::MappingQ<dim, spacedim>::InternalData *>(
                &mapping_data) != nullptr),
             ExcInternalError());
      const typename dealii::MappingQ<dim, spacedim>::InternalData &data =
        static_cast<
          const typename dealii::MappingQ<dim, spacedim>::InternalData &>(
          mapping_data);

      switch (mapping_kind)
        {
          case mapping_covariant:
            {
              Assert(data.update_each & update_contravariant_transformation,
                     typename FEValuesBase<dim>::ExcAccessToUninitializedField(
                       "update_covariant_transformation"));

              for (unsigned int i = 0; i < output.size(); ++i)
                output[i] = apply_transformation(data.covariant[i], input[i]);

              return;
            }
          default:
            Assert(false, ExcNotImplemented());
        }
    }
  } // namespace MappingQImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
