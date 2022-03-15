// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2021 by the deal.II authors
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

#include <deal.II/base/function_tools.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/tools/roots.hpp>

#include <algorithm>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  namespace internal
  {
    namespace QuadratureGeneratorImplementation
    {
      template <int dim>
      void
      tensor_point_with_1D_quadrature(const Point<dim - 1> &point,
                                      const double          weight,
                                      const Quadrature<1> & quadrature1D,
                                      const double          start,
                                      const double          end,
                                      const unsigned int    component_in_dim,
                                      ExtendableQuadrature<dim> &quadrature)
      {
        Assert(start < end,
               ExcMessage("Interval start must be less than interval end."));

        const double length = end - start;
        for (unsigned int j = 0; j < quadrature1D.size(); ++j)
          {
            const double x = start + length * quadrature1D.point(j)[0];
            quadrature.push_back(dealii::internal::create_higher_dim_point(
                                   point, component_in_dim, x),
                                 length * weight * quadrature1D.weight(j));
          }
      }



      /**
       * For each (point, weight) in lower create a dim-dimensional quadrature
       * using tensor_point_with_1D_quadrature and add the results to @p quadrature.
       */
      template <int dim>
      void
      add_tensor_product(const Quadrature<dim - 1> &lower,
                         const Quadrature<1> &      quadrature1D,
                         const double               start,
                         const double               end,
                         const unsigned int         component_in_dim,
                         ExtendableQuadrature<dim> &quadrature)
      {
        for (unsigned int j = 0; j < lower.size(); ++j)
          {
            tensor_point_with_1D_quadrature(lower.point(j),
                                            lower.weight(j),
                                            quadrature1D,
                                            start,
                                            end,
                                            component_in_dim,
                                            quadrature);
          }
      }



      template <int dim>
      Definiteness
      pointwise_definiteness(
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &               functions,
        const Point<dim> &point)
      {
        Assert(functions.size() > 0,
               ExcMessage(
                 "The incoming vector must contain at least one function."));

        const int sign_of_first =
          boost::math::sign(functions[0].get().value(point));

        if (sign_of_first == 0)
          return Definiteness::indefinite;

        for (unsigned int j = 1; j < functions.size(); ++j)
          {
            const int sign = boost::math::sign(functions[j].get().value(point));

            if (sign != sign_of_first)
              return Definiteness::indefinite;
          }
        // If we got here all functions have the same sign.
        if (sign_of_first < 0)
          return Definiteness::negative;
        else
          return Definiteness::positive;
      }



      /**
       * Given the incoming lower and upper bounds on the value of a function
       * $[L, U]$, return the minimum/maximum of $[L, U]$ and the function
       * values at the vertices. That is, this function returns
       *
       * $[\min(L, L_f), \max(U, U_f)]$,
       *
       * where $L_f = \min_{v} f(x_v)$, $U_f = \max_{v} f(x_v)|$,
       * and $x_v$ is a vertex.
       *
       * It is assumed that the incoming function is scalar valued.
       */
      template <int dim>
      void
      take_min_max_at_vertices(const Function<dim> &      function,
                               const BoundingBox<dim> &   box,
                               std::pair<double, double> &value_bounds)
      {
        const ReferenceCell &cube = ReferenceCells::get_hypercube<dim>();
        for (unsigned int i = 0; i < cube.n_vertices(); ++i)
          {
            const double vertex_value = function.value(box.vertex(i));

            value_bounds.first  = std::min(value_bounds.first, vertex_value);
            value_bounds.second = std::max(value_bounds.second, vertex_value);
          }
      }



      /**
       * Estimate bounds on each of the functions in the incoming vector over
       * the incoming box.
       *
       * Bounds on the functions value and the gradient components are first
       * computed using FunctionTools::taylor_estimate_function_bounds.
       * In addition, the function value is checked for min/max at the at
       * the vertices of the box. The gradient is not checked at the box
       * vertices.
       */
      template <int dim>
      void
      estimate_function_bounds(
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                               functions,
        const BoundingBox<dim> &          box,
        std::vector<FunctionBounds<dim>> &all_function_bounds)
      {
        all_function_bounds.clear();
        all_function_bounds.reserve(functions.size());
        for (const Function<dim> &function : functions)
          {
            FunctionBounds<dim> bounds;
            FunctionTools::taylor_estimate_function_bounds<dim>(
              function, box, bounds.value, bounds.gradient);
            take_min_max_at_vertices(function, box, bounds.value);

            all_function_bounds.push_back(bounds);
          }
      }



      template <int dim>
      std::pair<double, double>
      find_extreme_values(const std::vector<FunctionBounds<dim>> &bounds)
      {
        Assert(bounds.size() > 0, ExcMessage("The incoming vector is empty."));

        std::pair<double, double> extremes = bounds[0].value;
        for (unsigned int i = 1; i < bounds.size(); ++i)
          {
            extremes.first  = std::min(extremes.first, bounds[i].value.first);
            extremes.second = std::max(extremes.second, bounds[i].value.second);
          }

        return extremes;
      }



      /**
       * Return true if the incoming function bounds correspond to a function
       * which is indefinite, i.e., that is not negative or positive definite.
       */
      inline bool
      is_indefinite(const std::pair<double, double> &function_bounds)
      {
        if (function_bounds.first > 0)
          return false;
        if (function_bounds.second < 0)
          return false;
        return true;
      }



      /**
       * Return a lower bound, $L_a$, on the absolute value of a function,
       * $f(x)$:
       *
       * $L_a \leq |f(x)|$,
       *
       * by estimating it from the incoming lower and upper bounds:
       * $L \leq f(x) \leq U$.
       *
       * By rewriting the lower and upper bounds as
       * $F - C \leq f(x) \leq F + C$,
       * where $L = F - C$, $U = F + C$ (or $F = (U + L)/2$, $C = (U - L)/2$),
       * we get $|f(x) - F| \leq C$.
       * Using the inverse triangle inequality gives
       * $|F| - |f(x)| \leq |f(x) - F| \leq C$.
       * Thus, $L_a = |F| - C$.
       *
       * Note that the returned value can be negative. This is used to indicate
       * "how far away" a function is from being definite.
       */
      inline double
      lower_bound_on_abs(const std::pair<double, double> &function_bounds)
      {
        Assert(function_bounds.first <= function_bounds.second,
               ExcMessage("Function bounds reversed, max < min."));

        return 0.5 * (std::abs(function_bounds.second + function_bounds.first) -
                      (function_bounds.second - function_bounds.first));
      }



      HeightDirectionData::HeightDirectionData()
      {
        direction    = numbers::invalid_unsigned_int;
        min_abs_dfdx = 0;
      }



      template <int dim>
      std_cxx17::optional<HeightDirectionData>
      find_best_height_direction(
        const std::vector<FunctionBounds<dim>> &all_function_bounds)
      {
        // Minimum (taken over the indefinite functions) on the lower bound on
        // each component of the gradient.
        std_cxx17::optional<std::array<double, dim>> min_lower_abs_grad;

        for (const FunctionBounds<dim> &bounds : all_function_bounds)
          {
            if (is_indefinite(bounds.value))
              {
                // For the first indefinite function we find, we write the lower
                // bounds on each gradient component to min_lower_abs_grad.
                if (!min_lower_abs_grad)
                  {
                    min_lower_abs_grad.emplace();
                    for (unsigned int d = 0; d < dim; ++d)
                      {
                        (*min_lower_abs_grad)[d] =
                          lower_bound_on_abs(bounds.gradient[d]);
                      }
                  }
                else
                  {
                    for (unsigned int d = 0; d < dim; ++d)
                      {
                        (*min_lower_abs_grad)[d] =
                          std::min((*min_lower_abs_grad)[d],
                                   lower_bound_on_abs(bounds.gradient[d]));
                      }
                  }
              }
          }

        if (min_lower_abs_grad)
          {
            const auto max_element =
              std::max_element(min_lower_abs_grad->begin(),
                               min_lower_abs_grad->end());

            HeightDirectionData data;
            data.direction =
              std::distance(min_lower_abs_grad->begin(), max_element);
            data.min_abs_dfdx = *max_element;

            return data;
          }

        return std_cxx17::optional<HeightDirectionData>();
      }



      /**
       * Return true if there are exactly two incoming FunctionBounds and
       * they corresponds to one function being positive definite and
       * one being negative definite. Return false otherwise.
       */
      template <int dim>
      inline bool
      one_positive_one_negative_definite(
        const std::vector<FunctionBounds<dim>> &all_function_bounds)
      {
        if (all_function_bounds.size() != 2)
          return false;
        else
          {
            const FunctionBounds<dim> &bounds0 = all_function_bounds.at(0);
            const FunctionBounds<dim> &bounds1 = all_function_bounds.at(1);

            if (bounds0.value.first > 0 && bounds1.value.second < 0)
              return true;
            if (bounds1.value.first > 0 && bounds0.value.second < 0)
              return true;
            return false;
          }
      }



      /**
       * Transform the points and weights of the incoming quadrature,
       * unit_quadrature, from unit space to the incoming box and add these to
       * quadrature.
       *
       * Note that unit_quadrature should be a quadrature over [0,1]^dim.
       */
      template <int dim>
      void
      map_quadrature_to_box(const Quadrature<dim> &    unit_quadrature,
                            const BoundingBox<dim> &   box,
                            ExtendableQuadrature<dim> &quadrature)
      {
        for (unsigned int i = 0; i < unit_quadrature.size(); ++i)
          {
            const Point<dim> point = box.unit_to_real(unit_quadrature.point(i));
            const double     weight = unit_quadrature.weight(i) * box.volume();

            quadrature.push_back(point, weight);
          }
      }



      /**
       * For each of the incoming dim-dimensional functions, create the
       * restriction to the top and bottom of the incoming BoundingBox and add
       * these two (dim-1)-dimensional functions to @p restrictions. Here, top and bottom is
       * meant with respect to the incoming @p direction. For each function, the
       * "bottom-restriction" will be added before the "top-restriction"
       *
       * @note @p restrictions will be cleared, so after this function
       * restrictions.size() == 2 * functions.size().
       */
      template <int dim>
      void
      restrict_to_top_and_bottom(
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                                                     functions,
        const BoundingBox<dim> &                                box,
        const unsigned int                                      direction,
        std::vector<Functions::CoordinateRestriction<dim - 1>> &restrictions)
      {
        AssertIndexRange(direction, dim);

        restrictions.clear();
        restrictions.reserve(2 * functions.size());

        const double bottom = box.lower_bound(direction);
        const double top    = box.upper_bound(direction);

        for (const auto &function : functions)
          {
            restrictions.push_back(Functions::CoordinateRestriction<dim - 1>(
              function, direction, bottom));
            restrictions.push_back(Functions::CoordinateRestriction<dim - 1>(
              function, direction, top));
          }
      }



      /**
       * Restrict each of the incoming @p functions to @p point,
       * while keeping the coordinate direction @p open_direction open,
       * and add the restriction to @p restrictions.
       *
       * @note @p restrictions will be cleared, so after this function
       * restrictions.size() == functions.size().
       */
      template <int dim>
      void
      restrict_to_point(
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                                                functions,
        const Point<dim - 1> &                             point,
        const unsigned int                                 open_direction,
        std::vector<Functions::PointRestriction<dim - 1>> &restrictions)
      {
        AssertIndexRange(open_direction, dim);

        restrictions.clear();
        restrictions.reserve(functions.size());
        for (const auto &function : functions)
          {
            restrictions.push_back(Functions::PointRestriction<dim - 1>(
              function, open_direction, point));
          }
      }



      /**
       * Let $\{ y_0, ..., y_{n+1} \}$ be such that $[y_0, y_{n+1}]$ is the
       * @p interval and $\{ y_1, ..., y_n \}$ are the @p roots. In each
       * subinterval, $[y_i, y_{i+1}]$, distribute point according to the
       * 1D-quadrature rule $\{(x_q, w_q)\}_q$ (@p quadrature1D).
       * Take the tensor product with the quadrature point $(x, w)$
       * (@p point, @p weight) to create dim-dimensional quadrature points
       * @f[
       * X_q = x_I \times (y_i + (y_{i+1} - y_i) x_q),
       * W_q = w_I (y_{i+1} - y_i) w_q,
       * @f]
       * and add these points to @p q_partitioning.
       */
      template <int dim>
      void
      distribute_points_between_roots(
        const Quadrature<1> &      quadrature1D,
        const BoundingBox<1> &     interval,
        const std::vector<double> &roots,
        const Point<dim - 1> &     point,
        const double               weight,
        const unsigned int         height_function_direction,
        const std::vector<std::reference_wrapper<const Function<1>>>
          &                             level_sets,
        const AdditionalQGeneratorData &additional_data,
        QPartitioning<dim> &            q_partitioning)
      {
        // Make this int to avoid a warning signed/unsigned comparision.
        const int n_roots = roots.size();

        // The number of intervals are roots.size() + 1
        for (int i = -1; i < n_roots; ++i)
          {
            // Start and end point of the subinterval.
            const double start = i < 0 ? interval.lower_bound(0) : roots[i];
            const double end =
              i + 1 < n_roots ? roots[i + 1] : interval.upper_bound(0);

            const double length = end - start;
            // It might be that the end points of the subinterval are roots.
            // If this is the case then the subinterval has length zero.
            // Don't distribute points on the subinterval if it is shorter than
            // some tolerance.
            if (length > additional_data.min_distance_between_roots)
              {
                // All points on the interval belong to the same region in
                // the QPartitioning. Determine the quadrature we should add
                // the points to.
                const Point<1>     center(start + 0.5 * length);
                const Definiteness definiteness =
                  pointwise_definiteness(level_sets, center);
                ExtendableQuadrature<dim> &target_quadrature =
                  q_partitioning.quadrature_by_definiteness(definiteness);

                tensor_point_with_1D_quadrature(point,
                                                weight,
                                                quadrature1D,
                                                start,
                                                end,
                                                height_function_direction,
                                                target_quadrature);
              }
          }
      }



      RootFinder::AdditionalData::AdditionalData(
        const double       tolerance,
        const unsigned int max_recursion_depth,
        const unsigned int max_iterations)
        : tolerance(tolerance)
        , max_recursion_depth(max_recursion_depth)
        , max_iterations(max_iterations)
      {}



      RootFinder::RootFinder(const AdditionalData &data)
        : additional_data(data)
      {}



      void
      RootFinder::find_roots(
        const std::vector<std::reference_wrapper<const Function<1>>> &functions,
        const BoundingBox<1> &                                        interval,
        std::vector<double> &                                         roots)
      {
        for (const Function<1> &function : functions)
          {
            const unsigned int recursion_depth = 0;
            find_roots(function, interval, recursion_depth, roots);
          }
        // Sort and make sure no roots are duplicated
        std::sort(roots.begin(), roots.end());

        const auto roots_are_equal = [this](const double &a, const double &b) {
          return std::abs(a - b) < additional_data.tolerance;
        };
        roots.erase(unique(roots.begin(), roots.end(), roots_are_equal),
                    roots.end());
      }



      void
      RootFinder::find_roots(const Function<1> &   function,
                             const BoundingBox<1> &interval,
                             const unsigned int    recursion_depth,
                             std::vector<double> & roots)
      {
        // Compute function values at end points.
        const double left_value  = function.value(interval.vertex(0));
        const double right_value = function.value(interval.vertex(1));

        // If we have a sign change we solve for the root.
        if (boost::math::sign(left_value) != boost::math::sign(right_value))
          {
            const auto lambda = [&function](const double x) {
              return function.value(Point<1>(x));
            };

            const auto stopping_criteria = [this](const double &a,
                                                  const double &b) {
              return std::abs(a - b) < additional_data.tolerance;
            };

            boost::uintmax_t iterations = additional_data.max_iterations;

            const std::pair<double, double> root_bracket =
              boost::math::tools::toms748_solve(lambda,
                                                interval.lower_bound(0),
                                                interval.upper_bound(0),
                                                left_value,
                                                right_value,
                                                stopping_criteria,
                                                iterations);

            const double root = .5 * (root_bracket.first + root_bracket.second);
            roots.push_back(root);
          }
        else
          {
            // Compute bounds on the incoming function to check if there are
            // roots. If the function is positive or negative on the whole
            // interval we do nothing.
            std::pair<double, double>                value_bounds;
            std::array<std::pair<double, double>, 1> gradient_bounds;
            FunctionTools::taylor_estimate_function_bounds<1>(function,
                                                              interval,
                                                              value_bounds,
                                                              gradient_bounds);

            // Since we already know the function values at the interval ends we
            // might as well check these for min/max too.
            const double function_min =
              std::min(std::min(left_value, right_value), value_bounds.first);

            // If the functions is positive there are no roots.
            if (function_min > 0)
              return;

            const double function_max =
              std::max(std::max(left_value, right_value), value_bounds.second);

            // If the function is negative there are no roots.
            if (function_max < 0)
              return;

            // If we can't say that the function is strictly positive/negative
            // we split the interval in half. We can't split forever, so if we
            // have reached the max recursion, we stop looking for roots.
            if (recursion_depth < additional_data.max_recursion_depth)
              {
                find_roots(function,
                           interval.child(0),
                           recursion_depth + 1,
                           roots);
                find_roots(function,
                           interval.child(1),
                           recursion_depth + 1,
                           roots);
              }
          }
      }



      template <int dim>
      ExtendableQuadrature<dim>::ExtendableQuadrature(
        const Quadrature<dim> &quadrature)
        : Quadrature<dim>(quadrature)
      {}



      template <int dim>
      void
      ExtendableQuadrature<dim>::push_back(const Point<dim> &point,
                                           const double      weight)
      {
        this->quadrature_points.push_back(point);
        this->weights.push_back(weight);
      }



      template <int dim>
      ExtendableQuadrature<dim> &
      QPartitioning<dim>::quadrature_by_definiteness(
        const Definiteness definiteness)
      {
        switch (definiteness)
          {
            case Definiteness::negative:
              return negative;
            case Definiteness::positive:
              return positive;
            default:
              return indefinite;
          }
      }



      /**
       * Takes a (dim-1)-dimensional point from the cross-section (orthogonal
       * to direction) of the box. Creates the two dim-dimensional points, which
       * are the projections from the cross section to the faces of the box and
       * returns the point closest to the zero-contour of the incoming level set
       * function.
       */
      template <int dim>
      Point<dim>
      face_projection_closest_zero_contour(const Point<dim - 1> &  point,
                                           const unsigned int      direction,
                                           const BoundingBox<dim> &box,
                                           const Function<dim> &   level_set)
      {
        const Point<dim> bottom_point =
          dealii::internal::create_higher_dim_point(point,
                                                    direction,
                                                    box.lower_bound(direction));
        const double bottom_value = level_set.value(bottom_point);

        const Point<dim> top_point =
          dealii::internal::create_higher_dim_point(point,
                                                    direction,
                                                    box.upper_bound(direction));
        const double top_value = level_set.value(top_point);

        // The end point closest to the zero-contour is the one with smallest
        // absolute value.
        return std::abs(bottom_value) < std::abs(top_value) ? bottom_point :
                                                              top_point;
      }



      template <int dim, int spacedim>
      UpThroughDimensionCreator<dim, spacedim>::UpThroughDimensionCreator(
        const hp::QCollection<1> &      q_collection1D,
        const AdditionalQGeneratorData &additional_data)
        : q_collection1D(&q_collection1D)
        , additional_data(additional_data)
        , root_finder(
            RootFinder::AdditionalData(additional_data.root_finder_tolerance,
                                       additional_data.max_root_finder_splits))
      {
        q_index = 0;
      }



      template <int dim, int spacedim>
      void
      UpThroughDimensionCreator<dim, spacedim>::generate(
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                        level_sets,
        const BoundingBox<dim> &   box,
        const Quadrature<dim - 1> &low_dim_quadrature,
        const unsigned int         height_function_direction,
        QPartitioning<dim> &       q_partitioning)
      {
        const Quadrature<1> &quadrature1D = (*q_collection1D)[q_index];

        for (unsigned int q = 0; q < low_dim_quadrature.size(); ++q)
          {
            const Point<dim - 1> &point  = low_dim_quadrature.point(q);
            const double          weight = low_dim_quadrature.weight(q);
            restrict_to_point(level_sets,
                              point,
                              height_function_direction,
                              point_restrictions);

            // We need a vector of references to do the recursive call.
            const std::vector<std::reference_wrapper<const Function<1>>>
              restrictions(point_restrictions.begin(),
                           point_restrictions.end());

            const BoundingBox<1> bounds_in_direction =
              box.bounds(height_function_direction);

            roots.clear();
            root_finder.find_roots(restrictions, bounds_in_direction, roots);

            distribute_points_between_roots(quadrature1D,
                                            bounds_in_direction,
                                            roots,
                                            point,
                                            weight,
                                            height_function_direction,
                                            restrictions,
                                            additional_data,
                                            q_partitioning);

            if (dim == spacedim)
              create_surface_point(point,
                                   weight,
                                   level_sets,
                                   box,
                                   height_function_direction,
                                   q_partitioning.surface);
          }

        point_restrictions.clear();
      }



      template <int dim, int spacedim>
      void
      UpThroughDimensionCreator<dim, spacedim>::create_surface_point(
        const Point<dim - 1> &point,
        const double          weight,
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                             level_sets,
        const BoundingBox<dim> &        box,
        const unsigned int              height_function_direction,
        ImmersedSurfaceQuadrature<dim> &surface_quadrature)
      {
        AssertIndexRange(roots.size(), 2);
        Assert(level_sets.size() == 1, ExcInternalError());


        const Function<dim> &level_set = level_sets.at(0);

        Point<dim> surface_point;
        if (roots.size() == 1)
          {
            surface_point = dealii::internal::create_higher_dim_point(
              point, height_function_direction, roots[0]);
          }
        else
          {
            // If we got here, we have missed roots in the lower dimensional
            // algorithm. This is a rare event but can happen if the
            // zero-contour has a high curvature. The problem is that the
            // incoming point has been incorrectly added to the indefinite
            // quadrature in QPartitioning<dim-1>. Since we missed a root on
            // this box, we will likely miss it on the neighboring box too. If
            // this happens, the point will NOT be in the indefinite quadrature
            // on the neighbor. The best thing we can do is to compute the
            // surface point by projecting the lower dimensional point to the
            // face closest to the zero-contour. We should add a surface point
            // because the neighbor will not.
            surface_point = face_projection_closest_zero_contour(
              point, height_function_direction, box, level_set);
          }

        const Tensor<1, dim> gradient = level_set.gradient(surface_point);
        Tensor<1, dim>       normal   = gradient;
        normal *= 1. / normal.norm();

        // Note that gradient[height_function_direction] is non-zero
        // because of the implicit function theorem.
        const double surface_weight =
          weight * gradient.norm() /
          std::abs(gradient[height_function_direction]);
        surface_quadrature.push_back(surface_point, surface_weight, normal);
      }



      template <int dim, int spacedim>
      void
      UpThroughDimensionCreator<dim, spacedim>::set_1D_quadrature(
        unsigned int q_index)
      {
        AssertIndexRange(q_index, q_collection1D->size());
        this->q_index = q_index;
      }



      template <int dim, int spacedim>
      QGeneratorBase<dim, spacedim>::QGeneratorBase(
        const hp::QCollection<1> &      q_collection1D,
        const AdditionalQGeneratorData &additional_data)
        : additional_data(additional_data)
        , q_collection1D(&q_collection1D)
      {
        q_index = 0;
      }



      template <int dim, int spacedim>
      QGenerator<dim, spacedim>::QGenerator(
        const hp::QCollection<1> &      q_collection1D,
        const AdditionalQGeneratorData &additional_data)
        : QGeneratorBase<dim, spacedim>(q_collection1D, additional_data)
        , low_dim_algorithm(q_collection1D, additional_data)
        , up_through_dimension_creator(q_collection1D, additional_data)
      {
        for (unsigned int i = 0; i < q_collection1D.size(); ++i)
          tensor_products.push_back(Quadrature<dim>(q_collection1D[i]));
      }



      template <int dim, int spacedim>
      void
      QGeneratorBase<dim, spacedim>::clear_quadratures()
      {
        q_partitioning = QPartitioning<dim>();
      }



      template <int dim, int spacedim>
      const QPartitioning<dim> &
      QGeneratorBase<dim, spacedim>::get_quadratures() const
      {
        return q_partitioning;
      }



      template <int dim, int spacedim>
      void
      QGenerator<dim, spacedim>::generate(
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                     level_sets,
        const BoundingBox<dim> &box,
        const unsigned int      n_box_splits)
      {
        std::vector<FunctionBounds<dim>> all_function_bounds;
        estimate_function_bounds(level_sets, box, all_function_bounds);

        const std::pair<double, double> extreme_values =
          find_extreme_values(all_function_bounds);

        if (extreme_values.first > this->additional_data.limit_to_be_definite)
          {
            map_quadrature_to_box(tensor_products[this->q_index],
                                  box,
                                  this->q_partitioning.positive);
          }
        else if (extreme_values.second <
                 -(this->additional_data.limit_to_be_definite))
          {
            map_quadrature_to_box(tensor_products[this->q_index],
                                  box,
                                  this->q_partitioning.negative);
          }
        else if (one_positive_one_negative_definite(all_function_bounds))
          {
            map_quadrature_to_box(tensor_products[this->q_index],
                                  box,
                                  this->q_partitioning.indefinite);
          }
        else
          {
            const std_cxx17::optional<HeightDirectionData> data =
              find_best_height_direction(all_function_bounds);

            // Check larger than a constant to avoid that min_abs_dfdx is only
            // larger by 0 by floating point precision.
            if (data && data->min_abs_dfdx >
                          this->additional_data.lower_bound_implicit_function)
              {
                create_low_dim_quadratures(data->direction,
                                           level_sets,
                                           box,
                                           n_box_splits);
                create_high_dim_quadratures(data->direction, level_sets, box);
              }
            else if (n_box_splits < this->additional_data.max_box_splits)
              {
                split_box_and_recurse(level_sets, box, data, n_box_splits);
              }
            else
              {
                // We can't split the box recursively forever. Use the midpoint
                // method as a last resort.
                use_midpoint_method(level_sets, box);
              }
          }
      }



      /**
       * Return the coordinate direction of the largest side of the box.
       * If two or more sides have the same length the returned std::optional
       * will be non-set.
       */
      template <int dim>
      std_cxx17::optional<unsigned int>
      direction_of_largest_extent(const BoundingBox<dim> &box)
      {
        // Get the side lengths for each direction and sort them.
        std::array<std::pair<double, unsigned int>, dim> side_lengths;
        for (unsigned int i = 0; i < dim; ++i)
          {
            side_lengths[i].first  = box.side_length(i);
            side_lengths[i].second = i;
          }
        // Sort is lexicographic, so this sorts based on side length first.
        std::sort(side_lengths.begin(), side_lengths.end());

        // Check if the two largest side lengths have the same length. This
        // function isn't called in 1D, so the (dim - 2)-element exists.
        if (boost::math::epsilon_difference(side_lengths[dim - 1].first,
                                            side_lengths[dim - 2].first) < 100)
          return std_cxx17::optional<unsigned int>();

        return side_lengths.back().second;
      }



      /**
       * Return the coordinate direction that the box should be split in,
       * assuming that the box should be split it half.
       *
       * If the box is larger in one coordante direction, this direction is
       * returned. If the box have the same extent in all directions, we choose
       * the coordinate direction which is closest to being a height-function
       * direction. That is, the direction $i$ that has a least negative
       * estimate of $|\partial_i \psi_j|$. As a last resort, we choose the
       * direction 0, if @p height_direction_data non-set.
       */
      template <int dim>
      unsigned int
      compute_split_direction(
        const BoundingBox<dim> &                        box,
        const std_cxx17::optional<HeightDirectionData> &height_direction_data)
      {
        const std_cxx17::optional<unsigned int> direction =
          direction_of_largest_extent(box);

        if (direction)
          return *direction;

        // This direction is closest to being a height direction, so
        // we split in this direction.
        if (height_direction_data)
          return height_direction_data->direction;

        // We have to choose some direction, we might aswell take 0.
        return 0;
      }



      /**
       * Split the incoming box in half with respect to the incoming coordinate
       * direction and return the left half.
       */
      template <int dim>
      inline BoundingBox<dim>
      left_half(const BoundingBox<dim> &box, const unsigned int direction)
      {
        AssertIndexRange(direction, dim);

        // Move the upper corner half a side-length to the left.
        std::pair<Point<dim>, Point<dim>> corners = box.get_boundary_points();
        corners.second[direction] -= .5 * box.side_length(direction);

        return BoundingBox<dim>(corners);
      }



      /**
       * Split the incoming box in half with respect to the incoming coordinate
       * direction and return the right half.
       */
      template <int dim>
      inline BoundingBox<dim>
      right_half(const BoundingBox<dim> &box, const unsigned int direction)
      {
        AssertIndexRange(direction, dim);

        // Move the lower corner half a side-length to the right.
        std::pair<Point<dim>, Point<dim>> corners = box.get_boundary_points();
        corners.first[direction] += .5 * box.side_length(direction);

        return BoundingBox<dim>(corners);
      }



      template <int dim, int spacedim>
      void
      QGenerator<dim, spacedim>::split_box_and_recurse(
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                                             level_sets,
        const BoundingBox<dim> &                        box,
        const std_cxx17::optional<HeightDirectionData> &direction_data,
        const unsigned int                              n_box_splits)
      {
        if (this->additional_data.split_in_half)
          {
            const unsigned int direction =
              compute_split_direction(box, direction_data);

            const BoundingBox<dim> left_box  = left_half(box, direction);
            const BoundingBox<dim> right_box = right_half(box, direction);

            generate(level_sets, left_box, n_box_splits + 1);
            generate(level_sets, right_box, n_box_splits + 1);
          }
        else
          {
            for (unsigned int i = 0;
                 i < GeometryInfo<dim>::max_children_per_cell;
                 ++i)
              {
                generate(level_sets, box.child(i), n_box_splits + 1);
              }
          }
      }



      template <int dim, int spacedim>
      void
      QGenerator<dim, spacedim>::create_low_dim_quadratures(
        const unsigned int height_function_direction,
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                     level_sets,
        const BoundingBox<dim> &box,
        const unsigned int      n_box_splits)
      {
        std::vector<Functions::CoordinateRestriction<dim - 1>>
          face_restrictions;
        restrict_to_top_and_bottom(level_sets,
                                   box,
                                   height_function_direction,
                                   face_restrictions);

        // We need a vector of references to do the recursive call.
        const std::vector<std::reference_wrapper<const Function<dim - 1>>>
          restrictions(face_restrictions.begin(), face_restrictions.end());

        const BoundingBox<dim - 1> cross_section =
          box.cross_section(height_function_direction);

        low_dim_algorithm.clear_quadratures();
        low_dim_algorithm.generate(restrictions, cross_section, n_box_splits);
      }



      template <int dim, int spacedim>
      void
      QGenerator<dim, spacedim>::create_high_dim_quadratures(
        const unsigned int height_function_direction,
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                     level_sets,
        const BoundingBox<dim> &box)
      {
        const QPartitioning<dim - 1> &low_dim_quadratures =
          low_dim_algorithm.get_quadratures();

        const Quadrature<1> &quadrature1D =
          (*this->q_collection1D)[this->q_index];

        add_tensor_product(low_dim_quadratures.negative,
                           quadrature1D,
                           box.lower_bound(height_function_direction),
                           box.upper_bound(height_function_direction),
                           height_function_direction,
                           this->q_partitioning.negative);

        add_tensor_product(low_dim_quadratures.positive,
                           quadrature1D,
                           box.lower_bound(height_function_direction),
                           box.upper_bound(height_function_direction),
                           height_function_direction,
                           this->q_partitioning.positive);

        up_through_dimension_creator.generate(level_sets,
                                              box,
                                              low_dim_quadratures.indefinite,
                                              height_function_direction,
                                              this->q_partitioning);
      }



      template <int dim, int spacedim>
      void
      QGenerator<dim, spacedim>::use_midpoint_method(
        const std::vector<std::reference_wrapper<const Function<dim>>>
          &                     level_sets,
        const BoundingBox<dim> &box)
      {
        const Point<dim>   center = box.center();
        const Definiteness definiteness =
          pointwise_definiteness(level_sets, center);

        ExtendableQuadrature<dim> &quadrature =
          this->q_partitioning.quadrature_by_definiteness(definiteness);

        quadrature.push_back(center, box.volume());
      }



      template <int dim, int spacedim>
      void
      QGenerator<dim, spacedim>::set_1D_quadrature(const unsigned int q_index)
      {
        AssertIndexRange(q_index, this->q_collection1D->size());

        this->q_index = q_index;
        low_dim_algorithm.set_1D_quadrature(q_index);
        up_through_dimension_creator.set_1D_quadrature(q_index);
      }



      template <int spacedim>
      QGenerator<1, spacedim>::QGenerator(
        const hp::QCollection<1> &      q_collection1D,
        const AdditionalQGeneratorData &additional_data)
        : QGeneratorBase<1, spacedim>(q_collection1D, additional_data)
        , root_finder(
            RootFinder::AdditionalData(additional_data.root_finder_tolerance,
                                       additional_data.max_root_finder_splits))
      {
        Assert(q_collection1D.size() > 0,
               ExcMessage("Incoming quadrature collection is empty."));
      }



      template <int spacedim>
      void
      QGenerator<1, spacedim>::generate(
        const std::vector<std::reference_wrapper<const Function<1>>>
          &                   level_sets,
        const BoundingBox<1> &box,
        const unsigned int    n_box_splits)
      {
        (void)n_box_splits;

        roots.clear();
        root_finder.find_roots(level_sets, box, roots);

        const Quadrature<1> &quadrature1D =
          (*this->q_collection1D)[this->q_index];

        distribute_points_between_roots(quadrature1D,
                                        box,
                                        roots,
                                        zero_dim_point,
                                        unit_weight,
                                        direction,
                                        level_sets,
                                        this->additional_data,
                                        this->q_partitioning);

        if (spacedim == 1)
          this->create_surface_points(level_sets);
      }



      template <int spacedim>
      void
      QGenerator<1, spacedim>::create_surface_points(
        const std::vector<std::reference_wrapper<const Function<1>>>
          &level_sets)
      {
        Assert(level_sets.size() == 1, ExcInternalError());

        for (const double root : roots)
          {
            // A surface integral in 1D is just a point evaluation,
            // so the weight is always 1.
            const double   weight = 1;
            const Point<1> point(root);

            Tensor<1, 1> normal        = level_sets[0].get().gradient(point);
            const double gradient_norm = normal.norm();
            Assert(
              gradient_norm > 1e-11,
              ExcMessage(
                "The level set function has a gradient almost equal to 0."));
            normal *= 1. / gradient_norm;

            this->q_partitioning.surface.push_back(point, weight, normal);
          }
      }



      template <int spacedim>
      void
      QGenerator<1, spacedim>::set_1D_quadrature(const unsigned int q_index)
      {
        AssertIndexRange(q_index, this->q_collection1D->size());
        this->q_index = q_index;
      }
    } // namespace QuadratureGeneratorImplementation
  }   // namespace internal



  AdditionalQGeneratorData::AdditionalQGeneratorData(
    const unsigned int max_box_splits,
    const double       lower_bound_implicit_function,
    const double       min_distance_between_roots,
    const double       limit_to_be_definite,
    const double       root_finder_tolerance,
    const unsigned int max_root_finder_splits,
    bool               split_in_half)
    : max_box_splits(max_box_splits)
    , lower_bound_implicit_function(lower_bound_implicit_function)
    , min_distance_between_roots(min_distance_between_roots)
    , limit_to_be_definite(limit_to_be_definite)
    , root_finder_tolerance(root_finder_tolerance)
    , max_root_finder_splits(max_root_finder_splits)
    , split_in_half(split_in_half)
  {}



  template <int dim>
  QuadratureGenerator<dim>::QuadratureGenerator(
    const hp::QCollection<1> &q_collection,
    const AdditionalData &    additional_data)
    : q_generator(q_collection, additional_data)
  {
    Assert(q_collection.size() > 0,
           ExcMessage("Incoming hp::QCollection<1> is empty."));
  }



  template <int dim>
  void
  QuadratureGenerator<dim>::generate(const Function<dim> &   level_set,
                                     const BoundingBox<dim> &box)
  {
    Assert(level_set.n_components == 1,
           ExcMessage(
             "The incoming function should be a scalar level set function,"
             " it should have one component."));
    Assert(box.volume() > 0, ExcMessage("Incoming box has zero volume."));

    q_generator.clear_quadratures();

    std::vector<std::reference_wrapper<const Function<dim>>> level_sets;
    level_sets.push_back(level_set);

    const unsigned int n_box_splits = 0;
    q_generator.generate(level_sets, box, n_box_splits);

    // With a single level set function, the "indefinite" quadrature should be
    // zero. If you call generate() with a ZeroFunction nothing good can be
    // done. You will end up here.
    Assert(
      q_generator.get_quadratures().indefinite.size() == 0,
      ExcMessage(
        "Generation of quadrature rules failed. This can mean that the level"
        "set function is degenerate in some way, e.g. oscillating extremely"
        "rapidly."));
  }



  template <int dim>
  const Quadrature<dim> &
  QuadratureGenerator<dim>::get_inside_quadrature() const
  {
    return q_generator.get_quadratures().negative;
  }



  template <int dim>
  const Quadrature<dim> &
  QuadratureGenerator<dim>::get_outside_quadrature() const
  {
    return q_generator.get_quadratures().positive;
  }



  template <int dim>
  const ImmersedSurfaceQuadrature<dim> &
  QuadratureGenerator<dim>::get_surface_quadrature() const
  {
    return q_generator.get_quadratures().surface;
  }


  template <int dim>
  void
  QuadratureGenerator<dim>::set_1D_quadrature(const unsigned int q_index)
  {
    q_generator.set_1D_quadrature(q_index);
  }



  template <int dim>
  FaceQuadratureGenerator<dim>::FaceQuadratureGenerator(
    const hp::QCollection<1> &quadratures1D,
    const AdditionalData &    additional_data)
    : quadrature_generator(quadratures1D, additional_data)
  {}



  template <int dim>
  void
  FaceQuadratureGenerator<dim>::generate(const Function<dim> &   level_set,
                                         const BoundingBox<dim> &box,
                                         const unsigned int      face_index)
  {
    AssertIndexRange(face_index, GeometryInfo<dim>::faces_per_cell);

    // We restrict the level set function to the face, by locking the coordinate
    // that is constant over the face. This will be the same as the direction of
    // the face normal.
    const unsigned int face_normal_direction =
      GeometryInfo<dim>::unit_normal_direction[face_index];

    const Point<dim> vertex0 =
      box.vertex(GeometryInfo<dim>::face_to_cell_vertices(face_index, 0));
    const double coordinate_value = vertex0(face_normal_direction);

    const Functions::CoordinateRestriction<dim - 1> face_restriction(
      level_set, face_normal_direction, coordinate_value);

    // Reuse the lower dimensional QuadratureGenerator on the face.
    const BoundingBox<dim - 1> cross_section =
      box.cross_section(face_normal_direction);
    quadrature_generator.generate(face_restriction, cross_section);

    // We need the dim-dimensional normals of the zero-contour.
    // Recompute these.
    const ImmersedSurfaceQuadrature<dim - 1, dim - 1>
      &surface_quadrature_wrong_normal =
        quadrature_generator.get_surface_quadrature();

    std::vector<Tensor<1, dim>> normals;
    normals.reserve(surface_quadrature_wrong_normal.size());
    for (unsigned int i = 0; i < surface_quadrature_wrong_normal.size(); ++i)
      {
        const Point<dim> point = dealii::internal::create_higher_dim_point(
          surface_quadrature_wrong_normal.point(i),
          face_normal_direction,
          coordinate_value);

        Tensor<1, dim> normal = level_set.gradient(point);
        normal /= normal.norm();
        normals.push_back(normal);
      }
    surface_quadrature = ImmersedSurfaceQuadrature<dim - 1, dim>(
      surface_quadrature_wrong_normal.get_points(),
      surface_quadrature_wrong_normal.get_weights(),
      normals);
  }



  template <int dim>
  void
  FaceQuadratureGenerator<dim>::set_1D_quadrature(const unsigned int q_index)
  {
    quadrature_generator.set_1D_quadrature(q_index);
  }



  template <int dim>
  const Quadrature<dim - 1> &
  FaceQuadratureGenerator<dim>::get_inside_quadrature() const
  {
    return quadrature_generator.get_inside_quadrature();
  }


  template <int dim>
  const Quadrature<dim - 1> &
  FaceQuadratureGenerator<dim>::get_outside_quadrature() const
  {
    return quadrature_generator.get_outside_quadrature();
  }



  template <int dim>
  const ImmersedSurfaceQuadrature<dim - 1, dim> &
  FaceQuadratureGenerator<dim>::get_surface_quadrature() const
  {
    return surface_quadrature;
  }
} // namespace NonMatching
#include "quadrature_generator.inst"
DEAL_II_NAMESPACE_CLOSE
