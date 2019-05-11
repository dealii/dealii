// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_fe_field_convolution_function_templates_h
#define dealii_fe_field_convolution_function_templates_h


#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/meshworker/scratch_data.h>

#include <deal.II/numerics/fe_field_convolution_function.h>
#include <deal.II/numerics/vector_tools.h>

#include <boost/geometry.hpp>

namespace bgi = boost::geometry::index;

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    FEFieldConvolutionFunction(const DoFHandlerType &                   mydh,
                               const GridTools::Cache<dim, spacedim> &  cache,
                               const VectorType &                       myv,
                               Functions::CutOffFunctionBase<spacedim> &kernel,
                               const Quadrature<dim> &quadrature)
    : Function<spacedim, typename VectorType::value_type>(
        mydh.get_fe(0).n_components())
    , dh(&mydh, "FEFieldConvolutionFunction")
    , cache(&cache, "FEFieldConvolutionFunction")
    , data_vector(myv)
    , kernel(&kernel, "FEFieldConvolutionFunction")
    , scratch(cache.get_mapping(),
              dh->get_fe(),
              quadrature,
              update_values | update_JxW_values | update_quadrature_points)
  {}



  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    vector_value(const Point<spacedim> &                  p,
                 Vector<typename VectorType::value_type> &values) const
  {
    Assert(values.size() == this->n_components,
           ExcDimensionMismatch(values.size(), this->n_components));

    std::vector<Point<spacedim>>                         points({p});
    std::vector<Vector<typename VectorType::value_type>> value_list({values});
    vector_value_list(points, value_list);
    values = value_list[0];
  }



  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  typename VectorType::value_type
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::value(
    const Point<spacedim> &p,
    const unsigned int     comp) const
  {
    Vector<typename VectorType::value_type> values(this->n_components);
    vector_value(p, values);
    return values(comp);
  }



  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    vector_gradient(
      const Point<spacedim> &p,
      std::vector<Tensor<1, spacedim, typename VectorType::value_type>>
        &gradients) const
  {
    Assert(gradients.size() == this->n_components,
           ExcDimensionMismatch(gradients.size(), this->n_components));

    std::vector<Point<spacedim>> points({p});
    std::vector<
      std::vector<Tensor<1, spacedim, typename VectorType::value_type>>>
      gradient_list({gradients});
    vector_gradient_list(points, gradient_list);
    gradients = gradient_list[0];
  }



  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  Tensor<1, spacedim, typename VectorType::value_type>
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    gradient(const Point<spacedim> &p, const unsigned int comp) const
  {
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> grads(
      this->n_components);
    vector_gradient(p, grads);
    return grads[comp];
  }

  // Now the list versions
  // ==============================

  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    vector_value_list(
      const std::vector<Point<spacedim>> &                  points,
      std::vector<Vector<typename VectorType::value_type>> &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    using PointLeaf = std::pair<BoundingBox<spacedim>, unsigned int>;

    std::vector<PointLeaf> point_boxes(points.size());

    for (unsigned int i = 0; i < points.size(); ++i)
      {
        const auto &          p = points[i];
        BoundingBox<spacedim> box({p, p});
        box.extend(kernel->get_radius());
        point_boxes[i] = {box, i};
        // Reset also all values to zero
        values[i] = 0;
      }

    const auto  points_tree = pack_rtree(point_boxes);
    const auto &cells_tree  = cache->get_cell_bounding_boxes_rtree();

    std::vector<bool> visited_cells(cells_tree.size(), false);

    AssertDimension(cells_tree.size(),
                    dh->get_triangulation().n_active_cells());

    std::vector<Vector<typename VectorType::value_type>> field_values;
    std::vector<double>                                  kernel_values;

    for (const auto &box_and_id : points_tree)
      for (const auto &box_and_cell :
           cells_tree |
             bgi::adaptors::queried(bgi::intersects(box_and_id.first)))
        if (box_and_cell.second->is_locally_owned() &&
            visited_cells[box_and_cell.second->index()] == false)
          {
            const auto cell =
              typename DoFHandlerType::cell_iterator(*box_and_cell.second, dh);
            const auto &fev      = scratch.reinit(cell);
            const auto &q_points = scratch.get_quadrature_points();
            const auto &JxW      = scratch.get_JxW_values();

            field_values.resize(q_points.size(),
                                Vector<typename VectorType::value_type>(
                                  dh->get_fe().n_components()));

            kernel_values.resize(q_points.size());

            fev.get_function_values(data_vector, field_values);
            // We now gather all other points that are effectively at distance
            // less than or equal to radius from this cell
            for (const auto &box_and_id :
                 points_tree |
                   bgi::adaptors::queried(bgi::intersects(box_and_cell.first)))
              {
                const auto &point_id = box_and_id.second;
                const auto &p        = points[point_id];
                kernel->set_center(p);
                kernel->value_list(q_points, kernel_values);
                for (unsigned int q = 0; q < q_points.size(); ++q)
                  for (unsigned int i = 0; i < dh->get_fe().n_components(); ++i)
                    {
                      values[point_id][i] +=
                        JxW[q] * field_values[q][i] * kernel_values[q];
                    }
              }
            // Mark this cell as visited. All overlapping points have been taken
            // care of, and we can skip this cell for all future intersections
            visited_cells[box_and_cell.second->index()] = true;
          }
  }



  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    value_list(const std::vector<Point<spacedim>> &          points,
               std::vector<typename VectorType::value_type> &values,
               const unsigned int                            component) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    // Simply forward everything to the vector_value_list()
    // function. This requires a temporary object, but everything we
    // do here is so expensive that that really doesn't make any
    // difference any more.
    std::vector<Vector<typename VectorType::value_type>> vvalues(
      points.size(),
      Vector<typename VectorType::value_type>(this->n_components));

    vector_value_list(points, vvalues);

    for (unsigned int q = 0; q < points.size(); ++q)
      values[q] = vvalues[q](component);
  }



  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    vector_gradient_list(
      const std::vector<Point<spacedim>> &,
      std::vector<
        std::vector<Tensor<1, spacedim, typename VectorType::value_type>>> &)
      const
  {
    Assert(false, ExcNotImplemented());
  }



  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    gradient_list(
      const std::vector<Point<spacedim>> &                               points,
      std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &values,
      const unsigned int component) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    // Simply forward everything to the vector_gradient_list()
    // function. This requires a temporary object, but everything we
    // do here is so expensive that that really doesn't make any
    // difference any more.
    std::vector<
      std::vector<Tensor<1, spacedim, typename VectorType::value_type>>>
      vvalues(points.size(),
              std::vector<Tensor<1, spacedim, typename VectorType::value_type>>(
                this->n_components));

    vector_gradient_list(points, vvalues);

    for (unsigned int q = 0; q < points.size(); ++q)
      values[q] = vvalues[q][component];
  }



  template <int dim, int spacedim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldConvolutionFunction<dim, spacedim, DoFHandlerType, VectorType>::
    set_radius(const double &value)
  {
    kernel->set_radius(value);
  }
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
