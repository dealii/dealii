// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2020 by the deal.II authors
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

#ifndef dealii_fe_field_function_templates_h
#define dealii_fe_field_function_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools_common.h>

#include <tuple>



DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  template <int dim, typename DoFHandlerType, typename VectorType>
  FEFieldFunction<dim, DoFHandlerType, VectorType>::FEFieldFunction(
    const DoFHandlerType &mydh,
    const VectorType &    myv,
    const Mapping<dim> &  mymapping)
    : Function<dim, typename VectorType::value_type>(
        mydh.get_fe(0).n_components())
    , dh(&mydh, "FEFieldFunction")
    , data_vector(myv)
    , mapping(mymapping)
    , cache(dh->get_triangulation(), mymapping)
    , cell_hint(dh->end())
  {}



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::set_active_cell(
    const typename DoFHandlerType::active_cell_iterator &newcell)
  {
    cell_hint.get() = newcell;
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::vector_value(
    const Point<dim> &                       p,
    Vector<typename VectorType::value_type> &values) const
  {
    Assert(values.size() == this->n_components,
           ExcDimensionMismatch(values.size(), this->n_components));
    typename DoFHandlerType::active_cell_iterator cell = cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    std_cxx17::optional<Point<dim>> qp = get_reference_coordinates(cell, p);
    if (!qp)
      {
        const std::pair<typename dealii::internal::
                          ActiveCellIterator<dim, dim, DoFHandlerType>::type,
                        Point<dim>>
          my_pair = GridTools::find_active_cell_around_point(mapping, *dh, p);
        AssertThrow(!my_pair.first->is_artificial(),
                    VectorTools::ExcPointNotAvailableHere());

        cell = my_pair.first;
        qp   = my_pair.second;
      }

    cell_hint.get() = cell;

    // check that the current cell is available:
    AssertThrow(!cell->is_artificial(),
                VectorTools::ExcPointNotAvailableHere());

    // Now we can find out about the point
    Quadrature<dim> quad = *qp;
    FEValues<dim>   fe_v(mapping, cell->get_fe(), quad, update_values);
    fe_v.reinit(cell);
    std::vector<Vector<typename VectorType::value_type>> vvalues(
      1, Vector<typename VectorType::value_type>(values.size()));
    fe_v.get_function_values(data_vector, vvalues);
    values = vvalues[0];
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  typename VectorType::value_type
  FEFieldFunction<dim, DoFHandlerType, VectorType>::value(
    const Point<dim> & p,
    const unsigned int comp) const
  {
    Vector<typename VectorType::value_type> values(this->n_components);
    vector_value(p, values);
    return values(comp);
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::vector_gradient(
    const Point<dim> &                                            p,
    std::vector<Tensor<1, dim, typename VectorType::value_type>> &gradients)
    const
  {
    using number = typename VectorType::value_type;
    Assert(gradients.size() == this->n_components,
           ExcDimensionMismatch(gradients.size(), this->n_components));
    typename DoFHandlerType::active_cell_iterator cell = cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    std_cxx17::optional<Point<dim>> qp = get_reference_coordinates(cell, p);
    if (!qp)
      {
        const std::pair<typename dealii::internal::
                          ActiveCellIterator<dim, dim, DoFHandlerType>::type,
                        Point<dim>>
          my_pair = GridTools::find_active_cell_around_point(mapping, *dh, p);
        AssertThrow(!my_pair.first->is_artificial(),
                    VectorTools::ExcPointNotAvailableHere());

        cell = my_pair.first;
        qp   = my_pair.second;
      }

    // check that the current cell is available:
    AssertThrow(!cell->is_artificial(),
                VectorTools::ExcPointNotAvailableHere());

    cell_hint.get() = cell;

    // Now we can find out about the point
    Quadrature<dim> quad = *qp;
    FEValues<dim>   fe_v(mapping, cell->get_fe(), quad, update_gradients);
    fe_v.reinit(cell);

    if (this->n_components == 1)
      {
        // the size of the @p gradients coincidentally coincides
        // with the number of quadrature points we evaluate the function at.
        fe_v.get_function_gradients(data_vector, gradients);
      }
    else
      {
        // Unfortunately we still need a temporary argument as we want to
        // evaluate a gradient of a (generally) multicomponent function at
        // a single quadrature point. Note that the first std::vector<> is
        // related to the number of quadrature points (always one here), whereas
        // the second to the number of components.
        std::vector<std::vector<Tensor<1, dim, number>>> vgrads(
          1, std::vector<Tensor<1, dim, number>>(this->n_components));
        fe_v.get_function_gradients(data_vector, vgrads);
        gradients = vgrads[0];
      }
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  Tensor<1, dim, typename VectorType::value_type>
  FEFieldFunction<dim, DoFHandlerType, VectorType>::gradient(
    const Point<dim> & p,
    const unsigned int comp) const
  {
    std::vector<Tensor<1, dim, typename VectorType::value_type>> grads(
      this->n_components);
    vector_gradient(p, grads);
    return grads[comp];
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::vector_laplacian(
    const Point<dim> &                       p,
    Vector<typename VectorType::value_type> &values) const
  {
    Assert(values.size() == this->n_components,
           ExcDimensionMismatch(values.size(), this->n_components));
    typename DoFHandlerType::active_cell_iterator cell = cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    std_cxx17::optional<Point<dim>> qp = get_reference_coordinates(cell, p);
    if (!qp)
      {
        const std::pair<typename dealii::internal::
                          ActiveCellIterator<dim, dim, DoFHandlerType>::type,
                        Point<dim>>
          my_pair = GridTools::find_active_cell_around_point(mapping, *dh, p);
        AssertThrow(!my_pair.first->is_artificial(),
                    VectorTools::ExcPointNotAvailableHere());

        cell = my_pair.first;
        qp   = my_pair.second;
      }

    // check that the current cell is available:
    AssertThrow(!cell->is_artificial(),
                VectorTools::ExcPointNotAvailableHere());

    cell_hint.get() = cell;

    // Now we can find out about the point
    Quadrature<dim> quad = *qp;
    FEValues<dim>   fe_v(mapping, cell->get_fe(), quad, update_hessians);
    fe_v.reinit(cell);
    std::vector<Vector<typename VectorType::value_type>> vvalues(
      1, Vector<typename VectorType::value_type>(values.size()));
    fe_v.get_function_laplacians(data_vector, vvalues);
    values = vvalues[0];
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  typename VectorType::value_type
  FEFieldFunction<dim, DoFHandlerType, VectorType>::laplacian(
    const Point<dim> & p,
    const unsigned int comp) const
  {
    Vector<typename VectorType::value_type> lap(this->n_components);
    vector_laplacian(p, lap);
    return lap[comp];
  }


  // Now the list versions
  // ==============================

  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::vector_value_list(
    const std::vector<Point<dim>> &                       points,
    std::vector<Vector<typename VectorType::value_type>> &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandlerType::active_cell_iterator> cells;
    std::vector<std::vector<Point<dim>>>                       qpoints;
    std::vector<std::vector<unsigned int>>                     maps;

    const unsigned int n_cells =
      compute_point_locations(points, cells, qpoints, maps);

    // Create quadrature collection
    hp::QCollection<dim> quadrature_collection;
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        // Number of quadrature points on this cell
        unsigned int nq = qpoints[i].size();
        // Construct a quadrature formula
        std::vector<double> ww(nq, 1. / nq);

        quadrature_collection.push_back(Quadrature<dim>(qpoints[i], ww));
      }

    // Now gather all the information we need
    hp::MappingCollection<dim> mapping_collection(mapping);
    hp::FEValues<dim>          fe_v(mapping_collection,
                           dh->get_fe_collection(),
                           quadrature_collection,
                           update_values);
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        AssertThrow(!cells[i]->is_artificial(),
                    VectorTools::ExcPointNotAvailableHere());

        fe_v.reinit(cells[i], i, 0);
        const unsigned int nq = qpoints[i].size();

        std::vector<Vector<typename VectorType::value_type>> vvalues(
          nq, Vector<typename VectorType::value_type>(this->n_components));
        fe_v.get_present_fe_values().get_function_values(data_vector, vvalues);

        for (unsigned int q = 0; q < nq; ++q)
          values[maps[i][q]] = vvalues[q];
      }
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::value_list(
    const std::vector<Point<dim>> &               points,
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



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::vector_gradient_list(
    const std::vector<Point<dim>> &points,
    std::vector<std::vector<Tensor<1, dim, typename VectorType::value_type>>>
      &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandlerType::active_cell_iterator> cells;
    std::vector<std::vector<Point<dim>>>                       qpoints;
    std::vector<std::vector<unsigned int>>                     maps;

    const unsigned int n_cells =
      compute_point_locations(points, cells, qpoints, maps);

    // Create quadrature collection
    hp::QCollection<dim> quadrature_collection;
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        // Number of quadrature points on this cell
        unsigned int nq = qpoints[i].size();
        // Construct a quadrature formula
        std::vector<double> ww(nq, 1. / nq);

        quadrature_collection.push_back(Quadrature<dim>(qpoints[i], ww));
      }

    // Now gather all the information we need
    hp::MappingCollection<dim> mapping_collection(mapping);
    hp::FEValues<dim>          fe_v(mapping_collection,
                           dh->get_fe_collection(),
                           quadrature_collection,
                           update_gradients);
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        AssertThrow(!cells[i]->is_artificial(),
                    VectorTools::ExcPointNotAvailableHere());

        fe_v.reinit(cells[i], i, 0);

        const unsigned int nq = qpoints[i].size();
        std::vector<
          std::vector<Tensor<1, dim, typename VectorType::value_type>>>
          vgrads(nq,
                 std::vector<Tensor<1, dim, typename VectorType::value_type>>(
                   this->n_components));
        fe_v.get_present_fe_values().get_function_gradients(data_vector,
                                                            vgrads);

        for (unsigned int q = 0; q < nq; ++q)
          {
            const unsigned int s = vgrads[q].size();
            values[maps[i][q]].resize(s);
            for (unsigned int l = 0; l < s; l++)
              values[maps[i][q]][l] = vgrads[q][l];
          }
      }
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::gradient_list(
    const std::vector<Point<dim>> &                               points,
    std::vector<Tensor<1, dim, typename VectorType::value_type>> &values,
    const unsigned int component) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    // Simply forward everything to the vector_gradient_list()
    // function. This requires a temporary object, but everything we
    // do here is so expensive that that really doesn't make any
    // difference any more.
    std::vector<std::vector<Tensor<1, dim, typename VectorType::value_type>>>
      vvalues(points.size(),
              std::vector<Tensor<1, dim, typename VectorType::value_type>>(
                this->n_components));

    vector_gradient_list(points, vvalues);

    for (unsigned int q = 0; q < points.size(); ++q)
      values[q] = vvalues[q][component];
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::vector_laplacian_list(
    const std::vector<Point<dim>> &                       points,
    std::vector<Vector<typename VectorType::value_type>> &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandlerType::active_cell_iterator> cells;
    std::vector<std::vector<Point<dim>>>                       qpoints;
    std::vector<std::vector<unsigned int>>                     maps;

    const unsigned int n_cells =
      compute_point_locations(points, cells, qpoints, maps);

    // Create quadrature collection
    hp::QCollection<dim> quadrature_collection;
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        // Number of quadrature points on this cell
        unsigned int nq = qpoints[i].size();
        // Construct a quadrature formula
        std::vector<double> ww(nq, 1. / nq);

        quadrature_collection.push_back(Quadrature<dim>(qpoints[i], ww));
      }

    // Now gather all the information we need
    hp::MappingCollection<dim> mapping_collection(mapping);
    hp::FEValues<dim>          fe_v(mapping_collection,
                           dh->get_fe_collection(),
                           quadrature_collection,
                           update_hessians);
    // Now gather all the information we need
    for (unsigned int i = 0; i < n_cells; ++i)
      {
        AssertThrow(!cells[i]->is_artificial(),
                    VectorTools::ExcPointNotAvailableHere());

        fe_v.reinit(cells[i], i, 0);

        const unsigned int nq = qpoints[i].size();
        std::vector<Vector<typename VectorType::value_type>> vvalues(
          nq, Vector<typename VectorType::value_type>(this->n_components));
        fe_v.get_present_fe_values().get_function_laplacians(data_vector,
                                                             vvalues);

        for (unsigned int q = 0; q < nq; ++q)
          values[maps[i][q]] = vvalues[q];
      }
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::laplacian_list(
    const std::vector<Point<dim>> &               points,
    std::vector<typename VectorType::value_type> &values,
    const unsigned int                            component) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    // Simply forward everything to the vector_gradient_list()
    // function. This requires a temporary object, but everything we
    // do here is so expensive that that really doesn't make any
    // difference any more.
    std::vector<Vector<typename VectorType::value_type>> vvalues(
      points.size(),
      Vector<typename VectorType::value_type>(this->n_components));

    vector_laplacian_list(points, vvalues);

    for (unsigned int q = 0; q < points.size(); ++q)
      values[q] = vvalues[q](component);
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  unsigned int
  FEFieldFunction<dim, DoFHandlerType, VectorType>::compute_point_locations(
    const std::vector<Point<dim>> &                             points,
    std::vector<typename DoFHandlerType::active_cell_iterator> &cells,
    std::vector<std::vector<Point<dim>>> &                      qpoints,
    std::vector<std::vector<unsigned int>> &                    maps) const
  {
    // Calling the GridTools routine and preparing output
    auto cell_qpoint_map =
      GridTools::compute_point_locations(cache, points, cell_hint.get());
    const auto &tria_cells = std::get<0>(cell_qpoint_map);
    cells.resize(tria_cells.size());
    unsigned int i = 0;
    for (const auto &c : tria_cells)
      cells[i++] = typename DoFHandlerType::cell_iterator(*c, dh);
    qpoints = std::get<1>(cell_qpoint_map);
    maps    = std::get<2>(cell_qpoint_map);
    return cells.size();
  }


  template <int dim, typename DoFHandlerType, typename VectorType>
  std_cxx17::optional<Point<dim>>
  FEFieldFunction<dim, DoFHandlerType, VectorType>::get_reference_coordinates(
    const typename DoFHandlerType::active_cell_iterator &cell,
    const Point<dim> &                                   point) const
  {
    try
      {
        Point<dim> qp = mapping.transform_real_to_unit_cell(cell, point);
        if (GeometryInfo<dim>::is_inside_unit_cell(qp))
          return qp;
        else
          return std_cxx17::optional<Point<dim>>();
      }
    catch (const typename Mapping<dim>::ExcTransformationFailed &)
      {
        // transformation failed, so
        // assume the point is
        // outside
        return std_cxx17::optional<Point<dim>>();
      }
  }

} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
