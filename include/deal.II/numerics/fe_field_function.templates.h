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
  template <int dim, typename VectorType, int spacedim>
  FEFieldFunction<dim, VectorType, spacedim>::FEFieldFunction(
    const DoFHandler<dim, spacedim> &mydh,
    const VectorType &               myv,
    const Mapping<dim> &             mymapping)
    : Function<dim, typename VectorType::value_type>(
        mydh.get_fe(0).n_components())
    , dh(&mydh, "FEFieldFunction")
    , data_vector(myv)
    , mapping(mymapping)
    , cache(dh->get_triangulation(), mymapping)
    , cell_hint(dh->end())
  {}



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::set_active_cell(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &newcell)
  {
    cell_hint.get() = newcell;
  }



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::vector_value(
    const Point<dim> &                       p,
    Vector<typename VectorType::value_type> &values) const
  {
    Assert(values.size() == this->n_components,
           ExcDimensionMismatch(values.size(), this->n_components));
    typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
      cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    std_cxx17::optional<Point<dim>> qp = get_reference_coordinates(cell, p);
    if (!qp)
      {
        const std::pair<
          typename dealii::internal::
            ActiveCellIterator<dim, dim, DoFHandler<dim, spacedim>>::type,
          Point<dim>>
          my_pair = GridTools::find_active_cell_around_point(mapping, *dh, p);
        AssertThrow(my_pair.first.state() == IteratorState::valid &&
                      !my_pair.first->is_artificial(),
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



  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  FEFieldFunction<dim, VectorType, spacedim>::value(
    const Point<dim> & p,
    const unsigned int comp) const
  {
    Vector<typename VectorType::value_type> values(this->n_components);
    vector_value(p, values);
    return values(comp);
  }



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::vector_gradient(
    const Point<dim> &                                            p,
    std::vector<Tensor<1, dim, typename VectorType::value_type>> &gradients)
    const
  {
    using number = typename VectorType::value_type;
    Assert(gradients.size() == this->n_components,
           ExcDimensionMismatch(gradients.size(), this->n_components));
    typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
      cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    std_cxx17::optional<Point<dim>> qp = get_reference_coordinates(cell, p);
    if (!qp)
      {
        const std::pair<
          typename dealii::internal::
            ActiveCellIterator<dim, dim, DoFHandler<dim, spacedim>>::type,
          Point<dim>>
          my_pair = GridTools::find_active_cell_around_point(mapping, *dh, p);
        AssertThrow(my_pair.first.state() == IteratorState::valid &&
                      !my_pair.first->is_artificial(),
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



  template <int dim, typename VectorType, int spacedim>
  Tensor<1, dim, typename VectorType::value_type>
  FEFieldFunction<dim, VectorType, spacedim>::gradient(
    const Point<dim> & p,
    const unsigned int comp) const
  {
    std::vector<Tensor<1, dim, typename VectorType::value_type>> grads(
      this->n_components);
    vector_gradient(p, grads);
    return grads[comp];
  }



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::vector_laplacian(
    const Point<dim> &                       p,
    Vector<typename VectorType::value_type> &values) const
  {
    Assert(values.size() == this->n_components,
           ExcDimensionMismatch(values.size(), this->n_components));
    typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
      cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    std_cxx17::optional<Point<dim>> qp = get_reference_coordinates(cell, p);
    if (!qp)
      {
        const std::pair<
          typename dealii::internal::
            ActiveCellIterator<dim, dim, DoFHandler<dim, spacedim>>::type,
          Point<dim>>
          my_pair = GridTools::find_active_cell_around_point(mapping, *dh, p);
        AssertThrow(my_pair.first.state() == IteratorState::valid &&
                      !my_pair.first->is_artificial(),
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



  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  FEFieldFunction<dim, VectorType, spacedim>::laplacian(
    const Point<dim> & p,
    const unsigned int comp) const
  {
    Vector<typename VectorType::value_type> lap(this->n_components);
    vector_laplacian(p, lap);
    return lap[comp];
  }


  // Now the list versions
  // ==============================

  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::vector_value_list(
    const std::vector<Point<dim>> &                       points,
    std::vector<Vector<typename VectorType::value_type>> &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator> cells;
    std::vector<std::vector<Point<dim>>>   qpoints;
    std::vector<std::vector<unsigned int>> maps;

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



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::value_list(
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



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::vector_gradient_list(
    const std::vector<Point<dim>> &points,
    std::vector<std::vector<Tensor<1, dim, typename VectorType::value_type>>>
      &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator> cells;
    std::vector<std::vector<Point<dim>>>   qpoints;
    std::vector<std::vector<unsigned int>> maps;

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
            for (unsigned int l = 0; l < s; ++l)
              values[maps[i][q]][l] = vgrads[q][l];
          }
      }
  }



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::gradient_list(
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



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::vector_laplacian_list(
    const std::vector<Point<dim>> &                       points,
    std::vector<Vector<typename VectorType::value_type>> &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator> cells;
    std::vector<std::vector<Point<dim>>>   qpoints;
    std::vector<std::vector<unsigned int>> maps;

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



  template <int dim, typename VectorType, int spacedim>
  void
  FEFieldFunction<dim, VectorType, spacedim>::laplacian_list(
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



  template <int dim, typename VectorType, int spacedim>
  unsigned int
  FEFieldFunction<dim, VectorType, spacedim>::compute_point_locations(
    const std::vector<Point<dim>> &points,
    std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      &                                     cells,
    std::vector<std::vector<Point<dim>>> &  qpoints,
    std::vector<std::vector<unsigned int>> &maps) const
  {
    // Calling the GridTools routine and preparing output
    auto cell_qpoint_map =
      GridTools::compute_point_locations_try_all(cache,
                                                 points,
                                                 cell_hint.get());
    const auto &tria_cells = std::get<0>(cell_qpoint_map);
    cells.resize(tria_cells.size());
    unsigned int i = 0;
    for (const auto &c : tria_cells)
      cells[i++] = typename DoFHandler<dim, spacedim>::cell_iterator(*c, dh);
    qpoints = std::get<1>(cell_qpoint_map);
    maps    = std::get<2>(cell_qpoint_map);

    // Ensure that we found all points
    AssertThrow(std::get<3>(cell_qpoint_map).empty(),
                VectorTools::ExcPointNotAvailableHere());
    return cells.size();
  }


  template <int dim, typename VectorType, int spacedim>
  std_cxx17::optional<Point<dim>>
  FEFieldFunction<dim, VectorType, spacedim>::get_reference_coordinates(
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const Point<dim> &                                              point) const
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



  template <int dim, class VECTOR>
  RefSpaceFEFieldFunction<dim, VECTOR>::RefSpaceFEFieldFunction(
    const DoFHandler<dim> &dof_handler,
    const VECTOR &         dof_values)
    : Function<dim>(dof_handler.get_fe_collection().n_components())
    , dof_handler(&dof_handler)
    , global_dof_values(&dof_values)
  {
    Assert(dof_handler.n_dofs() == dof_values.size(),
           ExcDimensionMismatch(dof_handler.n_dofs(), dof_values.size()));
  }



  template <int dim, class VECTOR>
  void
  RefSpaceFEFieldFunction<dim, VECTOR>::set_active_cell(
    const typename Triangulation<dim>::active_cell_iterator &cell)
  {
    Assert(
      &cell->get_triangulation() == &dof_handler->get_triangulation(),
      ExcMessage(
        "The incoming cell must belong to the triangulation associated with "
        "the DoFHandler passed to the constructor."));

    const typename DoFHandler<dim>::active_cell_iterator dof_handler_cell(
      &dof_handler->get_triangulation(),
      cell->level(),
      cell->index(),
      dof_handler);

    // Save the element and the local dof values, since this is what we need to
    // evaluate the function.
    element = &dof_handler_cell->get_fe();

    local_dof_indices.resize(element->dofs_per_cell);
    dof_handler_cell->get_dof_indices(local_dof_indices);

    local_dof_values.resize(element->dofs_per_cell);

    for (unsigned int i = 0; i < local_dof_indices.size(); i++)
      local_dof_values[i] =
        dealii::internal::ElementAccess<VECTOR>::get(*global_dof_values,
                                                     local_dof_indices[i]);
  }



  template <int dim, class VECTOR>
  void
  RefSpaceFEFieldFunction<dim, VECTOR>::assert_cell_is_set() const
  {
    // If set cell hasn't been called the size of local_dof_values will be zero.
    Assert(
      local_dof_values.size() > 0,
      ExcMessage(
        "The set_active_cell function has to be called before calling this function."))
  }



  template <int dim, class VECTOR>
  double
  RefSpaceFEFieldFunction<dim, VECTOR>::value(
    const Point<dim> & point,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    assert_cell_is_set();

    double value = 0;
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      value += local_dof_values[i] *
               element->shape_value_component(i, point, component);

    return value;
  }



  template <int dim, class VECTOR>
  Tensor<1, dim>
  RefSpaceFEFieldFunction<dim, VECTOR>::gradient(
    const Point<dim> & point,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    assert_cell_is_set();

    Tensor<1, dim> gradient;
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      gradient += local_dof_values[i] *
                  element->shape_grad_component(i, point, component);

    return gradient;
  }



  template <int dim, class VECTOR>
  SymmetricTensor<2, dim>
  RefSpaceFEFieldFunction<dim, VECTOR>::hessian(
    const Point<dim> & point,
    const unsigned int component) const
  {
    AssertIndexRange(component, this->n_components);
    assert_cell_is_set();

    Tensor<2, dim> hessian;
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      hessian += local_dof_values[i] *
                 element->shape_grad_grad_component(i, point, component);

    return symmetrize(hessian);
  }

} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
