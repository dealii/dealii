// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/error_estimator.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

DEAL_II_NAMESPACE_OPEN


template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const Mapping<1, spacedim>    &mapping,
  const DoFHandler<1, spacedim> &dof_handler,
  const Quadrature<0>           &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                           &neumann_bc,
  const ReadVector<Number> &solution,
  Vector<float>            &error,
  const ComponentMask      &component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  // just pass on to the other function
  std::vector<const ReadVector<Number> *> solutions(1, &solution);
  std::vector<Vector<float> *>            errors(1, &error);
  ArrayView<Vector<float> *>              error_view = make_array_view(errors);
  estimate(mapping,
           dof_handler,
           quadrature,
           neumann_bc,
           make_array_view(solutions),
           error_view,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const DoFHandler<1, spacedim> &dof_handler,
  const Quadrature<0>           &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                           &neumann_bc,
  const ReadVector<Number> &solution,
  Vector<float>            &error,
  const ComponentMask      &component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  const auto reference_cell = ReferenceCells::Line;
  estimate(reference_cell.template get_default_linear_mapping<1, spacedim>(),
           dof_handler,
           quadrature,
           neumann_bc,
           solution,
           error,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const DoFHandler<1, spacedim> &dof_handler,
  const Quadrature<0>           &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &neumann_bc,
  const ArrayView<const ReadVector<Number> *> &solutions,
  ArrayView<Vector<float> *>                  &errors,
  const ComponentMask                         &component_mask,
  const Function<spacedim>                    *coefficients,
  const unsigned int                           n_threads,
  const types::subdomain_id                    subdomain_id,
  const types::material_id                     material_id,
  const Strategy                               strategy)
{
  const auto reference_cell = ReferenceCells::Line;
  estimate(reference_cell.template get_default_linear_mapping<1, spacedim>(),
           dof_handler,
           quadrature,
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const hp::MappingCollection<1, spacedim> &mapping,
  const DoFHandler<1, spacedim>            &dof_handler,
  const hp::QCollection<0>                 &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                           &neumann_bc,
  const ReadVector<Number> &solution,
  Vector<float>            &error,
  const ComponentMask      &component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  // just pass on to the other function
  std::vector<const ReadVector<Number> *> solutions(1, &solution);
  std::vector<Vector<float> *>            errors(1, &error);
  ArrayView<Vector<float> *>              error_view = make_array_view(errors);
  estimate(mapping,
           dof_handler,
           quadrature,
           neumann_bc,
           make_array_view(solutions),
           error_view,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const Mapping<1, spacedim>    &mapping,
  const DoFHandler<1, spacedim> &dof_handler,
  const hp::QCollection<0>      &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                           &neumann_bc,
  const ReadVector<Number> &solution,
  Vector<float>            &error,
  const ComponentMask      &component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  // DEPRECATED
  // just pass on to the other function
  std::vector<const ReadVector<Number> *>  solutions(1, &solution);
  std::vector<Vector<float> *>             errors(1, &error);
  ArrayView<Vector<float> *>               error_view = make_array_view(errors);
  const hp::MappingCollection<1, spacedim> mapping_collection(mapping);
  estimate(mapping_collection,
           dof_handler,
           quadrature,
           neumann_bc,
           make_array_view(solutions),
           error_view,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const DoFHandler<1, spacedim> &dof_handler,
  const hp::QCollection<0>      &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                           &neumann_bc,
  const ReadVector<Number> &solution,
  Vector<float>            &error,
  const ComponentMask      &component_mask,
  const Function<spacedim> *coefficients,
  const unsigned int        n_threads,
  const types::subdomain_id subdomain_id,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  const auto reference_cell = ReferenceCells::Line;
  const hp::MappingCollection<1, spacedim> mapping(
    reference_cell.template get_default_linear_mapping<1, spacedim>());
  estimate(mapping,
           dof_handler,
           quadrature,
           neumann_bc,
           solution,
           error,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const DoFHandler<1, spacedim> &dof_handler,
  const hp::QCollection<0>      &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &neumann_bc,
  const ArrayView<const ReadVector<Number> *> &solutions,
  ArrayView<Vector<float> *>                  &errors,
  const ComponentMask                         &component_mask,
  const Function<spacedim>                    *coefficients,
  const unsigned int                           n_threads,
  const types::subdomain_id                    subdomain_id,
  const types::material_id                     material_id,
  const Strategy                               strategy)
{
  const auto reference_cell = ReferenceCells::Line;
  const hp::MappingCollection<1, spacedim> mapping(
    reference_cell.template get_default_linear_mapping<1, spacedim>());
  estimate(mapping,
           dof_handler,
           quadrature,
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const hp::MappingCollection<1, spacedim> &mapping,
  const DoFHandler<1, spacedim>            &dof_handler,
  const hp::QCollection<0> &,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &neumann_bc,
  const ArrayView<const ReadVector<Number> *> &solutions,
  ArrayView<Vector<float> *>                  &errors,
  const ComponentMask                         &component_mask,
  const Function<spacedim>                    *coefficient,
  const unsigned int,
  const types::subdomain_id subdomain_id_,
  const types::material_id  material_id,
  const Strategy            strategy)
{
  AssertThrow(strategy == cell_diameter_over_24, ExcNotImplemented());
  using number                     = Number;
  types::subdomain_id subdomain_id = numbers::invalid_subdomain_id;
  if (const auto *triangulation = dynamic_cast<
        const parallel::DistributedTriangulationBase<1, spacedim> *>(
        &dof_handler.get_triangulation()))
    {
      Assert((subdomain_id_ == numbers::invalid_subdomain_id) ||
               (subdomain_id_ == triangulation->locally_owned_subdomain()),
             ExcMessage(
               "For distributed Triangulation objects and associated "
               "DoFHandler objects, asking for any subdomain other than the "
               "locally owned one does not make sense."));
      subdomain_id = triangulation->locally_owned_subdomain();
    }
  else
    {
      subdomain_id = subdomain_id_;
    }

  const unsigned int n_components       = dof_handler.get_fe(0).n_components();
  const unsigned int n_solution_vectors = solutions.size();

  // sanity checks
  Assert(neumann_bc.find(numbers::internal_face_boundary_id) ==
           neumann_bc.end(),
         ExcMessage("You are not allowed to list the special boundary "
                    "indicator for internal boundaries in your boundary "
                    "value map."));

  for (const auto &boundary_function : neumann_bc)
    {
      (void)boundary_function;
      Assert(boundary_function.second->n_components == n_components,
             ExcInvalidBoundaryFunction(boundary_function.first,
                                        boundary_function.second->n_components,
                                        n_components));
    }

  Assert(component_mask.represents_n_components(n_components),
         ExcInvalidComponentMask());
  Assert(component_mask.n_selected_components(n_components) > 0,
         ExcInvalidComponentMask());

  Assert((coefficient == nullptr) ||
           (coefficient->n_components == n_components) ||
           (coefficient->n_components == 1),
         ExcInvalidCoefficient());

  Assert(solutions.size() > 0, ExcNoSolutions());
  Assert(solutions.size() == errors.size(),
         ExcIncompatibleNumberOfElements(solutions.size(), errors.size()));
  for (unsigned int n = 0; n < solutions.size(); ++n)
    Assert(solutions[n]->size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(solutions[n]->size(), dof_handler.n_dofs()));

  Assert((coefficient == nullptr) ||
           (coefficient->n_components == n_components) ||
           (coefficient->n_components == 1),
         ExcInvalidCoefficient());

  for (const auto &boundary_function : neumann_bc)
    {
      (void)boundary_function;
      Assert(boundary_function.second->n_components == n_components,
             ExcInvalidBoundaryFunction(boundary_function.first,
                                        boundary_function.second->n_components,
                                        n_components));
    }

  // reserve one slot for each cell and set it to zero
  for (unsigned int n = 0; n < n_solution_vectors; ++n)
    (*errors[n]).reinit(dof_handler.get_triangulation().n_active_cells());

  // fields to get the gradients on the present and the neighbor cell.
  //
  // for the neighbor gradient, we need several auxiliary fields, depending on
  // the way we get it (see below)
  std::vector<std::vector<std::vector<Tensor<1, spacedim, number>>>>
    gradients_here(n_solution_vectors,
                   std::vector<std::vector<Tensor<1, spacedim, number>>>(
                     2,
                     std::vector<Tensor<1, spacedim, number>>(n_components)));
  std::vector<std::vector<std::vector<Tensor<1, spacedim, number>>>>
    gradients_neighbor(gradients_here);
  std::vector<Vector<typename ProductType<number, double>::type>>
    grad_dot_n_neighbor(n_solution_vectors,
                        Vector<typename ProductType<number, double>::type>(
                          n_components));

  // reserve some space for coefficient values at one point.  if there is no
  // coefficient, then we fill it by unity once and for all and don't set it
  // any more
  Vector<double> coefficient_values(n_components);
  if (coefficient == nullptr)
    for (unsigned int c = 0; c < n_components; ++c)
      coefficient_values(c) = 1;

  const QTrapezoid<1>      quadrature;
  const hp::QCollection<1> q_collection(quadrature);
  const QGauss<0>          face_quadrature(1);
  const hp::QCollection<0> q_face_collection(face_quadrature);

  const hp::FECollection<1, spacedim> &fe = dof_handler.get_fe_collection();


  hp::FEValues<1, spacedim>     fe_values(mapping,
                                      fe,
                                      q_collection,
                                      update_gradients);
  hp::FEFaceValues<1, spacedim> fe_face_values(
    /*mapping,*/ fe, q_face_collection, update_normal_vectors);

  // loop over all cells and do something on the cells which we're told to
  // work on. note that the error indicator is only a sum over the two
  // contributions from the two vertices of each cell.
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (((subdomain_id == numbers::invalid_subdomain_id) ||
         (cell->subdomain_id() == subdomain_id)) &&
        ((material_id == numbers::invalid_material_id) ||
         (cell->material_id() == material_id)))
      {
        for (unsigned int n = 0; n < n_solution_vectors; ++n)
          (*errors[n])(cell->active_cell_index()) = 0;

        fe_values.reinit(cell);
        for (unsigned int s = 0; s < n_solution_vectors; ++s)
          fe_values.get_present_fe_values().get_function_gradients(
            *solutions[s], gradients_here[s]);

        // loop over the two points bounding this line. n==0 is left point,
        // n==1 is right point
        for (unsigned int n = 0; n < 2; ++n)
          {
            // find left or right active neighbor
            auto neighbor = cell->neighbor(n);
            if (neighbor.state() == IteratorState::valid)
              while (neighbor->has_children())
                neighbor = neighbor->child(n == 0 ? 1 : 0);

            fe_face_values.reinit(cell, n);
            Tensor<1, spacedim> normal =
              fe_face_values.get_present_fe_values().get_normal_vectors()[0];

            if (neighbor.state() == IteratorState::valid)
              {
                fe_values.reinit(neighbor);

                for (unsigned int s = 0; s < n_solution_vectors; ++s)
                  fe_values.get_present_fe_values().get_function_gradients(
                    *solutions[s], gradients_neighbor[s]);

                fe_face_values.reinit(neighbor, n == 0 ? 1 : 0);
                Tensor<1, spacedim> neighbor_normal =
                  fe_face_values.get_present_fe_values()
                    .get_normal_vectors()[0];

                // extract the gradient in normal direction of all the
                // components.
                for (unsigned int s = 0; s < n_solution_vectors; ++s)
                  for (unsigned int c = 0; c < n_components; ++c)
                    grad_dot_n_neighbor[s](c) =
                      -(gradients_neighbor[s][n == 0 ? 1 : 0][c] *
                        neighbor_normal);
              }
            else if (neumann_bc.find(n) != neumann_bc.end())
              // if Neumann b.c., then fill the gradients field which will be
              // used later on.
              {
                if (n_components == 1)
                  {
                    const Number v =
                      neumann_bc.find(n)->second->value(cell->vertex(n));

                    for (unsigned int s = 0; s < n_solution_vectors; ++s)
                      grad_dot_n_neighbor[s](0) = v;
                  }
                else
                  {
                    Vector<Number> v(n_components);
                    neumann_bc.find(n)->second->vector_value(cell->vertex(n),
                                                             v);

                    for (unsigned int s = 0; s < n_solution_vectors; ++s)
                      grad_dot_n_neighbor[s] = v;
                  }
              }
            else
              // fill with zeroes.
              for (unsigned int s = 0; s < n_solution_vectors; ++s)
                grad_dot_n_neighbor[s] = 0;

            // if there is a coefficient, then evaluate it at the present
            // position. if there is none, reuse the preset values.
            if (coefficient != nullptr)
              {
                if (coefficient->n_components == 1)
                  {
                    const double c_value = coefficient->value(cell->vertex(n));
                    for (unsigned int c = 0; c < n_components; ++c)
                      coefficient_values(c) = c_value;
                  }
                else
                  coefficient->vector_value(cell->vertex(n),
                                            coefficient_values);
              }


            for (unsigned int s = 0; s < n_solution_vectors; ++s)
              for (unsigned int component = 0; component < n_components;
                   ++component)
                if (component_mask[component] == true)
                  {
                    // get gradient here
                    const typename ProductType<number, double>::type
                      grad_dot_n_here =
                        gradients_here[s][n][component] * normal;

                    const typename ProductType<number, double>::type jump =
                      ((grad_dot_n_here - grad_dot_n_neighbor[s](component)) *
                       coefficient_values(component));
                    (*errors[s])(cell->active_cell_index()) +=
                      numbers::NumberTraits<
                        typename ProductType<number,
                                             double>::type>::abs_square(jump) *
                      cell->diameter();
                  }
          }

        for (unsigned int s = 0; s < n_solution_vectors; ++s)
          (*errors[s])(cell->active_cell_index()) =
            std::sqrt((*errors[s])(cell->active_cell_index()));
      }
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const Mapping<1, spacedim>    &mapping,
  const DoFHandler<1, spacedim> &dof_handler,
  const hp::QCollection<0>      &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &neumann_bc,
  const ArrayView<const ReadVector<Number> *> &solutions,
  ArrayView<Vector<float> *>                  &errors,
  const ComponentMask                         &component_mask,
  const Function<spacedim>                    *coefficients,
  const unsigned int                           n_threads,
  const types::subdomain_id                    subdomain_id,
  const types::material_id                     material_id,
  const Strategy                               strategy)
{
  // DEPRECATED
  // just pass on to the other function
  const hp::MappingCollection<1, spacedim> mapping_collection(mapping);
  estimate(mapping_collection,
           dof_handler,
           quadrature,
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



template <int spacedim>
template <typename Number>
void
KellyErrorEstimator<1, spacedim>::estimate(
  const Mapping<1, spacedim>    &mapping,
  const DoFHandler<1, spacedim> &dof_handler,
  const Quadrature<0>           &quadrature,
  const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &neumann_bc,
  const ArrayView<const ReadVector<Number> *> &solutions,
  ArrayView<Vector<float> *>                  &errors,
  const ComponentMask                         &component_mask,
  const Function<spacedim>                    *coefficients,
  const unsigned int                           n_threads,
  const types::subdomain_id                    subdomain_id,
  const types::material_id                     material_id,
  const Strategy                               strategy)
{
  const hp::MappingCollection<1, spacedim> mapping_collection(mapping);
  const hp::QCollection<0>                 quadrature_collection(quadrature);
  estimate(mapping_collection,
           dof_handler,
           quadrature_collection,
           neumann_bc,
           solutions,
           errors,
           component_mask,
           coefficients,
           n_threads,
           subdomain_id,
           material_id,
           strategy);
}



// explicit instantiations
#include "numerics/error_estimator_1d.inst"


DEAL_II_NAMESPACE_CLOSE
