// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/distributed/cell_weights.h>

#include <deal.II/dofs/dof_accessor.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN


namespace parallel
{
  template <int dim, int spacedim>
  CellWeights<dim, spacedim>::CellWeights(
    const DoFHandler<dim, spacedim> &dof_handler,
    const WeightingFunction         &weighting_function)
  {
    reinit(dof_handler, weighting_function);
  }



  template <int dim, int spacedim>
  CellWeights<dim, spacedim>::~CellWeights()
  {
    connection.disconnect();
  }



  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::reinit(
    const DoFHandler<dim, spacedim> &dof_handler,
    const typename CellWeights<dim, spacedim>::WeightingFunction
      &weighting_function)
  {
    connection.disconnect();
    connection = dof_handler.get_triangulation().signals.weight.connect(
      make_weighting_callback(dof_handler, weighting_function));
  }



  // ---------- static member functions ----------

  // ---------- selection of weighting functions ----------

  template <int dim, int spacedim>
  typename CellWeights<dim, spacedim>::WeightingFunction
  CellWeights<dim, spacedim>::constant_weighting(const unsigned int factor)
  {
    return [factor](const typename DoFHandler<dim, spacedim>::cell_iterator &,
                    const FiniteElement<dim, spacedim> &) -> unsigned int {
      return factor;
    };
  }



  template <int dim, int spacedim>
  typename CellWeights<dim, spacedim>::WeightingFunction
  CellWeights<dim, spacedim>::ndofs_weighting(
    const std::pair<float, float> &coefficients)
  {
    return [coefficients](
             const typename DoFHandler<dim, spacedim>::cell_iterator &,
             const FiniteElement<dim, spacedim> &future_fe) -> unsigned int {
      const float result =
        std::trunc(coefficients.first *
                   std::pow(future_fe.n_dofs_per_cell(), coefficients.second));

      Assert(result >= 0. &&
               result <=
                 static_cast<float>(std::numeric_limits<unsigned int>::max()),
             ExcMessage(
               "Cannot cast determined weight for this cell to unsigned int!"));

      return static_cast<unsigned int>(result);
    };
  }



  template <int dim, int spacedim>
  typename CellWeights<dim, spacedim>::WeightingFunction
  CellWeights<dim, spacedim>::ndofs_weighting(
    const std::vector<std::pair<float, float>> &coefficients)
  {
    return [coefficients](
             const typename DoFHandler<dim, spacedim>::cell_iterator &,
             const FiniteElement<dim, spacedim> &future_fe) -> unsigned int {
      float result = 0;
      for (const auto &pair : coefficients)
        result +=
          pair.first * std::pow(future_fe.n_dofs_per_cell(), pair.second);
      result = std::trunc(result);

      Assert(result >= 0. &&
               result <=
                 static_cast<float>(std::numeric_limits<unsigned int>::max()),
             ExcMessage(
               "Cannot cast determined weight for this cell to unsigned int!"));

      return static_cast<unsigned int>(result);
    };
  }



  // ---------- handling callback functions ----------

  template <int dim, int spacedim>
  std::function<unsigned int(
    const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellStatus                                                    status)>
  CellWeights<dim, spacedim>::make_weighting_callback(
    const DoFHandler<dim, spacedim> &dof_handler,
    const typename CellWeights<dim, spacedim>::WeightingFunction
      &weighting_function)
  {
    const parallel::TriangulationBase<dim, spacedim> *tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &(dof_handler.get_triangulation()));

    Assert(
      tria != nullptr,
      ExcMessage(
        "parallel::CellWeights requires a parallel::TriangulationBase object."));

    return [&dof_handler, tria, weighting_function](
             const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                             &cell,
             const CellStatus status) -> unsigned int {
      return CellWeights<dim, spacedim>::weighting_callback(
        cell, status, dof_handler, *tria, weighting_function);
    };
  }



  template <int dim, int spacedim>
  unsigned int
  CellWeights<dim, spacedim>::weighting_callback(
    const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell_,
    const CellStatus                                                    status,
    const DoFHandler<dim, spacedim>                  &dof_handler,
    const parallel::TriangulationBase<dim, spacedim> &triangulation,
    const typename CellWeights<dim, spacedim>::WeightingFunction
      &weighting_function)
  {
    // Check if we are still working with the correct combination of
    // Triangulation and DoFHandler.
    Assert(&triangulation == &(dof_handler.get_triangulation()),
           ExcMessage(
             "Triangulation associated with the DoFHandler has changed!"));

    // Skip if the DoFHandler has not been initialized yet.
    if (dof_handler.get_fe_collection().empty())
      return 0;

    // Convert cell type from Triangulation to DoFHandler to be able
    // to access the information about the degrees of freedom.
    const typename DoFHandler<dim, spacedim>::cell_iterator cell(*cell_,
                                                                 &dof_handler);

    // Determine which FiniteElement object will be present on this cell after
    // refinement and will thus specify the number of degrees of freedom.
    unsigned int fe_index = numbers::invalid_unsigned_int;
    switch (status)
      {
        case CellStatus::cell_will_persist:
        case CellStatus::cell_will_be_refined:
        case CellStatus::cell_invalid:
          fe_index = cell->future_fe_index();
          break;

        case CellStatus::children_will_be_coarsened:
          if constexpr (running_in_debug_mode())
            {
              for (const auto &child : cell->child_iterators())
                Assert(child->is_active() && child->coarsen_flag_set(),
                       typename dealii::Triangulation<
                         dim>::ExcInconsistentCoarseningFlags());
            }

          fe_index = dealii::internal::hp::DoFHandlerImplementation::
            dominated_future_fe_on_children<dim, spacedim>(cell);
          break;

        default:
          DEAL_II_ASSERT_UNREACHABLE();
          break;
      }

    // Return the cell weight determined by the function of choice.
    return weighting_function(cell, dof_handler.get_fe(fe_index));
  }
} // namespace parallel


// explicit instantiations
#include "distributed/cell_weights.inst"

DEAL_II_NAMESPACE_CLOSE
