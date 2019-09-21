// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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


#include <deal.II/distributed/cell_weights.h>

#include <deal.II/dofs/dof_accessor.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  template <int dim, int spacedim>
  CellWeights<dim, spacedim>::CellWeights(
    const hp::DoFHandler<dim, spacedim> &dof_handler)
    : dof_handler(&dof_handler, typeid(*this).name())
  {
    triangulation = (dynamic_cast<parallel::TriangulationBase<dim, spacedim> *>(
      const_cast<dealii::Triangulation<dim, spacedim> *>(
        &(this->dof_handler->get_triangulation()))));

    if (triangulation != nullptr)
      {
        // Choose default weighting.
        register_constant_weighting();

        // Add callback function to the cell_weight signal and store its
        // connection.
        tria_listener = triangulation->signals.cell_weight.connect(
          [this](
            const typename Triangulation<dim, spacedim>::cell_iterator &cell_,
            const typename Triangulation<dim, spacedim>::CellStatus status) {
            return this->weight_callback(cell_, status);
          });
      }
    else
      Assert(
        triangulation != nullptr,
        ExcMessage(
          "parallel::CellWeights requires a parallel::TriangulationBase object."));
  }


  template <int dim, int spacedim>
  CellWeights<dim, spacedim>::~CellWeights()
  {
    tria_listener.disconnect();
  }



  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::register_constant_weighting(
    const unsigned int factor)
  {
    weighting_function =
      [factor](const FiniteElement<dim, spacedim> &,
               const typename hp::DoFHandler<dim, spacedim>::cell_iterator &)
      -> unsigned int { return factor; };
  }


  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::register_ndofs_weighting(
    const unsigned int factor)
  {
    weighting_function =
      [factor](const FiniteElement<dim, spacedim> &active_fe,
               const typename hp::DoFHandler<dim, spacedim>::cell_iterator &)
      -> unsigned int { return factor * active_fe.dofs_per_cell; };
  }


  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::register_ndofs_squared_weighting(
    const unsigned int factor)
  {
    weighting_function =
      [factor](const FiniteElement<dim, spacedim> &active_fe,
               const typename hp::DoFHandler<dim, spacedim>::cell_iterator &)
      -> unsigned int {
      return factor * active_fe.dofs_per_cell * active_fe.dofs_per_cell;
    };
  }


  template <int dim, int spacedim>
  void
  CellWeights<dim, spacedim>::register_custom_weighting(
    const std::function<unsigned int(
      const FiniteElement<dim, spacedim> &,
      const typename hp::DoFHandler<dim, spacedim>::cell_iterator &)>
      custom_function)
  {
    weighting_function = custom_function;
  }



  template <int dim, int spacedim>
  unsigned int
  CellWeights<dim, spacedim>::weight_callback(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell_,
    const typename Triangulation<dim, spacedim>::CellStatus     status)
  {
    // Check if we are still working with the correct combination of
    // Triangulation and DoFHandler.
    Assert(&(*triangulation) == &(dof_handler->get_triangulation()),
           ExcMessage(
             "Triangulation associated with the DoFHandler has changed!"));

    // Convert cell type from Triangulation to DoFHandler to be able
    // to access the information about the degrees of freedom.
    const typename hp::DoFHandler<dim, spacedim>::cell_iterator cell(
      *cell_, dof_handler);

    // Determine which FiniteElement object will be present on this cell after
    // refinement and will thus specify the number of degrees of freedom.
    unsigned int fe_index = numbers::invalid_unsigned_int;
    switch (status)
      {
        case Triangulation<dim, spacedim>::CELL_PERSIST:
        case Triangulation<dim, spacedim>::CELL_REFINE:
        case Triangulation<dim, spacedim>::CELL_INVALID:
          fe_index = cell->future_fe_index();
          break;

        case Triangulation<dim, spacedim>::CELL_COARSEN:
          {
            std::set<unsigned int> fe_indices_children;
            for (unsigned int child_index = 0; child_index < cell->n_children();
                 ++child_index)
              {
                const auto child = cell->child(child_index);
                Assert(child->active() && child->coarsen_flag_set(),
                       typename dealii::Triangulation<
                         dim>::ExcInconsistentCoarseningFlags());

                fe_indices_children.insert(child->future_fe_index());
              }
            Assert(!fe_indices_children.empty(), ExcInternalError());

            fe_index =
              dof_handler->get_fe_collection().find_dominated_fe_extended(
                fe_indices_children, /*codim=*/0);

            Assert(fe_index != numbers::invalid_unsigned_int,
                   typename dealii::hp::FECollection<
                     dim>::ExcNoDominatedFiniteElementAmongstChildren());
          }
          break;

        default:
          Assert(false, ExcInternalError());
          break;
      }

    // Return the cell weight determined by the function of choice.
    return weighting_function(dof_handler->get_fe(fe_index), cell);
  }
} // namespace parallel


// explicit instantiations
#include "cell_weights.inst"

DEAL_II_NAMESPACE_CLOSE
