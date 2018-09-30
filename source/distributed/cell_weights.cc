// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


#include <deal.II/base/config.h>

#include <deal.II/distributed/cell_weights.h>

#include <deal.II/dofs/dof_accessor.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    CellWeights<dim, spacedim>::CellWeights(
      const hp::DoFHandler<dim, spacedim> &dof_handler)
      : dof_handler(&dof_handler, typeid(*this).name())
    {
      triangulation =
        (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
          const_cast<dealii::Triangulation<dim, spacedim> *>(
            &(this->dof_handler->get_triangulation()))));

      if (triangulation != nullptr)
        {
          // Choose default weighting.
          register_constant_weighting();

          // Add callback function to the cell_weight signal and store its
          // connection.
          tria_listener = triangulation->signals.cell_weight.connect(
            std::bind(&CellWeights<dim, spacedim>::weight_callback,
                      std::ref(*this),
                      std::placeholders::_1,
                      std::placeholders::_2));
        }
      else
        Assert(
          triangulation != nullptr,
          ExcMessage(
            "parallel::distributed::CellWeights requires a parallel::distributed::Triangulation object."));
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
          case parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_PERSIST:
          case parallel::distributed::Triangulation<dim, spacedim>::CELL_REFINE:
          case parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_INVALID:
            fe_index = cell->active_fe_index();
            break;

          case parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_COARSEN:
            {
              std::set<unsigned int> fe_indices_children;
              for (unsigned int child_index = 0;
                   child_index < GeometryInfo<dim>::max_children_per_cell;
                   ++child_index)
                fe_indices_children.insert(
                  cell->child(child_index)->active_fe_index());

              fe_index = dof_handler->get_fe().find_least_face_dominating_fe(
                fe_indices_children);

              Assert(fe_index != numbers::invalid_unsigned_int,
                     ExcMessage(
                       "No FiniteElement has been found in your FECollection "
                       "that dominates all children of a cell you are trying "
                       "to coarsen!"));
            }
            break;

          default:
            Assert(false, ExcInternalError());
            break;
        }

      // Return the cell weight determined by the function of choice.
      return weighting_function(dof_handler->get_fe(fe_index), cell);
    }
  } // namespace distributed
} // namespace parallel


// explicit instantiations
#include "cell_weights.inst"

DEAL_II_NAMESPACE_CLOSE
