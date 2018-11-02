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

#ifdef DEAL_II_WITH_P4EST

#  include <deal.II/distributed/active_fe_indices_transfer.h>
#  include <deal.II/distributed/tria.h>

#  include <deal.II/dofs/dof_accessor.h>

#  include <deal.II/hp/dof_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    ActiveFEIndicesTransfer<dim, spacedim>::ActiveFEIndicesTransfer(
      const hp::DoFHandler<dim, spacedim> &dof_handler)
      : dof_handler(&dof_handler, typeid(*this).name())
      , handle(numbers::invalid_unsigned_int)
    {
      Assert(
        (dynamic_cast<
           const parallel::distributed::Triangulation<dim, spacedim> *>(
           &(this->dof_handler)->get_triangulation()) != nullptr),
        ExcMessage(
          "parallel::distributed::ActiveFEIndicesTransfer requires a parallel::distributed::Triangulation object."));
    }



    template <int dim, int spacedim>
    void
    ActiveFEIndicesTransfer<dim, spacedim>::prepare_for_transfer()
    {
      parallel::distributed::Triangulation<dim, spacedim> *tria =
        (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
          const_cast<dealii::Triangulation<dim, spacedim> *>(
            &dof_handler->get_triangulation())));
      Assert(tria != nullptr, ExcInternalError());

      handle = tria->register_data_attach(
        std::bind(&ActiveFEIndicesTransfer<dim, spacedim>::pack_callback,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2),
        /*returns_variable_size_data=*/false);
    }



    template <int dim, int spacedim>
    void
    ActiveFEIndicesTransfer<dim, spacedim>::unpack()
    {
      // TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim, spacedim> *tria =
        (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
          const_cast<dealii::Triangulation<dim, spacedim> *>(
            &dof_handler->get_triangulation())));
      Assert(tria != nullptr, ExcInternalError());

      tria->notify_ready_to_unpack(
        handle,
        std::bind(&ActiveFEIndicesTransfer<dim, spacedim>::unpack_callback,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3));
    }



    template <int dim, int spacedim>
    void
    ActiveFEIndicesTransfer<dim, spacedim>::deserialize()
    {
      // For deserialization, we need to register this object
      // to the triangulation first to get a valid handle for
      // data access.
      prepare_for_transfer();
      unpack();
    }



    template <int dim, int spacedim>
    std::vector<char>
    ActiveFEIndicesTransfer<dim, spacedim>::pack_callback(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell_,
      const typename Triangulation<dim, spacedim>::CellStatus     status)
    {
      const typename hp::DoFHandler<dim, spacedim>::cell_iterator cell(
        *cell_, dof_handler);

      unsigned int fe_index = numbers::invalid_unsigned_int;

      switch (status)
        {
          case parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_PERSIST:
          case parallel::distributed::Triangulation<dim, spacedim>::CELL_REFINE:
            fe_index = cell->active_fe_index();
            break;

          case parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_COARSEN:
            // In this case, the callback function will be called on the parent
            // cell which shall store the packed information. We need to choose
            // from its children which ID to store.
            {
              std::set<unsigned int> fe_indices_children;
              for (unsigned int child_index = 0;
                   child_index < GeometryInfo<dim>::max_children_per_cell;
                   ++child_index)
                {
                  typename hp::DoFHandler<dim, spacedim>::cell_iterator child =
                    cell->child(child_index);

                  fe_indices_children.insert(child->active_fe_index());
                }

              fe_index =
                dof_handler->get_fe_collection()
                  .find_least_dominating_fe_in_collection(fe_indices_children,
                                                          /*codim=*/0);

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

      return Utilities::pack(fe_index, /*allow_compression=*/false);
    }



    template <int dim, int spacedim>
    void
    ActiveFEIndicesTransfer<dim, spacedim>::unpack_callback(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell_,
      const typename Triangulation<dim, spacedim>::CellStatus     status,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &data_range)
    {
      typename hp::DoFHandler<dim, spacedim>::cell_iterator cell(*cell_,
                                                                 dof_handler);

      const unsigned int fe_index =
        Utilities::unpack<unsigned int>(data_range.begin(),
                                        data_range.end(),
                                        /*allow_compression=*/false);

      Assert(fe_index <= dof_handler->get_fe().size(), ExcInternalError());

      switch (status)
        {
          case parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_PERSIST:
          case parallel::distributed::Triangulation<dim,
                                                    spacedim>::CELL_COARSEN:
            cell->set_active_fe_index(fe_index);
            break;

          case parallel::distributed::Triangulation<dim, spacedim>::CELL_REFINE:
            // In this case, the callback function will be called on the parent
            // cell which stores the packed information. We need to distribute
            // it on its children.
            for (unsigned int child_index = 0;
                 child_index < GeometryInfo<dim>::max_children_per_cell;
                 ++child_index)
              cell->child(child_index)->set_active_fe_index(fe_index);
            break;

          default:
            Assert(false, ExcInternalError());
            break;
        }
    }
  } // namespace distributed
} // namespace parallel


// explicit instantiations
#  include "active_fe_indices_transfer.inst"

DEAL_II_NAMESPACE_CLOSE

#endif
