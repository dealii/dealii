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

#ifndef dealii_distributed_active_fe_indices_transfer_h
#define dealii_distributed_active_fe_indices_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/hp/dof_handler.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    /**
     * This class transfers each cell's `active_fe_index` of a hp::FECollection
     * attached to a hp::DoFHandler while refining and/or coarsening a
     * distributed grid and handles the necessary communication.
     *
     * This class therefore does for the `active_fe_index` information of each
     * cell what parallel::distributed::SolutionTransfer does for the values
     * of degrees of freedom defined on a parallel::distributed::Triangulation.
     *
     * If refinement is involved in the data transfer process, the children of
     * a refined cell inherit the `active_fe_index` from their parent. If
     * cells get coarsened into one, the latter will get the least dominating
     * `active_fe_index` amongst its children, as determined by the function
     * hp::FECollection::find_least_face_dominating_fe_in_collection().
     *
     * @note If you use more than one object to attach data to a
     * parallel::distributed::Triangulation at the same time (e.g. a
     * parallel::distributed::SolutionTransfer object), the calls to
     * parallel::distributed::ActiveFEIndicesTransfer::prepare_for_transfer(),
     * parallel::distributed::SolutionTransfer::prepare_for_coarsening_and_refinement()
     * and parallel::distributed::SolutionTransfer::prepare_serialization(),
     * as well as parallel::distributed::ActiveFEIndicesTransfer::unpack() and
     * parallel::distributed::SolutionTransfer::interpolate(), or
     * parallel::distributed::ActiveFEIndicesTransfer::deserialize() and
     * parallel::distributed::SolutionTransfer::deserialize() for serialization,
     * need to be in the same order, respectively.
     *
     * <h3>Transferring each cell's active_fe_index</h3>
     *
     * The following code snippet demonstrates how to transfer all active FE
     * indices across refinement/coarsening of the triangulation that is
     * registered to the hp::DoFHandler member object of this class.
     *
     * @code
     * parallel::distributed::ActiveFEIndicesTransfer<dim, spacedim>
     *   feidx_trans(hp_dof_handler);
     * // flag some cells for refinement and coarsening, e.g.
     * GridRefinement::refine_and_coarsen_fixed_fraction(tria,
     *                                                   error_indicators,
     *                                                   0.3,
     *                                                   0.05);
     *
     * // prepare the triangulation,
     * tria.prepare_coarsening_and_refinement();
     *
     * // prepare the SolutionTransfer object for coarsening and refinement
     * // and give the solution vector that we intend to interpolate later,
     * feidx_trans.prepare_for_transfer();
     *
     * // actually execute the refinement,
     * tria.execute_coarsening_and_refinement();
     *
     * // unpack all active fe indices,
     * feidx_trans.unpack();
     *
     * // redistribute dofs for further use,
     * // using the active FE indices just restored
     * hp_dof_handler.distribute_dofs(fe_collection);
     * @endcode
     *
     *
     * <h3>Use for serialization</h3>
     *
     * This class can be used to serialize and later deserialize a distributed
     * mesh with attached data to separate files. If you use more than one
     * hp::DoFHandler and therefore more than one
     * parallel::distributed::ActiveFEIndicesTransfer object, they need to be
     * serialized and deserialized in the same order.
     *
     * For serialization, the following code snippet saves not only the
     * triangulation itself, but also the active FE indices of a hp::DoFHandler
     * that works on this triangulation:
     * @code
     * parallel::distributed::ActiveFEIndicesTransfer<dim, spacedim>
     *   feidx_trans(hp_dof_handler);
     * feidx_trans.prepare_for_transfer();
     *
     * triangulation.save(filename);
     * @endcode
     *
     * Later, during deserialization, both the triangulation and all active FE
     * indices of the hp::DoFHandler can be restored as follows:
     * @code
     * //[create coarse mesh...]
     * triangulation.load(filename);
     *
     * parallel::distributed::ActiveFEIndicesTransfer<dim, spacedim>
     *   feidx_trans(hp_dof_handler);
     * feidx_trans.deserialize();
     *
     * // distribute dofs for further use,
     * // using the active FE indices just restored
     * hp_dof_handler.distribute_dofs(fe_collection);
     * @endcode
     *
     * @note See documentation of parallel::distributed::SolutionTransfer for
     * matching code snippets in both cases.
     *
     * @ingroup distributed
     * @author Marc Fehling, 2018
     */
    template <int dim, int spacedim = dim>
    class ActiveFEIndicesTransfer
    {
    public:
      /**
       * Constructor.
       *
       * @param[in] dof The hp::DoFHandler on which all
       *   operations will happen. At the time when this constructor
       *   is called, the hp::DoFHandler still points to the triangulation
       *   before the refinement in question happens.
       */
      ActiveFEIndicesTransfer(const hp::DoFHandler<dim, spacedim> &dof_handler);

      /**
       * Prepare the current object for coarsening and refinement or
       * serialization.
       */
      void
      prepare_for_transfer();

      /**
       * Unpack the information previously stored in this object before
       * the mesh was refined or coarsened onto the current set of cells.
       */
      void
      unpack();

      /**
       * Execute the deserialization of the stored information.
       * This needs to be done after calling Triangulation::load().
       */
      void
      deserialize();

    private:
      /**
       * Pointer to the hp::DoFHandler to work with.
       */
      SmartPointer<const hp::DoFHandler<dim, spacedim>,
                   ActiveFEIndicesTransfer<dim, spacedim>>
        dof_handler;

      /**
       * The handle that the parallel::distributed::Triangulation has
       * assigned to this object with which we can access our memory
       * offset and our pack function.
       */
      unsigned int handle;

      /**
       * A callback function used to pack the data on the current mesh into
       * objects that can later be retrieved after refinement, coarsening and
       * repartitioning.
       */
      std::vector<char>
      pack_callback(
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const typename Triangulation<dim, spacedim>::CellStatus     status);

      /**
       * A callback function used to unpack the data on the current mesh that
       * has been packed up previously on the mesh before refinement,
       * coarsening and repartitioning.
       */
      void
      unpack_callback(
        const typename Triangulation<dim, spacedim>::cell_iterator &cell,
        const typename Triangulation<dim, spacedim>::CellStatus     status,
        const boost::iterator_range<std::vector<char>::const_iterator>
          &data_range);
    };
  } // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
