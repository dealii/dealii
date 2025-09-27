// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_distributed_cell_data_transfer_h
#define dealii_distributed_cell_data_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/base/observer_pointer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/numerics/adaptation_strategies.h>

#include <boost/range/iterator_range.hpp>

#include <functional>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    /**
     * Transfer data that is associated with each active cell (like error
     * indicators) while refining and/or coarsening a distributed triangulation
     * and handle the necessary communication.
     *
     * This class therefore does for cell-related information what
     * parallel::distributed::SolutionTransfer does for the values of degrees of
     * freedom defined on a parallel::distributed::Triangulation.
     *
     * This class has been designed to operate on any kind of datatype that is
     * serializable. A non-distributed container (like Vector or `std::vector`)
     * has to be provided, which holds the cell-wise data in the same order as
     * active cells are traversed. In other words, each entry corresponds to the
     * cell with the same index CellAccessor::active_cell_index(), and the
     * container has to be of size Triangulation::n_active_cells().
     *
     * <h3>Transferring cell-wise data</h3>
     *
     * The following code snippet demonstrates how to transfer cell-related data
     * across refinement/coarsening of the registered triangulation.
     *
     * @code
     * // prepare the triangulation,
     * triangulation.prepare_coarsening_and_refinement();
     *
     * // prepare the CellDataTransfer object for coarsening and refinement
     * // and give the cell data vector that we intend to unpack later,
     * Vector<double> data_to_transfer(triangulation.n_active_cells());
     * //[fill data_to_transfer with cell-wise values...]
     *
     * parallel::distributed::CellDataTransfer<dim, spacedim, Vector<double>>
     *   cell_data_trans(triangulation);
     * cell_data_trans.prepare_for_coarsening_and_refinement(data_to_transfer);
     *
     * // actually execute the refinement,
     * triangulation.execute_coarsening_and_refinement();
     *
     * // unpack transferred data,
     * Vector<double> transferred_data(triangulation.n_active_cells());
     * cell_data_trans.unpack(transferred_data);
     *
     * @endcode
     *
     *
     * <h3>Use for serialization</h3>
     *
     * This class can be used to serialize and later deserialize a distributed
     * mesh with attached data to separate files.
     *
     * For serialization, the following code snippet saves not only the
     * triangulation itself, but also the cell-wise data attached:
     * @code
     * Vector<double> data_to_transfer(triangulation.n_active_cells());
     * //[fill data_to_transfer with cell-wise values...]
     *
     * parallel::distributed::CellDataTransfer<dim, spacedim, Vector<double>>
     *   cell_data_trans(triangulation);
     * cell_data_trans.prepare_for_serialization(data_to_transfer);
     *
     * triangulation.save(filename);
     * @endcode
     *
     * Later, during deserialization, both triangulation and data can be
     * restored as follows:
     * @code
     * //[create coarse mesh...]
     * triangulation.load(filename);
     *
     * parallel::distributed::CellDataTransfer<dim, spacedim, Vector<double>>
     *   cell_data_trans(triangulation);
     * Vector<double> transferred_data(triangulation.n_active_cells());
     * cell_data_trans.deserialize(transferred_data);
     * @endcode
     *
     * @note If you use more than one object to transfer data via the
     * parallel::distributed::Triangulation::register_data_attach() and
     * parallel::distributed::Triangulation::notify_ready_for_unpack() interface
     * with the aim of serialization, the calls to the corresponding
     * prepare_for_serialization() and deserialize() functions need to happen in
     * the same order, respectively. Classes relying on this interface are e.g.
     * parallel::distributed::CellDataTransfer,
     * parallel::distributed::SolutionTransfer, and Particles::ParticleHandler.
     *
     * @note See the documentation of parallel::distributed::SolutionTransfer for
     * matching code snippets for both transfer and serialization.
     *
     * @ingroup distributed
     */
    template <int dim, int spacedim = dim, typename VectorType = Vector<double>>
    class CellDataTransfer
    {
    private:
      /**
       * An alias that defines the data type of provided container template.
       */
      using value_type = typename VectorType::value_type;

    public:
      /**
       * Constructor.
       *
       * @param[in] triangulation The triangulation on which all operations will
       *   happen. At the time when this constructor is called, the refinement
       *   in question has not happened yet.
       * @param[in] transfer_variable_size_data Specify whether your VectorType
       *   container stores values that differ in size. A varying amount of data
       *   may be packed per cell, if for example the underlying ValueType of
       *   the VectorType container is a container itself.
       * @param[in] refinement_strategy %Function deciding how data will be
       *   stored on refined cells from its parent cell.
       * @param[in] coarsening_strategy %Function deciding which data to store
       *   on a cell whose children will get coarsened into.
       */
      CellDataTransfer(
        const parallel::distributed::Triangulation<dim, spacedim>
                                         &triangulation,
        const bool                        transfer_variable_size_data = false,
        const std::function<std::vector<value_type>(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                          &parent,
          const value_type parent_value)> refinement_strategy =
          &dealii::AdaptationStrategies::Refinement::
            preserve<dim, spacedim, value_type>,
        const std::function<value_type(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator
                                        &parent,
          const std::vector<value_type> &children_values)> coarsening_strategy =
          &dealii::AdaptationStrategies::Coarsening::
            check_equality<dim, spacedim, value_type>);

      /**
       * Prepare the current object for coarsening and refinement.
       *
       * It registers the data transfer of @p in on the underlying triangulation.
       * @p in includes data to be interpolated onto the new (refined and/or
       * coarsened) grid. See documentation of this class for more information
       * on how to use this functionality.
       *
       * This function can be called only once for the specified container
       * until data transfer has been completed. If multiple vectors shall be
       * transferred via this class, use the function below.
       */
      void
      prepare_for_coarsening_and_refinement(const VectorType &in);

      /**
       * Same as the function above, only for a list of vectors.
       */
      void
      prepare_for_coarsening_and_refinement(
        const std::vector<const VectorType *> &all_in);

      /**
       * Prepare the serialization of the given vector.
       *
       * The serialization is done by Triangulation::save(). See documentation
       * of this class for more information on how to use this functionality.
       *
       * This function can be called only once for the specified container
       * until data transfer has been completed. If multiple vectors shall be
       * transferred via this class, use the function below.
       */
      void
      prepare_for_serialization(const VectorType &in);

      /**
       * Same as the function above, only for a list of vectors.
       */
      void
      prepare_for_serialization(const std::vector<const VectorType *> &all_in);

      /**
       * Unpack the information previously stored in this object before
       * the mesh was refined or coarsened onto the current set of cells.
       */
      void
      unpack(VectorType &out);

      /**
       * Same as the function above, only for a list of vectors.
       */
      void
      unpack(std::vector<VectorType *> &all_out);

      /**
       * Execute the deserialization of the stored information.
       * This needs to be done after calling Triangulation::load().
       */
      void
      deserialize(VectorType &out);

      /**
       * Same as the function above, only for a list of vectors.
       */
      void
      deserialize(std::vector<VectorType *> &all_out);

    private:
      /**
       * Pointer to the triangulation to work with.
       */
      ObserverPointer<const parallel::distributed::Triangulation<dim, spacedim>,
                      CellDataTransfer<dim, spacedim, VectorType>>
        triangulation;

      /**
       * Specifies if size of data to transfer varies from cell to cell.
       */
      const bool transfer_variable_size_data;

      /**
       * %Function deciding how data will be stored on refined cells from its
       * parent cell.
       */
      const std::function<std::vector<value_type>(
        const typename Triangulation<dim, spacedim>::cell_iterator &parent,
        const value_type parent_value)>
        refinement_strategy;

      /**
       * %Function deciding on how to process data from children to be stored on
       * the parent cell.
       */
      const std::function<value_type(
        const typename Triangulation<dim, spacedim>::cell_iterator &parent,
        const std::vector<value_type> &children_values)>
        coarsening_strategy;

      /**
       * A vector that stores pointers to all the vectors we are supposed to
       * copy over from the old to the new mesh.
       */
      std::vector<const VectorType *> input_vectors;

      /**
       * The handle that triangulation has assigned to this object
       * with which we can access our memory offset and our pack function.
       */
      unsigned int handle;

      /**
       * Registers the pack_callback() function to the triangulation
       * and stores the returning handle.
       */
      void
      register_data_attach();

      /**
       * A callback function used to pack the data on the current mesh into
       * objects that can later be retrieved after refinement, coarsening and
       * repartitioning.
       */
      std::vector<char>
      pack_callback(const typename parallel::distributed::
                      Triangulation<dim, spacedim>::cell_iterator &cell,
                    const CellStatus                               status);

      /**
       * A callback function used to unpack the data on the current mesh that
       * has been packed up previously on the mesh before refinement,
       * coarsening and repartitioning.
       */
      void
      unpack_callback(
        const typename parallel::distributed::Triangulation<dim, spacedim>::
          cell_iterator &cell,
        const CellStatus status,
        const boost::iterator_range<std::vector<char>::const_iterator>
                                  &data_range,
        std::vector<VectorType *> &all_out);
    };
  } // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_distributed_cell_data_transfer_h */
