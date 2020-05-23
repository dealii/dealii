// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#ifndef dealii_cell_data_transfer_h
#define dealii_cell_data_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/grid/tria.h>

#include <deal.II/numerics/adaptation_strategies.h>

#include <algorithm>
#include <functional>
#include <set>


DEAL_II_NAMESPACE_OPEN

/**
 * Transfer data that is associated with each active cell (like error
 * indicators) while refining and/or coarsening a triangulation.
 *
 * This class therefore does for cell-related information what
 * SolutionTransfer does for the values of degrees of freedom defined on a
 * Triangulation.
 *
 * A non-distributed container (like Vector or `std::vector`) has to be
 * provided, which holds the cell-wise data in the same order as active cells
 * are traversed. In other words, each entry corresponds to the cell with the
 * same index CellAccessor::active_cell_index(), and the container has to be of
 * size Triangulation::n_active_cells().
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
 * CellDataTransfer<dim, spacedim, Vector<double>>
 *   cell_data_trans(triangulation);
 * cell_data_trans.prepare_for_coarsening_and_refinement();
 *
 * // actually execute the refinement,
 * triangulation.execute_coarsening_and_refinement();
 *
 * // unpack transferred data,
 * Vector<double> transferred_data(triangulation.n_active_cells());
 * cell_data_trans.unpack(data_to_transfer, transferred_data);
 * @endcode
 *
 * When using a parallel::shared::Triangulation, we need to ensure that we have
 * the global data available in our local vector before refinement happened. We
 * can achieve this as follows:
 *
 * @code
 * Vector<double> data_to_transfer(triangulation.n_active_cells());
 * //[fill data_to_transfer with cell-wise values...]
 *
 * PETScWrappers::MPI::Vector
 * distributed_data_to_transfer(mpi_communicator,
 *                              triangulation.n_active_cells(),
 *                              triangulation.n_locally_owned_active_cells());
 * for (const auto &cell : triangulation.active_cell_iterators())
 *   if (cell->is_locally_owned())
 *     {
 *       const unsigned int index = cell->active_cell_index();
 *       distributed_data_to_transfer(index) = data_to_transfer(index);
 *     }
 * distributed_data_to_transfer.compress(VectorOperation::insert);
 *
 * data_to_transfer = distributed_data_to_transfer;
 * @endcode
 *
 * For the parallel distributed case, a designated class
 * parallel::distributed::CellDataTransfer is available. Please refer to this
 * particular class when using a parallel::distributed::Triangulation.
 *
 * @note See the documentation of SolutionTransfer for matching code snippets
 *   for transfer.
 *
 * @ingroup numerics
 * @author Marc Fehling, 2019 - 2020
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
   * @param[in] refinement_strategy Function deciding how data will be stored on
   *   refined cells from its parent cell.
   * @param[in] coarsening_strategy Function deciding which data to store on
   *   a cell whose children will get coarsened into.
   */
  CellDataTransfer(
    const Triangulation<dim, spacedim> &triangulation,
    const std::function<std::vector<value_type>(
      const typename Triangulation<dim, spacedim>::cell_iterator &parent,
      const value_type parent_value)>   refinement_strategy =
      &AdaptationStrategies::Refinement::preserve<dim, spacedim, value_type>,
    const std::function<value_type(
      const typename Triangulation<dim, spacedim>::cell_iterator &parent,
      const std::vector<value_type> &children_values)> coarsening_strategy =
      &AdaptationStrategies::Coarsening::
        check_equality<dim, spacedim, value_type>);

  /**
   * Prepare the current object for coarsening and refinement.
   *
   * Stores the active_cell_indices of all active cells on the associated
   * triangulation and attribute them to either persisting, refined or coarsened
   * cells.
   */
  void
  prepare_for_coarsening_and_refinement();

  /**
   * Transfer the information from the previous mesh to the updated one.
   *
   * Data from the previous mesh supplied by @p in will be transferred to the updated
   * mesh and stored in @p out. @p out has to provide enough space to hold the
   * transferred data, i.e. has to be of size `triangulation.n_active_cells()`.
   */
  void
  unpack(const VectorType &in, VectorType &out);

private:
  /**
   * Pointer to the triangulation to work with.
   */
  SmartPointer<const Triangulation<dim, spacedim>,
               CellDataTransfer<dim, spacedim, VectorType>>
    triangulation;

  /**
   * Function deciding how data will be stored on refined cells from its parent
   * cell.
   */
  const std::function<std::vector<value_type>(
    const typename Triangulation<dim, spacedim>::cell_iterator &parent,
    const value_type                                            parent_value)>
    refinement_strategy;

  /**
   * Function deciding on how to process data from children to be stored on the
   * parent cell.
   */
  const std::function<value_type(
    const typename Triangulation<dim, spacedim>::cell_iterator &parent,
    const std::vector<value_type> &children_indices)>
    coarsening_strategy;

  /**
   * Container to temporarily store the iterator and active cell index
   * of cells that persist.
   */
  std::map<const typename Triangulation<dim, spacedim>::cell_iterator,
           const unsigned int>
    persisting_cells_active_index;

  /**
   * Container to temporarily store the iterator and active cell index
   * of cells that will be refined.
   */
  std::map<const typename Triangulation<dim, spacedim>::cell_iterator,
           const unsigned int>
    refined_cells_active_index;

  /**
   * Container to temporarily store the iterator of parent cells that will
   * remain after coarsening along with the active cell indices of the
   * corresponding children cells.
   */
  std::map<const typename Triangulation<dim, spacedim>::cell_iterator,
           const std::set<unsigned int>>
    coarsened_cells_active_index;

  /**
   * Number of active cells on the initial triangulation that has not been
   * refined yet.
   *
   * It will be set in prepare_for_coarsening_and_refinement() and used to
   * validate user inputs after refinement happened (only in debug mode).
   */
  unsigned int n_active_cells_pre;
};


DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_cell_data_transfer_h */
