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

#ifndef dealii_distributed_cell_weights_h
#define dealii_distributed_cell_weights_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <boost/signals2/connection.hpp>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  /**
   * Anytime a parallel::TriangulationBase is repartitioned, either upon request
   * or by refinement/coarsening, cells will be distributed amongst all
   * subdomains to achieve an equally balanced workload. By default, if you
   * don't do anything specific, the workload is assumed to be the same for
   * every cell of the triangulation, and in that case a well-balanced mesh
   * will have about the same number of cells on each process. But, if the
   * workload per cell varies, which is in general the case for DoFHandler
   * objects that hp-capabilities (i.e., where you use different finite elements
   * on different cells), we can take that into account by introducing
   * individual weights for different cells.
   *
   * This class allows computing these weights for load balancing by
   * consulting the FiniteElement that is associated with each cell of
   * a DoFHandler. One can choose from predefined weighting algorithms provided
   * by this class or provide a custom one.
   *
   * If the associated DoFHandler has not been initialized yet, i.e., its
   * hp::FECollection is empty, all cell weights will be evaluated as zero.
   *
   * This class offers two different ways of connecting the chosen weighting
   * function to the corresponding signal of the linked
   * parallel::TriangulationBase. The recommended way involves creating an
   * object of this class which will automatically take care of registering the
   * weighting function upon creation and de-registering it once destroyed. An
   * object of this class needs to exist for every DoFHandler associated with
   * the Triangulation we work on to achieve satisfying work balancing results.
   * The connected weighting function may be changed anytime using the
   * CellWeights::reinit() function. The following code snippet demonstrates how
   * to achieve each cell being weighted by its current number of degrees of
   * freedom.
   * @code
   * parallel::CellWeights<dim, spacedim> cell_weights(
   *   hp_dof_handler,
   *   parallel::CellWeights<dim, spacedim>::ndofs_weighting({1, 1}));
   * @endcode
   *
   * On the other hand, you are also able to take care of handling the signal
   * connection manually by using the static member function of this class. In
   * this case, an analogous code example looks as follows.
   * @code
   * boost::signals2::connection connection =
   *   hp_dof_handler.get_triangulation().signals.weight.connect(
   *     parallel::CellWeights<dim, spacedim>::make_weighting_callback(
   *       hp_dof_handler,
   *       parallel::CellWeights<dim, spacedim>::ndofs_weighting({1, 1}));
   * @endcode
   *
   * The use of this class is demonstrated in step-75.
   *
   * @note See Triangulation::Signals::weight for more information on
   * weighting and load balancing.
   *
   * @note Be aware that this class connects the weight function to the
   * Triangulation during this class's constructor. If the Triangulation
   * associated with the DoFHandler changes during the lifetime of the
   * latter via DoFHandler::reinit(), an assertion will be triggered in
   * the weight_callback() function. Use CellWeights::reinit() to deregister the
   * weighting function on the old Triangulation and connect it to the new one.
   *
   * @note A hp::FECollection needs to be attached to your DoFHandler object via
   * DoFHandler::distribute_dofs() <em>before</em> the
   * Triangulation::Signals::weight signal will be triggered. Otherwise,
   * your DoFHandler does not know many degrees of freedom your cells have. In
   * other words, you need to call DoFHandler::distribute_dofs() once before you
   * call
   * parallel::distributed::Triangulation::execute_coarsening_and_refinement(),
   * parallel::distributed::Triangulation::refine_global(),
   * parallel::distributed::Triangulation::repartition(), or
   * GridTools::partition_triangulation() for the very first time.
   *
   * @ingroup distributed
   */
  template <int dim, int spacedim = dim>
  class CellWeights
  {
  public:
    /**
     * An alias that defines the characteristics of a function that can be used
     * for weighting cells during load balancing.
     *
     * Such weighting functions take as arguments an iterator to a cell and the
     * future finite element that will be assigned to it after repartitioning.
     * They return an unsigned integer which is interpreted as the cell's
     * weight or, in other words, the additional computational load associated
     * with it.
     */
    using WeightingFunction = std::function<
      unsigned int(const typename DoFHandler<dim, spacedim>::cell_iterator &,
                   const FiniteElement<dim, spacedim> &)>;

    /**
     * Constructor.
     *
     * No weighting function will be connected yet. Please call reinit().
     */
    CellWeights() = default;

    /**
     * Constructor.
     *
     * @param[in] dof_handler The DoFHandler which will be used to
     *    determine each cell's finite element.
     * @param[in] weighting_function The function that determines each
     *    cell's weight during load balancing.
     */
    CellWeights(const DoFHandler<dim, spacedim> &dof_handler,
                const WeightingFunction         &weighting_function);

    /**
     * Destructor.
     *
     * Disconnects the function previously connected to the weighting signal.
     */
    ~CellWeights();

    /**
     * Connect a different @p weighting_function to the Triangulation
     * associated with the @p dof_handler.
     *
     * Disconnects the function previously connected to the weighting signal.
     */
    void
    reinit(const DoFHandler<dim, spacedim> &dof_handler,
           const WeightingFunction         &weighting_function);

    /**
     * Converts a @p weighting_function to a different type that qualifies as
     * a callback function, which can be connected to a weighting signal of a
     * Triangulation.
     *
     * This function does <b>not</b> connect the converted function to the
     * Triangulation associated with the @p dof_handler.
     */
    static std::function<unsigned int(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
      const CellStatus status)>
    make_weighting_callback(const DoFHandler<dim, spacedim> &dof_handler,
                            const WeightingFunction &weighting_function);

    /**
     * @name Selection of weighting functions
     * @{
     */

    /**
     * Choose a constant weight @p factor on each cell.
     */
    static WeightingFunction
    constant_weighting(const unsigned int factor = 1);

    /**
     * The pair of floating point numbers $(a,b)$ provided via
     * @p coefficients determines the weight $w_K$ of each cell $K$ with
     * $n_K$ degrees of freedom in the following way: \f[ w_K =
     * a \, n_K^b \f]
     *
     * The right hand side will be rounded to the nearest integer since cell
     * weights are required to be integers.
     */
    static WeightingFunction
    ndofs_weighting(const std::pair<float, float> &coefficients);

    /**
     * The container @p coefficients provides pairs of floating point numbers
     * $(a_i, b_i)$ that determine the weight $w_K$ of each cell
     * $K$ with $n_K$ degrees of freedom in the following way: \f[ w_K =
     * \sum_i a_i \, n_K^{b_i} \f]
     *
     * The right hand side will be rounded to the nearest integer since cell
     * weights are required to be integers.
     */
    static WeightingFunction
    ndofs_weighting(const std::vector<std::pair<float, float>> &coefficients);

    /**
     * @}
     */

  private:
    /**
     * A connection to the corresponding `weight` signal of the Triangulation
     * which is attached to the DoFHandler.
     */
    boost::signals2::connection connection;

    /**
     * A callback function that will be connected to the `weight` signal of
     * the @p triangulation, to which the @p dof_handler is attached. Ultimately
     * returns the weight for each cell, determined by the @p weighting_function
     * provided as a parameter. Returns zero if @p dof_handler has not been
     * initialized yet.
     */
    static unsigned int
    weighting_callback(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
      const CellStatus                                  status,
      const DoFHandler<dim, spacedim>                  &dof_handler,
      const parallel::TriangulationBase<dim, spacedim> &triangulation,
      const WeightingFunction                          &weighting_function);
  };
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
