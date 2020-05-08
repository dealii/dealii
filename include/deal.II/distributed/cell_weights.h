// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_distributed_cell_weights_h
#define dealii_distributed_cell_weights_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/hp/dof_handler.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  /**
   * Anytime a parallel::TriangulationBase is repartitioned, either upon request
   * or by refinement/coarsening, cells will be distributed amongst all
   * subdomains to achieve an equally balanced workload. If the workload per
   * cell varies, which is in general the case for hp::DoFHandler objects, we
   * can take that into account by introducing individual weights for
   * different cells.
   *
   * This class allows computing these weights for load balancing by
   * consulting the FiniteElement that is associated with each cell of
   * a hp::DoFHandler. One can choose from predefined weighting
   * algorithms provided by this class or provide a custom one.
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
   * freedom. We chose a factor of `1000` that corresponds to the initial weight
   * each cell is assigned to upon creation.
   * @code
   * parallel::CellWeights<dim, spacedim> cell_weights(
   *   hp_dof_handler,
   *   parallel::CellWeights<dim, spacedim>::ndofs_weighting({1000, 1}));
   * @endcode
   *
   * On the other hand, you are also able to take care of handling the signal
   * connection manually by using the static member function of this class. In
   * this case, an analogous code example looks as follows.
   * @code
   * boost::signals2::connection connection =
   *   hp_dof_handler.get_triangulation().signals.cell_weight.connect(
   *     parallel::CellWeights<dim, spacedim>::make_weighting_callback(
   *       hp_dof_handler,
   *       parallel::CellWeights<dim, spacedim>::ndofs_weighting(
   *         {1000, 1}));
   * @endcode
   *
   * @note See Triangulation::Signals::cell_weight for more information on
   * weighting and load balancing.
   *
   * @note Be aware that this class connects the weight function to the
   * Triangulation during this class's constructor. If the Triangulation
   * associated with the DoFHandler changes during the lifetime of the
   * latter via hp::DoFHandler::initialize(), an assertion will be triggered in
   * the weight_callback() function. Use CellWeights::reinit() to deregister the
   * weighting function on the old Triangulation and connect it to the new one.
   *
   * @ingroup distributed
   * @author Marc Fehling, 2018, 2019
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
    using WeightingFunction = std::function<unsigned int(
      const typename hp::DoFHandler<dim, spacedim>::cell_iterator &,
      const FiniteElement<dim, spacedim> &)>;

    /**
     * Constructor.
     *
     * @param[in] dof_handler The hp::DoFHandler which will be used to
     *    determine each cell's finite element.
     * @param[in] weighting_function The function that determines each
     *    cell's weight during load balancing.
     */
    CellWeights(const hp::DoFHandler<dim, spacedim> &dof_handler,
                const WeightingFunction &            weighting_function);

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
    reinit(const hp::DoFHandler<dim, spacedim> &dof_handler,
           const WeightingFunction &            weighting_function);

    /**
     * Converts a @p weighting_function to a different type that qualifies as
     * a callback function, which can be connected to a weighting signal of a
     * Triangulation.
     *
     * This function does <b>not</b> connect the converted function to the
     * Triangulation associated with the @p dof_handler.
     */
    static std::function<unsigned int(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename Triangulation<dim, spacedim>::CellStatus     status)>
    make_weighting_callback(const hp::DoFHandler<dim, spacedim> &dof_handler,
                            const WeightingFunction &weighting_function);

    /**
     * @name Selection of weighting functions
     * @{
     */

    /**
     * Choose a constant weight @p factor on each cell.
     */
    static WeightingFunction
    constant_weighting(const unsigned int factor = 1000);

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

    /**
     * @name Deprecated functions
     * @{
     */

    /**
     * Constructor.
     *
     * @param[in] dof_handler The hp::DoFHandler which will be used to
     *    determine each cell's finite element.
     */
    DEAL_II_DEPRECATED
    CellWeights(const hp::DoFHandler<dim, spacedim> &dof_handler);

    /**
     * Choose a constant weight @p factor on each cell.
     *
     * @deprecated Use CellWeights::constant_weighting() instead.
     */
    DEAL_II_DEPRECATED void
    register_constant_weighting(const unsigned int factor = 1000);

    /**
     * Choose a weight for each cell that is proportional to its number of
     * degrees of freedom with a factor @p factor.
     *
     * @deprecated Use CellWeights::ndofs_weighting() instead.
     */
    DEAL_II_DEPRECATED void
    register_ndofs_weighting(const unsigned int factor = 1000);

    /**
     * Choose a weight for each cell that is proportional to its number of
     * degrees of freedom <i>squared</i> with a factor @p factor.
     *
     * @deprecated Use CellWeights::ndofs_weighting() instead.
     */
    DEAL_II_DEPRECATED void
    register_ndofs_squared_weighting(const unsigned int factor = 1000);

    /**
     * Register a custom weight for each cell by providing a function as a
     * parameter.
     *
     * @param[in] custom_function The custom weighting function returning
     *    the weight of each cell as an unsigned integer. It is required
     *    to have two arguments, namely the FiniteElement that will be
     *    active on the particular cell, and the cell itself of type
     *    hp::DoFHandler::cell_iterator. We require both to make sure to
     *    get the right active FiniteElement on each cell in case that we
     *    coarsen the Triangulation.
     */
    DEAL_II_DEPRECATED void
    register_custom_weighting(
      const std::function<unsigned int(
        const FiniteElement<dim, spacedim> &,
        const typename hp::DoFHandler<dim, spacedim>::cell_iterator &)>
        custom_function);

    /**
     * @}
     */

  private:
    /**
     * @name Deprecated members
     * @{
     */

    /**
     * Pointer to the degree of freedom handler.
     */
    DEAL_II_DEPRECATED
    SmartPointer<const dealii::hp::DoFHandler<dim, spacedim>, CellWeights>
      dof_handler;

    /**
     * Pointer to the Triangulation associated with the degree of freedom
     * handler.
     *
     * We store both to make sure to always work on the correct combination of
     * both.
     */
    DEAL_II_DEPRECATED
    SmartPointer<const parallel::TriangulationBase<dim, spacedim>, CellWeights>
      triangulation;

    /**
     * @}
     */

    /**
     * A connection to the corresponding cell_weight signal of the Triangulation
     * which is attached to the DoFHandler.
     */
    boost::signals2::connection connection;

    /**
     * A callback function that will be connected to the cell_weight signal of
     * the @p triangulation, to which the @p dof_handler is attached. Ultimately
     * returns the weight for each cell, determined by the @p weighting_function
     * provided as a parameter.
     */
    static unsigned int
    weighting_callback(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename Triangulation<dim, spacedim>::CellStatus     status,
      const hp::DoFHandler<dim, spacedim> &                       dof_handler,
      const parallel::TriangulationBase<dim, spacedim> &          triangulation,
      const WeightingFunction &weighting_function);
  };
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
