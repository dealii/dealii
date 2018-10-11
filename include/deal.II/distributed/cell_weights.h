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

#ifndef dealii_distributed_cell_weights_h
#define dealii_distributed_cell_weights_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/hp/dof_handler.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  /**
   * Anytime a parallel::Triangulation is repartitioned, either upon request
   * or by refinement/coarsening, cells will be distributed amongst all
   * subdomains to achieve an equally balanced workload. If the workload per
   * cell varies, which is in general the case for hp::DoFHandler objects, we
   * can take that into account by introducing individual weights for
   * different cells.
   *
   * This class allows computing these weights for load balancing by
   * consulting the FiniteElement that is associated with each cell of
   * a hp::DoFHandler. One can choose from predefined weighting
   * algorithms provided by this class or provide a custom one. The
   * chosen weighting function will be connected to the corresponding
   * signal of the linked parallel::Triangulation via callback.
   *
   * An object of this class needs to exist for every DoFHandler associated
   * with the Triangulation we work on to achieve satisfying work balancing
   * results.
   *
   * A small code snippet follows explaining how to achieve each cell
   * being weighted by its current number of degrees of freedom.
   * @code
   * parallel::CellWeights<dim, spacedim> cell_weights(hp_dof_handler);
   * cell_weights.register_ndofs_weighting();
   * @endcode
   * The weighting function can be changed anytime. Even more ambitious
   * approaches are possible by submitting customized functions, e.g.
   * @code
   * cell_weights.register_custom_weighting(
   *   [](const FiniteElement<dim, spacedim> &active_fe,
   *      const typename hp::DoFHandler<dim, spacedim>::cell_iterator &)
   *   -> unsigned int {
   *   return 1000 * std::pow(active_fe.dofs_per_cell, 1.6);
   * });
   * @endcode
   * The returned value has to be an unsigned integer and is thus limited by
   * some large number. It is interpreted as the additional computational load
   * of each cell. See Triangulation::Signals::cell_weight for a discussion on
   * this topic.
   *
   * @note Be aware that this class connects the weight function to the
   * Triangulation during its construction. If the Triangulation
   * associated with the DoFHandler changes during the lifetime of the
   * latter, an assertion will be triggered in the weight_callback() function.
   * It is recommended to create a separate object in this case and to destroy
   * the previous one.
   *
   * @ingroup distributed
   * @author Marc Fehling, 2018
   */
  template <int dim, int spacedim = dim>
  class CellWeights
  {
  public:
    /**
     * Constructor.
     *
     * @param[in] dof_handler The hp::DoFHandler which will be used to
     *    determine each cell's finite element.
     */
    CellWeights(const dealii::hp::DoFHandler<dim, spacedim> &dof_handler);

    /**
     * Destructor.
     */
    ~CellWeights();

    /**
     * Choose a constant weight @p factor on each cell.
     */
    void
    register_constant_weighting(const unsigned int factor = 1000);

    /**
     * Choose a weight for each cell that is proportional to its number of
     * degrees of freedom with a factor @p factor.
     */
    void
    register_ndofs_weighting(const unsigned int factor = 1000);

    /**
     * Choose a weight for each cell that is proportional to its number of
     * degrees of freedom <i>squared</i> with a factor @p factor.
     */
    void
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
    void
    register_custom_weighting(
      const std::function<unsigned int(
        const FiniteElement<dim, spacedim> &,
        const typename hp::DoFHandler<dim, spacedim>::cell_iterator &)>
        custom_function);

  private:
    /**
     * Pointer to the degree of freedom handler.
     */
    SmartPointer<const dealii::hp::DoFHandler<dim, spacedim>, CellWeights>
      dof_handler;

    /**
     * Pointer to the Triangulation associated with the degree of freedom
     * handler.
     *
     * We store both to make sure to always work on the correct combination of
     * both.
     */
    SmartPointer<const parallel::Triangulation<dim, spacedim>, CellWeights>
      triangulation;

    /**
     * Function that will determine each cell's weight.
     *
     * Can be set using the register_constant_weighting(),
     * register_ndofs_weighting(), register_ndofs_squared_weighting(), and
     * register_custom_weighting() member functions.
     *
     * The function requires the active FiniteElement object on each cell
     * as an argument, as well as the cell itself of type
     * hp::DoFHandler::cell_iterator.
     */
    std::function<unsigned int(
      const FiniteElement<dim, spacedim> &,
      const typename hp::DoFHandler<dim, spacedim>::cell_iterator &)>
      weighting_function;

    /**
     * A connection to the Triangulation of the DoFHandler.
     */
    boost::signals2::connection tria_listener;

    /**
     * A callback function that will be attached to the cell_weight signal of
     * the Triangulation, that is a member of the DoFHandler. Ultimately
     * returns the weight for each cell, determined by the weighting_function
     * member.
     */
    unsigned int
    weight_callback(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename Triangulation<dim, spacedim>::CellStatus     status);
  };
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
