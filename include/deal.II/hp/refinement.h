// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_hp_refinement_h
#define dealii_hp_refinement_h


#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>

#include <functional>


DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename Number>
class Vector;

namespace hp
{
  template <int dim, int spacedim>
  class DoFHandler;
}
#endif


namespace hp
{
  /**
   * We supply adaptive methods to align computational ressources with the
   * complexity of the numerical solution. Error estimates are an appropriate
   * means of determining where adjustments need to be made.
   *
   * However with hp-adaptivity, we have two ways to realize these adjustments:
   * For irregular solutions, h-adaptive methods which dynamically assign cell
   * sizes tend to reduce the approximation error, while for smooth solutions
   * p-adaptive methods are better suited in which function spaces will be
   * selected dynamically. This namespace collects tools to decide which type
   * of adaptive methods to apply.
   *
   * <h3>Usage</h3>
   *
   * To successfully apply hp-adaptive methods, we recommend the following
   * workflow:
   * <ol>
   * <li> A suitable error estimate is the basis for any kind of adaptive method.
   * Similar to pure grid refinement, we will determine error estimates in the
   * usual way (i.e. KellyErrorEstimator) and mark cells for refinement or
   * coarsening (i.e. GridRefinement).
   *
   * Calling Triangulation::execute_coarsening_and_refinement() at this stage
   * will perform pure grid refinement as expected.
   *
   * <li> Once all refinement and coarsening flags have been distributed on the
   * mesh, we may determine if those qualify for p-adaptive methods.
   * Corresponding functions will set @p future_fe_indices on top of the
   * refinement and coarsening flags if they fulfil a certain criterion.
   *
   * In case of refinement, the superordinate element of the underlying
   * hp::FECollection will be assigned as the future finite element.
   * Correspondingly, the subordinate element will be selected for coarsening.
   *
   * Triangulation::execute_coarsening_and_refinement() will now supply both
   * h- and p-adaptive methods independently.
   *
   * <li> Right now, there may be cells scheduled for both h- and p-adaptation.
   * If we do not want to impose both methods at once, we need to decide which
   * one to pick for each cell individually and unambiguously. Since grid
   * refinement will be imposed by default and we only determine qualification
   * for p-adaptivity on top, we will always decide in favour of p-adaptive
   * methods.
   *
   * Calling Triangulation::execute_coarsening_and_refinement() will now perform
   * either h- or p-adaptive methods uniquely on each cell.
   *
   * <li> Up to this point, each cell knows its destiny in terms of adaptivity.
   * We can now move on to prepare all data structures to be transferred across
   * mesh changes. Previously set refinement and coarsening flags as well as
   * @p future_fe_indices will be used to update the data accordingly.
   * </ol>
   *
   * As an example, a realisation of pure p-adaptive methods would look like the
   * following:
   * @code
   * // step 1: flag cells for refinement or coarsening
   * Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
   * KellyErrorEstimator<dim>::estimate(hp_dof_handler,
   *                                    QGauss<dim-1> (quadrature_points),
   *                                    typename FunctionMap<dim>::type(),
   *                                    solution,
   *                                    estimated_error_per_cell);
   * GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
   *                                                   estimated_error_per_cell,
   *                                                   top_fraction,
   *                                                   bottom_fraction);
   *
   * // step 2: set future finite element indices on flagged cells
   * hp::Refinement::full_p_adaptivity(hp_dof_handler);
   *
   * // step 3: decide whether h- or p-adaptive methods will be supplied
   * hp::Refinement::force_p_over_h(hp_dof_handler);
   *
   * // step 4: prepare solutions to be transferred
   * ...
   *
   * triangulation.execute_coarsening_and_refinement();
   * @endcode
   *
   * @ingroup hp
   * @author Marc Fehling 2019
   */
  namespace Refinement
  {
    /**
     * An alias that defines the characteristics of a function that can be used
     * as a comparison criterion for deciding whether to perform h- or
     * p-adaptation.
     *
     * Such functions take two numbers as arguments: The first one corresponds
     * to the provided criterion, while the other one conforms to the reference.
     * The result of the comparision will be returned as a boolean.
     */
    template <typename Number>
    using ComparisonFunction =
      std::function<bool(const Number &, const Number &)>;

    /**
     * @name Setting p-adaptivity flags
     * @{
     */

    /**
     * Each cell flagged for h-refinement will also be flagged for p-refinement.
     * The same applies to coarsening.
     *
     * @note Preceeding calls of Triangulation::prepare_for_coarsening_and_refinement()
     *   may change refine and coarsen flags, which will ultimately change the
     *   results of this function.
     */
    template <int dim, int spacedim>
    void
    full_p_adaptivity(const hp::DoFHandler<dim, spacedim> &dof_handler);

    /**
     * Adapt which finite element to use on cells that have been specifically
     * flagged for p-adaptation via the parameter @p p_flags. Future finite
     * elements will only be assigned if cells have been flagged for refinement
     * and coarsening beforehand.
     *
     * Each entry of the parameter @p p_flags needs to correspond to an active
     * cell.
     *
     * @note Preceeding calls of Triangulation::prepare_for_coarsening_and_refinement()
     *   may change refine and coarsen flags, which will ultimately change the
     *   results of this function.
     */
    template <int dim, int spacedim>
    void
    p_adaptivity_from_flags(const hp::DoFHandler<dim, spacedim> &dof_handler,
                            const std::vector<bool> &            p_flags);

    /**
     * Adapt which finite element to use on cells whose criteria meet a certain
     * absolute threshold.
     *
     * For p-refinement and p-coarsening, two separate thresholds need to
     * provided via parameters @p p_refine_threshold and @p p_coarsen_threshold.
     *
     * We consider a cell for p-adaptivity if it is currently flagged for
     * refinement or coarsening and its criterion successfully compares to the
     * corresponding threshold. Let us be more specific on the default case: We
     * consider a cell for p-refinement if it is flagged for refinement and its
     * criterion is larger than the corresponding threshold. The same applies
     * for p-coarsening, but the cell's criterion must be lower than the
     * threshold. However, different compare function objects can be supplied
     * via the parameters @p compare_refine and @p compare_coarsen to impose
     * different decision strategies.
     *
     * Each entry of the parameter @p criteria needs to correspond to an active
     * cell.
     *
     * @note Preceeding calls of Triangulation::prepare_for_coarsening_and_refinement()
     *   may change refine and coarsen flags, which will ultimately change the
     *   results of this function.
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_absolute_threshold(
      const hp::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &               criteria,
      const Number                         p_refine_threshold,
      const Number                         p_coarsen_threshold,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_refine = std::greater<Number>(),
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen = std::less<Number>());

    /**
     * Adapt which finite element to use on cells whose criteria meet a certain
     * threshold relative to the overall range of criterion values.
     *
     * The threshold will be determined for refined and coarsened cells
     * separately based on the currently set refinement markers. For each class
     * of cells, we determine the maximal and minimal values of all criteria and
     * determine the threshold by linear interpolation between these limits.
     * Parameters @p p_refine_fraction and @p p_refine_coarsen are used as
     * interpolation factors, where `0` corresponds to the minimal and `1` to
     * the maximal value. By default, mean values are considered as thresholds.
     *
     * We consider a cell for p-adaptivity if it is currently flagged for
     * refinement or coarsening and its criterion successfully compares to the
     * corresponding threshold. Let us be more specific on the default case: We
     * consider a cell for p-refinement if it is flagged for refinement and its
     * criterion is larger than the corresponding threshold. The same applies
     * for p-coarsening, but the cell's criterion must be lower than the
     * threshold. However, different compare function objects can be supplied
     * via the parameters @p compare_refine and @p compare_coarsen to impose
     * different decision strategies.
     *
     * Each entry of the parameter @p criteria needs to correspond to an active
     * cell. Parameters @p p_refine_fraction and @p p_coarsen_fraction need to be
     * in the interval $[0,1]$.
     *
     * @note Preceeding calls of Triangulation::prepare_for_coarsening_and_refinement()
     *   may change refine and coarsen flags, which will ultimately change the
     *   results of this function.
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_relative_threshold(
      const hp::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &               criteria,
      const double                         p_refine_fraction  = 0.5,
      const double                         p_coarsen_fraction = 0.5,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_refine = std::greater<Number>(),
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen = std::less<Number>());

    /**
     * Adapt which finite element to use on cells based on the regularity of the
     * (unknown) analytical solution.
     *
     * With an approximation of the local Sobolev regularity index $k_K$,
     * we may assess to which finite element space our local solution on cell
     * $K$ belongs. Since the regularity index is only an estimate, we won't
     * use it to assign the finite element space directly, but rather consider
     * it as an indicator for adaptation. If a cell is flagged for refinement,
     * we will perform p-refinement once it satisfies
     * $k_K > p_{K,\text{super}}$, where $p_{K,\text{super}}$ is
     * the polynomial degree of the finite element superordinate to the
     * currently active element on cell $K$. In case of coarsening, the
     * criterion $k_K < p_{K,\text{sub}}$ has to be met, with
     * $p_{K,\text{sub}}$ the degree of the subordinate element.
     *
     * Each entry of the parameter @p sobolev_indices needs to correspond
     * to an active cell.
     *
     * For more theoretical details see
     * @code{.bib}
     * @article{Ainsworth1998,
     *  author    = {Ainsworth, Mark and Senior, Bill},
     *  title     = {An adaptive refinement strategy for hp-finite element
     *               computations},
     *  journal   = {{Applied Numerical Mathematics}},
     *  volume    = {26},
     *  number    = {1--2},
     *  pages     = {165--178},
     *  publisher = {Elsevier},
     *  year      = {1998},
     *  doi       = {10.1016/S0168-9274(97)00083-4}
     * }
     * @endcode
     *
     * @note Preceeding calls of Triangulation::prepare_for_coarsening_and_refinement()
     *   may change refine and coarsen flags, which will ultimately change the
     *   results of this function.
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_regularity(
      const hp::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &               sobolev_indices);

    /**
     * Adapt which finite element to use on each cell based on how its criterion
     * relates to a reference.
     *
     * We consider a cell for p-adaptivity if it is currently flagged for
     * refinement or coarsening and its criterion successfully compares to the
     * corresponding reference. Other than functions
     * p_adaptivity_from_absolute_threshold() and
     * p_adaptivity_from_relative_threshold(), compare function objects have to
     * be provided explicitly via the parameters @p compare_refine and
     * @p compare_coarsen.
     *
     * Each entry of the parameters @p criteria and @p references needs to
     * correspond to an active cell.
     *
     * @note Preceeding calls of Triangulation::prepare_for_coarsening_and_refinement()
     *   may change refine and coarsen flags, which will ultimately change the
     *   results of this function.
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_reference(
      const hp::DoFHandler<dim, spacedim> &                      dof_handler,
      const Vector<Number> &                                     criteria,
      const Vector<Number> &                                     references,
      const ComparisonFunction<typename identity<Number>::type> &compare_refine,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen);
    /**
     * @}
     */

    /**
     * @name Error prediction
     * @{
     */

    /**
     * Predict how the current @p error_indicators will adapt after refinement
     * and coarsening has happened on the provided @p dof_handler, and write its
     * results to @p predicted_errors. Each entry of @p error_indicators and
     * @p predicted_errors corresponds to an active cell on the underlying
     * Triangulation, thus each container has to be of size
     * Triangulation::n_active_cells(). The errors are interpreted to be
     * measured in the energy norm; this assumption enters the rate of
     * convergence that is used in the prediction. The `predicted_errors` output
     * argument has one entry per <i>current</i> cell, with the $2^d$ values for
     * each cell that will be coarsened away equal, and with the value stored on
     * a cell to be refined interpreted as applying to each of the future
     * children.
     *
     * For h-adaptation, we expect the local error $\eta_K$ on cell $K$ to be
     * proportional to $(h_K)^{p_K}$ in the energy norm, where $h_K$ denotes the
     * cell diameter and $p_K$ the polynomial degree of the currently assigned
     * finite element on cell $K$. Here, we assume that the finite element will
     * not change in the adaptation process so that $p_K = \text{const}$.
     * However during coarsening, the finite elements on siblings may be
     * different, and their parent cell will be assigned to their least
     * dominating finite element that belongs to its most general child. Thus,
     * we will always interpolate on an enclosing finite element space.
     * Additionaly assuming that the finite elements on the cells to be
     * coarsened are sufficient to represent the solution correctly (e.g. at
     * least quadratic basis functions for a quadratic solution), we are
     * confident to say that the error will not change by sole interpolation on
     * the larger finite element space.
     *
     * Further, the function assumes that the local error on a cell that will be
     * refined, will lead to errors on the $2^{dim}$ children that are all
     * equal, whereas local errors on siblings will be summed up on the parent
     * cell in case of coarsening. This assumption is often not satisfied in
     * practice: For example, if a cell is at a corner singularity, then the one
     * child cell that ends up closest to the singularity will inherit the
     * majority of the remaining error -- but this function can not know where
     * the singularity will be, and consequently assumes equal distribution.
     *
     * When transferring the predicted error to the coarsened mesh, make sure to
     * configure your CellDataTransfer object with CoarseningStrategies::sum()
     * as a coarsening strategy.
     *
     * For p-adaptation, the local error is expected to converge exponentially
     * with the polynomial degree of the assigned finite element. Each increase
     * or decrease of the degree will thus change its value by a user-defined
     * control parameter @p gamma_p. The assumption of exponential convergence
     * is only valid if both h- and p-adaptive methods are combined. An
     * exception is thrown if a cell is flagged for both h- and p-adaptation at
     * once.
     *
     * The prediction algorithm is formulated as follows with control parameters
     * @p gamma_p, @p gamma_h and @p gamma_n that may be used to influence
     * prediction for each adaptation type individually.
     * <table>
     *   <tr><th>Adaptation type <th colspan="2">Prediction formula
     *   <tr><td>no adaptation
     *       <td>$\eta_{K,\text{pred}} = \eta_{K} \, \gamma_\text{n}$
     *       <td>$\gamma_\text{n} \in (0,\infty)$
     *   <tr><td>p-adaptation
     *       <td>$\eta_{K,\text{pred}} = \eta_{K} \,
     *            \gamma_\text{p}^{(p_{K,\text{future}} - p_K)}$
     *       <td>$\gamma_\text{p} \in (0,1)$
     *   <tr><td>h-refinement
     *       <td>$\eta_{K_c,\text{pred}} = \eta_{K} \,
     *            \gamma_\text{h} \, 0.5^{p_K} \, 0.5^{\text{dim}}
     *            \quad \forall K_c \text{ children of } K$
     *       <td rowspan="2">$\gamma_\text{h} \in (0,\infty)$
     *   <tr><td>h-coarsening
     *       <td>$\eta_{K,\text{pred}} = \sum\limits_{K_c} \eta_{K_c} /
     *            (\gamma_\text{h} \, 0.5^{p_{K_c}})
     *            \quad \forall K_c \text{ children of } K$
     * </table>
     *
     * With these predicted error estimates, we are capable of adapting the
     * finite element on cells based on their refinement history or rather the
     * predicted change of their error estimates.
     *
     * If a cell is flagged for adaptation, we want to perform p-adaptation once
     * the associated error indicators $\eta_{K}$ on cell $K$ satisfy
     * $\eta_{K} < \eta_{K,\text{pred}}$, where the subscript $\text{pred}$
     * denotes the predicted error. This corresponds to our assumption of
     * smoothness being correct, else h-adaptation is applied. We achieve this
     * with the function hp::Refinement::p_adaptivity_from_criteria() and a
     * function object `std::less<Number>()` for both comparator parameters.
     *
     * For the very first adaptation step, the user needs to decide whether h-
     * or p-adaptation is supposed to happen. An h-step will be applied with
     * $\eta_{K,\text{pred}} = 0$, whereas $\eta_{K,\text{pred}} = \infty$
     * ensures a p-step. The latter may be realised with
     * `std::numeric_limits::infinity()`.
     *
     * The following code snippet demonstrates how to impose hp-adaptivity based
     * on refinement history in an application:
     * @code
     * // [initialisation...]
     * Vector<float> predicted_error_per_cell(triangulation.n_active_cells());
     * for(unsigned int i = 0; i < triangulation.n_active_cells(); ++i)
     *   predicted_error_per_cell[i] = std::numeric_limits<float>::max();
     *
     * // [during each refinement step...]
     * // set h-adaptivity flags
     * Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
     * KellyErrorEstimator::estimate(...);
     * GridRefinemet::refine_and_coarsen_fixed_fraction(...);
     *
     * // set p-adaptivity flags
     * hp::Refinement::p_adaptivity_from_reference(
     *   hp_dof_handler,
     *   estimated_error_per_cell,
     *   predicted_error_per_cell,
     *   std::less<float>(),
     *   std::less<float>());
     * hp::Refinement::{choose|force}_p_over_h(...);
     *
     * // predict error for the subsequent adaptation
     * triangulation.prepare_coarsening_and_refinement();
     * hp::Refinement::predict_error(
     *   hp_dof_handler,
     *   estimated_error_per_cell,
     *   predicted_error_per_cell);
     *
     * // perform adaptation
     * CellDataTransfer<dim, spacedim, Vector<float>>
     *   cell_data_transfer(triangulation);
     * cell_data_transfer.prepare_for_coarsening_and_refinement();
     *
     * triangulation.execute_coarsening_and_refinement();
     *
     * Vector<float> transferred_errors(triangulation.n_active_cells());
     * cell_data_transfer.unpack(predicted_error_per_cell, transferred_errors);
     * predicted_error_per_cell = std::move(transferred_errors);
     * @endcode
     *
     * For more theoretical details see
     * @code{.bib}
     * @article{Melenk2001,
     *  author    = {Melenk, Jens Markus and Wohlmuth, Barbara I.},
     *  title     = {{On residual-based a posteriori error estimation
     *                in hp-FEM}},
     *  journal   = {{Advances in Computational Mathematics}},
     *  volume    = {15},
     *  number    = {1},
     *  pages     = {311--331},
     *  publisher = {Springer US},
     *  year      = {2001},
     *  doi       = {10.1023/A:1014268310921}
     * }
     * @endcode
     *
     * @note This feature is currently only implemented for isotropic refinement.
     *
     * @note We want to predict the error by how adaptation will actually happen.
     *   Thus, this function needs to be called after
     *   Triangulation::prepare_for_coarsening_and_refinement().
     */
    template <int dim, typename Number, int spacedim>
    void
    predict_error(const hp::DoFHandler<dim, spacedim> &dof_handler,
                  const Vector<Number> &               error_indicators,
                  Vector<Number> &                     predicted_errors,
                  const double                         gamma_p = std::sqrt(0.1),
                  const double                         gamma_h = 1.,
                  const double                         gamma_n = 1.);

    /**
     * @}
     */

    /**
     * @name Decide between h- and p-adaptivity
     * @{
     */

    /**
     * Choose p-adaptivity over h-adaptivity in any case.
     *
     * Removes all refine and coarsen flags on cells that have a
     * @p future_fe_index assigned.
     *
     * @note Preceeding calls of Triangulation::prepare_for_coarsening_and_refinement()
     *   may change refine and coarsen flags, which will ultimately change the
     *   results of this function.
     */
    template <int dim, int spacedim>
    void
    force_p_over_h(const hp::DoFHandler<dim, spacedim> &dof_handler);

    /**
     * Choose p-adaptivity over h-adaptivity whenever it is invoked on all
     * related cells.
     *
     * In case of refinement, information about finite elements will be
     * inherited. Thus we will prefer p-refinement over h-refinement whenever
     * desired, i.e. clear the refine flag and supply a corresponding
     * @p future_fe_index.
     *
     * However for coarsening, we follow a different approach. Flagging a cell
     * for h-coarsening does not ultimately mean that it will be coarsened. Only
     * if a cell and all of its siblings are flagged, they will be merged into
     * their parent cell. If we consider p-coarsening on top, we must decide for
     * all siblings together how they will be coarsened. We distinguish between
     * three different cases:
     * <ol>
     * <li> Not all siblings flagged for coarsening: p-coarsening<br>
     *   We keep the @p future_fe_indices and clear the coarsen flags
     *   on all siblings.
     * <li> All siblings flagged for coarsening, but not all for
     *   p-adaptation: h-coarsening<br>
     *   We keep the coarsen flags and clear all @p future_fe_indices
     *   on all siblings.
     * <li> All siblings flagged for coarsening and p-adaptation: p-coarsening<br>
     *   We keep the @p future_fe_indices and clear the coarsen flags
     *   on all siblings.
     * </ol>
     *
     * @note The function Triangulation::prepare_coarsening_and_refinement()
     *   will clean up all h-coarsening flags if they are not shared among
     *   all siblings. In the hp case, we need to bring forward this decision:
     *   If the cell will not be coarsened, but qualifies for p-adaptivity,
     *   we have to set all flags accordingly. So this function anticipates
     *   the decision that Triangulation::prepare_coarsening_and_refinement()
     *   would have made later on.
     *
     * @note Preceeding calls of Triangulation::prepare_for_coarsening_and_refinement()
     *   may change refine and coarsen flags, which will ultimately change the
     *   results of this function.
     */
    template <int dim, int spacedim>
    void
    choose_p_over_h(const hp::DoFHandler<dim, spacedim> &dof_handler);

    /**
     * @}
     */
  } // namespace Refinement
} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_hp_refinement_h
