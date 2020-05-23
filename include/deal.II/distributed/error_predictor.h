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

#ifndef dealii_distributed_error_predictor_h
#define dealii_distributed_error_predictor_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <boost/range/iterator_range.hpp>

#include <vector>


DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename Number>
class Vector;

namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    class Triangulation;
  }
} // namespace parallel

namespace hp
{
  template <int dim, int spacedim>
  class DoFHandler;
}
#endif


namespace parallel
{
  namespace distributed
  {
    /**
     * Predict how current error indicators will change after refinement and
     * coarsening were to happen on the provided hp::DoFHandler in context of a
     * parallel::distributed::Triangulation.
     *
     * This algorithm follows the same logic as the error prediction algorithm
     * presented by @cite melenk2001hp and has been implemented in deal.II via
     * the function hp::Refinement::predict_error(). See the documentation of
     * this particular function for details on the algorithm and all featured
     * formulas.
     *
     * In theory, the combination of hp::Refinement::predict_error() to predict
     * how the error will change with  parallel::distributed::CellDataTransfer
     * to transfer data across meshes forms the intended way to apply the
     * prediction algorithm. However, p4est determines the details of grid
     * refinement in the parallel distributed case; consequently, it yields more
     * reliable and trustworthy results when we determine the predicted errors
     * during the adaptation process. This is achieved with this class, similar
     * to the parallel::distributed::SolutionTransfer and
     * parallel::distributed::CellDataTransfer classes.
     *
     * Before adaptation happens, one or multiple vectors containing error
     * estimates for every cell have to be registered with this class using the
     * prepare_for_coarsening_and_refinement() functions. Performing adaptation
     * on the triangulation initiates the calculation of the error prediction
     * based on how p4est performed grid adaptation. After the refinement
     * process, the predicted errors need to be unpacked on the updated mesh via
     * one of the unpack() functions.
     *
     * The following code snippet demonstrates how to impose hp-adaptivity based
     * on refinement history in a parallel distributed application:
     * @code
     * // [initialization...]
     * Vector<float> predicted_error_per_cell(triangulation.n_active_cells());
     * for(unsigned int i = 0; i < triangulation.n_active_cells(); ++i)
     *   predicted_error_per_cell[i] = std::numeric_limits<float>::infinity();
     *
     * // [during each refinement step...]
     * // set h-adaptivity flags
     * Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
     * KellyErrorEstimator<dim>::estimate(...);
     * parallel::distributed::GridRefinement::
     *   refine_and_coarsen_fixed_{number|fraction}(...);
     *
     * // set p-adaptivity flags
     * hp::Refinement::p_adaptivity_from_reference(
     *   hp_dof_handler,
     *   estimated_error_per_cell,
     *   predicted_error_per_cell
     *   std::less<float>(),
     *   std::less<float>());
     * hp::Refinement::{choose|force}_p_over_h(hp_dof_handler);
     *
     * // perform adaptation and predict error for the subsequent adaptation
     * parallel::distributed::ErrorPredictor<dim> predictor(hp_dof_handler);
     * predictor.prepare_for_coarsening_and_refinement(estimated_error_per_cell);
     * triangulation.execute_coarsening_and_refinement();
     *
     * // unpack predicted errors
     * predicted_error_per_cell.reinit(triangulation.n_active_cells());
     * predictor.unpack(predicted_error_per_cell);
     * @endcode
     *
     * For serialization of predicted errors, use the
     * parallel::distributed::CellDataTransfer class.
     *
     * @ingroup distributed
     * @author Marc Fehling, 2019 - 2020
     */
    template <int dim, int spacedim = dim>
    class ErrorPredictor
    {
    public:
      /**
       * Constructor.
       *
       * @param[in] dof The hp::DoFHandler on which all operations will
       *   happen. At the time when this constructor is called, the
       *   DoFHandler still points to the triangulation before the
       *   refinement in question happens.
       */
      ErrorPredictor(const hp::DoFHandler<dim, spacedim> &dof);

      /**
       * Prepare the current object for coarsening and refinement.
       * @p all_in includes all vectors that are to be interpolated onto the
       * new (refined and/or coarsened) grid.
       *
       * As specified in the error prediction algorithm in
       * hp::Refinement::predict_error(), the control parameters @p gamma_p,
       * @p gamma_h, @p gamma_n can be used to influence the predicted errors in
       * case of p-adaptation, h-adaptation, and no adapation, respectively.
       * Their default values comply to the studies of @cite melenk2001hp.
       */
      void
      prepare_for_coarsening_and_refinement(
        const std::vector<const Vector<float> *> &all_in,
        const double                              gamma_p = std::sqrt(0.4),
        const double                              gamma_h = 2.,
        const double                              gamma_n = 1.);

      /**
       * Same as the previous function but for only one discrete function to be
       * interpolated.
       */
      void
      prepare_for_coarsening_and_refinement(
        const Vector<float> &in,
        const double         gamma_p = std::sqrt(0.4),
        const double         gamma_h = 2.,
        const double         gamma_n = 1.);

      /**
       * Interpolate the data previously stored in this object before the mesh
       * was refined or coarsened onto the current set of cells. Do so for
       * each of the vectors provided to
       * prepare_for_coarsening_and_refinement() and write the result into the
       * given set of vectors.
       */
      void
      unpack(std::vector<Vector<float> *> &all_out);

      /**
       * Same as the previous function. It interpolates only one function. It
       * assumes the vectors having the right sizes (i.e.
       * <tt>in.size()==n_cells_old</tt>, <tt>out.size()==n_cells_refined</tt>)
       *
       * Multiple calling of this function is NOT allowed. Interpolating
       * several functions can be performed in one step by using the above
       * function.
       */
      void
      unpack(Vector<float> &out);

    private:
      /**
       * Pointer to the degree of freedom handler to work with.
       */
      SmartPointer<const hp::DoFHandler<dim, spacedim>,
                   ErrorPredictor<dim, spacedim>>
        dof_handler;

      /**
       * A vector that stores pointers to all the vectors we are supposed to
       * copy over from the old to the new mesh.
       */
      std::vector<const Vector<float> *> error_indicators;

      /**
       * The handle that the Triangulation has assigned to this object
       * with which we can access our memory offset and our pack function.
       */
      unsigned int handle;

      /**
       * Control parameters for the prediction algorithm.
       */
      double gamma_p, gamma_h, gamma_n;

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
          &                           data_range,
        std::vector<Vector<float> *> &all_out);


      /**
       * Registers the pack_callback() function to the
       * parallel::distributed::Triangulation that has been assigned to the
       * DoFHandler class member and stores the returning handle.
       */
      void
      register_data_attach();
    };
  } // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif
