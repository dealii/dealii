// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_function_h
#define dealii_fe_function_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_local_storage.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/vector.h>

#include <optional>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
namespace VectorTools
{
  class ExcPointNotAvailableHere;
}
#endif

namespace Functions
{
  /**
   * This is an interpolation function for the given dof handler and the given
   * solution vector. The points at which this function can be evaluated MUST
   * be inside the domain of the dof handler, but except from this, no other
   * requirement is given. This function is rather slow, as it needs to
   * construct a quadrature object for the point (or set of points) where you
   * want to evaluate your finite element function. In order to do so, it
   * needs to find out where the points lie.
   *
   * If you know in advance in which cell your points lie, you can accelerate
   * things a bit, by calling set_active_cell() before asking for values or
   * gradients of the function. If you don't do this, and your points don't
   * lie in the cell that is currently stored, the function
   * GridTools::find_active_cell_around_point is called to find out where the
   * point is.
   * You can specify an optional mapping to use when looking for points in
   * the grid. If you don't do so, this function uses a Q1 mapping.
   *
   * Once the FEFieldFunction knows where the points lie, it creates a
   * quadrature formula for those points, and calls
   * FEValues::get_function_values or FEValues::get_function_gradients with
   * the given quadrature points.
   *
   * If you only need the quadrature points but not the values of the finite
   * element function (you might want this for the adjoint interpolation), you
   * can also use the function compute_point_locations() alone.
   *
   * An example of how to use this function is the following:
   *
   * @code
   *
   * // Generate two triangulations
   * Triangulation<dim> tria_1;
   * Triangulation<dim> tria_2;
   *
   * // Read the triangulations from files, or build them up, or get them
   * // from some place. Assume that tria_2 is *entirely* included in tria_1.
   *
   * // Associate a dof handler and a solution to the first triangulation
   * DoFHandler<dim> dh1 (tria_1);
   * Vector<double> solution_1;
   *
   * // On this first domain, set up the various data structures,
   * // assemble matrices, solve the linear system, and get a Nobel
   * // prize for the work we have done here:
   * [...]
   *
   * // Then create a DoFHandler and solution vector for the second domain:
   * DoFHandler<dim> dh2 (tria_2);
   * Vector<double> solution_2;
   *
   * // Finally, project the solution on the first domain onto the
   * // second domain, assuming that this does not require querying
   * // values from outside the first domain:
   * Functions::FEFieldFunction<dim> fe_function_1 (dh_1, solution_1);
   * VectorTools::project (dh_2, constraints_2, quad,
   *                       fe_function_1, solution_2);
   *
   * // Alternatively, we could have also interpolated it:
   * Vector<double> solution_3;
   * VectorTools::interpolate (dh_2, fe_function_1, solution_3);
   * @endcode
   *
   * The snippet of code above will work assuming that the second
   * triangulation is entirely included in the first one.
   *
   * FEFieldFunction is designed to be an easy way to get the results of your
   * computations across different, possibly non matching, grids. No knowledge
   * of the location of the points is assumed in this class, which makes it
   * rely entirely on the GridTools::find_active_cell_around_point utility for
   * its job. However the class can be fed an "educated guess" of where the
   * points that will be computed actually are by using the
   * FEFieldFunction::set_active_cell method, so if you have a smart way to
   * tell where your points are, you will save a lot of computational time by
   * letting this class know.
   *
   *
   * <h3>Using FEFieldFunction with parallel triangulations</h3>
   *
   * When using this class with a parallel distributed triangulation object
   * and evaluating the solution at a particular point, not every processor
   * will own the cell at which the solution should be evaluated. Rather, it may
   * be that the cell in which this point is found is in fact a ghost or
   * artificial cell (see
   * @ref GlossArtificialCell
   * and
   * @ref GlossGhostCell).
   * The solution can be evaluated on ghost cells, but for artificial cells
   * we have no access to the solution there and
   * functions that evaluate the solution at such a point will trigger an
   * exception of type VectorTools::ExcPointNotAvailableHere.
   *
   * To deal with this situation, you will want to use code as follows when,
   * for example, evaluating the solution at the origin (here using a parallel
   * TrilinosWrappers vector to hold the solution):
   * @code
   *   Functions::FEFieldFunction<dim,TrilinosWrappers::MPI::Vector>
   *     solution_function (dof_handler, solution);
   *   Point<dim> origin = Point<dim>();
   *
   *   double solution_at_origin;
   *   bool   point_found = true;
   *   try
   *     {
   *       solution_at_origin = solution_function.value (origin);
   *     }
   *   catch (const VectorTools::ExcPointNotAvailableHere &)
   *     {
   *       point_found = false;
   *     }
   *
   *   if (point_found == true)
   *     ...do something...;
   * @endcode
   * The last part (the `...do something...`) part needs to be application
   * specific. If you really do need the information on this process, you need
   * to communicate with the process that actually owns the cell in which the
   * point at which you want the solution is located. Since you don't know
   * which process that may be, this will necessarily require a global
   * communication step. The FERemoteEvaluation class may be useful in this
   * context.
   *
   * Finally, there is the case of the value_list(), gradient_list(),
   * and similar functions of this class that take a list of points. For these
   * functions, it may be that some of the points are on locally owned cells,
   * and others are not -- in which case, again, you will receive an exception.
   * Because these functions do not return information about which of the points
   * can and which cannot be evaluated, the correct approach is to call the
   * `*_list()` function within a `try` block, and if it fails fall back to
   * calling value(), gradient(), or similar in a loop over the points you
   * wanted evaluated to find out for which point the call fails.
   *
   * @ingroup functions
   */
  template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
  class FEFieldFunction : public Function<dim, typename VectorType::value_type>
  {
  public:
    /**
     * Construct a vector function. A smart pointers is stored to the dof
     * handler, so you have to make sure that it make sense for the entire
     * lifetime of this object. The number of components of this functions is
     * equal to the number of components of the finite element object. If a
     * mapping is specified, that is what is used to find out where the points
     * lay. Otherwise the standard Q1 mapping is used.
     */
    FEFieldFunction(
      const DoFHandler<dim, spacedim> &dh,
      const VectorType                &data_vector,
      const Mapping<dim>              &mapping = StaticMappingQ1<dim>::mapping);

    /**
     * Set the current cell. If you know in advance where your points lie, you
     * can tell this object by calling this function. This will speed things
     * up a little.
     */
    void
    set_active_cell(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &newcell);

    /**
     * Get one vector value at the given point. It is inefficient to use
     * single points. If you need more than one at a time, use the
     * vector_value_list() function. For efficiency reasons, it is better if
     * all the points lie on the same cell. This is not mandatory, however it
     * does speed things up.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_value(
      const Point<dim>                        &p,
      Vector<typename VectorType::value_type> &values) const override;

    /**
     * Return the value of the function at the given point. Unless there is
     * only one component (i.e. the function is scalar), you should state the
     * component you want to have evaluated; it defaults to zero, i.e. the
     * first component. It is inefficient to use single points. If you need
     * more than one at a time, use the vector_value_list() function. For
     * efficiency reasons, it is better if all the points lie on the same
     * cell. This is not mandatory, however it does speed things up.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual typename VectorType::value_type
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * Set @p values to the point values of the specified component of the
     * function at the @p points. It is assumed that @p values already has the
     * right size, i.e. the same size as the points array. This is rather
     * efficient if all the points lie on the same cell. If this is not the
     * case, things may slow down a bit.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    value_list(const std::vector<Point<dim>>                &points,
               std::vector<typename VectorType::value_type> &values,
               const unsigned int component = 0) const override;


    /**
     * Set @p values to the point values of the function at the @p points. It
     * is assumed that @p values already has the right size, i.e. the same
     * size as the points array. This is rather efficient if all the points
     * lie on the same cell. If this is not the case, things may slow down a
     * bit.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<typename VectorType::value_type>>
                        &values) const override;

    /**
     * Return the gradient of all components of the function at the given
     * point.  It is inefficient to use single points. If you need more than
     * one at a time, use the vector_value_list() function. For efficiency
     * reasons, it is better if all the points lie on the same cell. This is
     * not mandatory, however it does speed things up.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_gradient(const Point<dim> &p,
                    std::vector<Tensor<1, dim, typename VectorType::value_type>>
                      &gradients) const override;

    /**
     * Return the gradient of the specified component of the function at the
     * given point. It is inefficient to use single points. If you need more
     * than one at a time, use the vector_value_list() function. For efficiency
     * reasons, it is better if all the points lie on the same cell. This is
     * not mandatory, however it does speed things up.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual Tensor<1, dim, typename VectorType::value_type>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;

    /**
     * Return the gradient of all components of the function at all the given
     * points. This is rather efficient if all the points lie on the same
     * cell. If this is not the case, things may slow down a bit.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &p,
      std::vector<std::vector<Tensor<1, dim, typename VectorType::value_type>>>
        &gradients) const override;

    /**
     * Return the gradient of the specified component of the function at all
     * the given points.  This is rather efficient if all the points lie on
     * the same cell. If this is not the case, things may slow down a bit.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    gradient_list(
      const std::vector<Point<dim>>                                &p,
      std::vector<Tensor<1, dim, typename VectorType::value_type>> &gradients,
      const unsigned int component = 0) const override;


    /**
     * Compute the Laplacian of a given component at point <tt>p</tt>.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual typename VectorType::value_type
    laplacian(const Point<dim>  &p,
              const unsigned int component = 0) const override;

    /**
     * Compute the Laplacian of all components at point <tt>p</tt> and store
     * them in <tt>values</tt>.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_laplacian(
      const Point<dim>                        &p,
      Vector<typename VectorType::value_type> &values) const override;

    /**
     * Compute the Laplacian of one component at a set of points.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>>                &points,
                   std::vector<typename VectorType::value_type> &values,
                   const unsigned int component = 0) const override;

    /**
     * Compute the Laplacians of all components at a set of points.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that lies on an artificial
     * cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_laplacian_list(const std::vector<Point<dim>> &points,
                          std::vector<Vector<typename VectorType::value_type>>
                            &values) const override;

    /**
     * Given a set of points located in the domain (or, in the case of
     * a parallel Triangulation, in the locally owned part of the domain
     * or on the ghost cells for the current processor), sort these
     * points into buckets for each of the cells on which at least
     * one of the points is located.
     *
     * This function fills three output vectors: @p cells, @p qpoints
     * and @p maps. The first is a list of the cells that contain the
     * points, the second is a list of quadrature points matching each
     * cell of the first list, and the third contains the index of the
     * given quadrature points, i.e., @p points[maps[3][4]] ends up as
     * the 5th quadrature point in the 4th cell.
     *
     * @return This function returns the number of cells that
     *   collectively contain the set of points give as @p
     *   points. This also equals the lengths of the output arrays.
     *
     * This function simply calls GridTools::compute_point_locations :
     * using the original function avoids computing a
     * new Cache at every function call.
     */
    unsigned int
    compute_point_locations(
      const std::vector<Point<dim>> &points,
      std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
                                             &cells,
      std::vector<std::vector<Point<dim>>>   &qpoints,
      std::vector<std::vector<unsigned int>> &maps) const;

  private:
    /**
     * Typedef holding the local cell_hint.
     */
    using cell_hint_t = Threads::ThreadLocalStorage<
      typename DoFHandler<dim, spacedim>::active_cell_iterator>;

    /**
     * Pointer to the dof handler.
     */
    ObserverPointer<const DoFHandler<dim, spacedim>,
                    FEFieldFunction<dim, VectorType, spacedim>>
      dh;

    /**
     * A reference to the actual data vector.
     */
    const VectorType &data_vector;

    /**
     * A reference to the mapping being used.
     */
    const Mapping<dim> &mapping;

    /**
     * The Cache object
     */
    GridTools::Cache<dim, spacedim> cache;

    /**
     * The latest cell hint.
     */
    mutable cell_hint_t cell_hint;

    /**
     * Given a cell, return the reference coordinates of the given point
     * within this cell if it indeed lies within the cell. Otherwise return an
     * uninitialized std::optional object.
     */
    std::optional<Point<dim>>
    get_reference_coordinates(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      const Point<dim> &point) const;
  };
} // namespace Functions


DEAL_II_NAMESPACE_CLOSE

#endif
