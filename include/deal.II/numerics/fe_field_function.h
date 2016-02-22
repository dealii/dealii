// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__fe_function_h
#define dealii__fe_function_h

#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_local_storage.h>

#include <deal.II/lac/vector.h>

#include <boost/optional.hpp>


DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  class ExcPointNotAvailableHere;
}

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
   * things a bit, by calling set_active_cell before asking for values or
   * gradients of the function. If you don't do this, and your points don't
   * lie in the cell that is currently stored, the function
   * GridTools::find_cell_around_point is called to find out where the point
   * is. You can specify an optional mapping to use when looking for points in
   * the grid. If you don't do so, this function uses a Q1 mapping.
   *
   * Once the FEFieldFunction knows where the points lie, it creates a
   * quadrature formula for those points, and calls
   * FEValues::get_function_values or FEValues::get_function_gradients with
   * the given quadrature points.
   *
   * If you only need the quadrature points but not the values of the finite
   * element function (you might want this for the adjoint interpolation), you
   * can also use the function @p compute_point_locations alone.
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
   * // Do the same with the second
   * DoFHandler<dim> dh2 (tria_2);
   * Vector<double> solution_2;
   *
   * // Setup the system, assemble matrices, solve problems and get the
   * // nobel prize on the first domain...
   *
   * // Now project it to the second domain
   * FEFieldFunction<dim> fe_function_1 (dh_1, solution_1);
   * VectorTools::project (dh_2, constraints_2, quad, fe_function_1, solution_2);
   *
   * // Or interpolate it...
   * Vector<double> solution_3;
   * VectorTools::interpolate (dh_2, fe_function_1, solution_3);
   *
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
   * <h3>Using FEFieldFunction with parallel::distributed::Triangulation</h3>
   *
   * When using this class with a parallel distributed triangulation object
   * and evaluating the solution at a particular point, not every processor
   * will own the cell at which the solution is evaluated. Rather, it may be
   * that the cell in which this point is found is in fact a ghost or
   * artificial cell (see
   * @ref GlossArtificialCell
   * and
   * @ref GlossGhostCell).
   * If the cell is artificial, we have no access to the solution there and
   * functions that evaluate the solution at such a point will trigger an
   * exception of type VectorTools::ExcPointNotAvailableHere. The same kind of
   * exception will also be produced if the cell is a ghost cell: On such
   * cells, one could in principle evaluate the solution, but it becomes
   * easier if we do not allow to do so because then there is exactly one
   * processor in a parallel distributed computation that can indeed evaluate
   * the solution. Consequently, it is clear which processor is responsible
   * for producing output if the point evaluation is done as a postprocessing
   * step.
   *
   * To deal with this situation, you will want to use code as follows when,
   * for example, evaluating the solution at the origin (here using a parallel
   * TrilinosWrappers vector to hold the solution):
   * @code
   *   Functions::FEFieldFunction<dim,DoFHandler<dim>,TrilinosWrappers::MPI::Vector>
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
   *
   * @ingroup functions
   * @author Luca Heltai, 2006, Markus Buerg, 2012, Wolfgang Bangerth, 2013
   */
  template <int dim,
            typename DoFHandlerType=DoFHandler<dim>,
            typename VectorType=Vector<double> >
  class FEFieldFunction :  public Function<dim, typename VectorType::value_type>
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
    FEFieldFunction (const DoFHandlerType &dh,
                     const VectorType     &data_vector,
                     const Mapping<dim>   &mapping = StaticMappingQ1<dim>::mapping);

    /**
     * Set the current cell. If you know in advance where your points lie, you
     * can tell this object by calling this function. This will speed things
     * up a little.
     */
    void set_active_cell (const typename DoFHandlerType::active_cell_iterator &newcell);

    /**
     * Get one vector value at the given point. It is inefficient to use
     * single points. If you need more than one at a time, use the
     * vector_value_list() function. For efficiency reasons, it is better if
     * all the points lie on the same cell. This is not mandatory, however it
     * does speed things up.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void vector_value (const Point<dim> &p,
                               Vector<typename VectorType::value_type>   &values) const;

    /**
     * Return the value of the function at the given point. Unless there is
     * only one component (i.e. the function is scalar), you should state the
     * component you want to have evaluated; it defaults to zero, i.e. the
     * first component. It is inefficient to use single points. If you need
     * more than one at a time, use the vector_value_list function. For
     * efficiency reasons, it is better if all the points lie on the same
     * cell. This is not mandatory, however it does speed things up.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual typename VectorType::value_type value (const Point< dim > &p,
                                                   const unsigned int  component = 0)    const;

    /**
     * Set @p values to the point values of the specified component of the
     * function at the @p points. It is assumed that @p values already has the
     * right size, i.e. the same size as the points array. This is rather
     * efficient if all the points lie on the same cell. If this is not the
     * case, things may slow down a bit.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void value_list (const std::vector<Point< dim > >     &points,
                             std::vector<typename VectorType::value_type > &values,
                             const unsigned int  component = 0)    const;


    /**
     * Set @p values to the point values of the function at the @p points. It
     * is assumed that @p values already has the right size, i.e. the same
     * size as the points array. This is rather efficient if all the points
     * lie on the same cell. If this is not the case, things may slow down a
     * bit.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void vector_value_list (const std::vector<Point< dim > >     &points,
                                    std::vector<Vector<typename VectorType::value_type> > &values) const;

    /**
     * Return the gradient of all components of the function at the given
     * point.  It is inefficient to use single points. If you need more than
     * one at a time, use the vector_value_list function. For efficiency
     * reasons, it is better if all the points lie on the same cell. This is
     * not mandatory, however it does speed things up.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_gradient (const Point< dim > &p,
                     std::vector< Tensor< 1, dim,typename VectorType::value_type > > &gradients) const;

    /**
     * Return the gradient of the specified component of the function at the
     * given point. It is inefficient to use single points. If you need more
     * than one at a time, use the vector_value_list function. For efficiency
     * reasons, it is better if all the points lie on the same cell. This is
     * not mandatory, however it does speed things up.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual Tensor<1,dim,typename VectorType::value_type> gradient(const Point< dim > &p,
        const unsigned int component = 0)const;

    /**
     * Return the gradient of all components of the function at all the given
     * points. This is rather efficient if all the points lie on the same
     * cell. If this is not the case, things may slow down a bit.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_gradient_list (const std::vector< Point< dim > > &p,
                          std::vector<
                          std::vector< Tensor< 1, dim,typename VectorType::value_type > > > &gradients) const;

    /**
     * Return the gradient of the specified component of the function at all
     * the given points.  This is rather efficient if all the points lie on
     * the same cell. If this is not the case, things may slow down a bit.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    gradient_list (const std::vector< Point< dim > > &p,
                   std::vector< Tensor< 1, dim,typename VectorType::value_type > > &gradients,
                   const unsigned int component=0) const;


    /**
     * Compute the Laplacian of a given component at point <tt>p</tt>.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual typename VectorType::value_type
    laplacian (const Point<dim>   &p,
               const unsigned int  component = 0) const;

    /**
     * Compute the Laplacian of all components at point <tt>p</tt> and store
     * them in <tt>values</tt>.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_laplacian (const Point<dim>   &p,
                      Vector<typename VectorType::value_type>     &values) const;

    /**
     * Compute the Laplacian of one component at a set of points.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    laplacian_list (const std::vector<Point<dim> > &points,
                    std::vector<typename VectorType::value_type>            &values,
                    const unsigned int              component = 0) const;

    /**
     * Compute the Laplacians of all components at a set of points.
     *
     * @note When using this function on a
     * parallel::distributed::Triangulation you may get an exception when
     * trying to evaluate the solution at a point that does not lie in a
     * locally owned cell (see
     * @ref GlossLocallyOwnedCell).
     * See the section in the general documentation of this class for more
     * information.
     */
    virtual void
    vector_laplacian_list (const std::vector<Point<dim> > &points,
                           std::vector<Vector<typename VectorType::value_type> >   &values) const;

    /**
     * Create quadrature rules. This function groups the points into blocks
     * that live in the same cell, and fills up three vectors: @p cells, @p
     * qpoints and @p maps. The first is a list of the cells that contain the
     * points, the second is a list of quadrature points matching each cell of
     * the first list, and the third contains the index of the given
     * quadrature points, i.e., @p points[maps[3][4]] ends up as the 5th
     * quadrature point in the 4th cell. This is where optimization would
     * help. This function returns the number of cells that contain the given
     * set of points.
     */
    unsigned int
    compute_point_locations
    (const std::vector<Point<dim> >                              &points,
     std::vector<typename DoFHandlerType::active_cell_iterator > &cells,
     std::vector<std::vector<Point<dim> > >                      &qpoints,
     std::vector<std::vector<unsigned int> >                     &maps) const;

    /**
     * @deprecated Use VectorTools::ExcPointNotAvailableHere instead.
     */
    typedef VectorTools::ExcPointNotAvailableHere ExcPointNotAvailableHere DEAL_II_DEPRECATED;

  private:
    /**
     * Typedef holding the local cell_hint.
     */
    typedef
    Threads::ThreadLocalStorage <typename DoFHandlerType::active_cell_iterator >
    cell_hint_t;

    /**
     * Pointer to the dof handler.
     */
    SmartPointer<const DoFHandlerType,FEFieldFunction<dim, DoFHandlerType, VectorType> > dh;

    /**
     * A reference to the actual data vector.
     */
    const VectorType &data_vector;

    /**
     * A reference to the mapping being used.
     */
    const Mapping<dim> &mapping;

    /**
     * The latest cell hint.
     */
    mutable cell_hint_t cell_hint;

    /**
     * Store the number of components of this function.
     */
    const unsigned int n_components;

    /**
     * Given a cell, return the reference coordinates of the given point
     * within this cell if it indeed lies within the cell. Otherwise return an
     * uninitialized boost::optional object.
     */
    boost::optional<Point<dim> >
    get_reference_coordinates (const typename DoFHandlerType::active_cell_iterator &cell,
                               const Point<dim>                                    &point) const;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif
