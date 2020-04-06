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

#ifndef dealii_fe_field_convolution_function_h
#define dealii_fe_field_convolution_function_h

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/scratch_data.h>


DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * This is an approximated interpolation function for the given dof handler
   * and the given solution vector, that actually performs a convolution with
   * a specified kernel in order to compute its values.
   *
   * In particular, the evaluation of this function at a point $x$ is computed
   * as:
   *
   * \f[
   * F(x) = \int_{\Omega} f(y) K^{\epsilon}(x-y) dy
   * \f]
   * where $f(y)$ is the finite element field function constructed using the
   * vector and the dof handler passed at construction time.
   *
   * This class also takes a CutOffFunctionBase object in input, that is used to
   * perform the convolution, and a radius, that is used to scale the
   * convolution kernel. This function may be called also for points
   * outside the domain, where it decays to zero according to what cut off
   * function you used to construct this object.
   *
   * An example usage is the following:
   * @code
   *
   * // Generate two triangulations
   * Triangulation<dim> tria_1;
   * Triangulation<dim> tria_2;
   *
   * // Read the triangulations from files, or build them up, or get them
   * // from some place.
   * ...
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
   * ...
   *
   * // Create a convolution kernel, and a convolution radius
   * Functions::CutOffFunctionC1<dim> kernel;
   * const double radius = 2*h;
   *
   * // Now build the convolution of the first function, projected on the
   * // second domain
   * FEFieldConvolutionFunction<dim> fe_function_1 (dh_1, solution_1,
   *                                                kernel, radius);
   * VectorTools::project (dh_2, constraints_2, quad,
   *                       fe_function_1, solution_2);
   *
   * // Or interpolated...
   * Vector<double> solution_3;
   * VectorTools::interpolate (dh_2, fe_function_1, solution_3);
   *
   * @endcode
   *
   * FEFieldConvolutionFunction is designed to be an easy way to get the results
   * of your computations across different, possibly non matching, grids. No
   * knowledge of the location of the points is assumed in this class.
   *
   * <h3>Using FEFieldConvolutionFunction with
   * parallel::distributed::Triangulation</h3>
   *
   * When using this class with a parallel distributed triangulation object
   * and evaluating the solution at a particular point, not every processor
   * will own the cell at which the solution is evaluated. Rather, it may be
   * that the cell in which this point is found is in fact a ghost or
   * artificial cell (see
   * @ref GlossArtificialCell
   * and
   * @ref GlossGhostCell).
   * The solution may be evaluated on non owned cells, but there the solution
   * will vanish at a distance greater than `radius` from the locally owned
   * domain. If you want to have a consistent solution on distributed
   * Triangulation objects, you have to make sure you sum the results from all
   * processors, i.e., using Utilities::MPI::sum(). In this way you are
   * guaranteed that the actual convolution is performed with a non-zero field
   * on only one processor, and the sum of the results across all processors
   * will produce the same value you'd get with a serial Triangulation.
   *
   * @ingroup functions
   * @author Luca Heltai, 2019
   */
  template <int dim,
            int spacedim            = dim,
            typename DoFHandlerType = DoFHandler<dim, spacedim>,
            typename VectorType     = Vector<double>>
  class FEFieldConvolutionFunction
    : public Function<spacedim, typename VectorType::value_type>
  {
  public:
    /**
     * Construct a (possibly vector-valued) function. A SmartPointer is stored to the dof
     * handler and to the cache object, so you have to make sure they remain
     * valid for the entire lifetime of this object. The number of
     * components of this functions is equal to the number of components of the
     * finite element object.
     */
    FEFieldConvolutionFunction(const DoFHandlerType &                 dh,
                               const GridTools::Cache<dim, spacedim> &cache,
                               const VectorType &data_vector,
                               Functions::CutOffFunctionBase<spacedim> &kernel,
                               const Quadrature<dim> &quadrature);

    /**
     * Get one vector value at the given point. It is inefficient to use
     * single points. If you need more than one at a time, use the
     * vector_value_list() function.
     */
    virtual void
    vector_value(
      const Point<spacedim> &                  p,
      Vector<typename VectorType::value_type> &values) const override;

    /**
     * Return the value of the function at the given point. Unless there is
     * only one component (i.e. the function is scalar), you should state the
     * component you want to have evaluated; it defaults to zero, i.e. the
     * first component. It is inefficient to use single points. If you need
     * more than one at a time, use the vector_value_list function.
     */
    virtual typename VectorType::value_type
    value(const Point<spacedim> &p,
          const unsigned int     component = 0) const override;

    /**
     * Set @p values to the point values of the specified component of the
     * function at the @p points. It is assumed that @p values already has the
     * right size, i.e. the same size as the points array.
     */
    virtual void
    value_list(const std::vector<Point<spacedim>> &          points,
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
    vector_value_list(const std::vector<Point<spacedim>> &points,
                      std::vector<Vector<typename VectorType::value_type>>
                        &values) const override;

    /**
     * Return the gradient of all components of the function at the given
     * point. Notice that by the definition of convolution, this will return
     * the convolution with the gradient of the kernel function.
     */
    virtual void
    vector_gradient(
      const Point<spacedim> &p,
      std::vector<Tensor<1, spacedim, typename VectorType::value_type>>
        &gradients) const override;

    /**
     * Return the gradient of the specified component of the function at the
     * given point. Notice that by the definition of convolution, this will
     * return the convolution with the gradient of the kernel function.
     */
    virtual Tensor<1, spacedim, typename VectorType::value_type>
    gradient(const Point<spacedim> &p,
             const unsigned int     component = 0) const override;

    /**
     * Return the gradient of all components of the function at all the given
     * points. Notice that by the definition of convolution, this will
     * return the convolution with the gradient of the kernel function.
     */
    virtual void
    vector_gradient_list(
      const std::vector<Point<spacedim>> &p,
      std::vector<
        std::vector<Tensor<1, spacedim, typename VectorType::value_type>>>
        &gradients) const override;

    /**
     * Return the gradient of the specified component of the function at all
     * the given points. Notice that by the definition of convolution, this will
     * return the convolution with the gradient of the kernel function.
     */
    virtual void
    gradient_list(
      const std::vector<Point<spacedim>> &p,
      std::vector<Tensor<1, spacedim, typename VectorType::value_type>>
        &                gradients,
      const unsigned int component = 0) const override;

    /**
     * Set the convolution radius to this new value.
     */
    void
    set_radius(const double &value);

  private:
    /**
     * Pointer to the dof handler.
     */
    const SmartPointer<const DoFHandlerType> dh;

    /**
     * Pointer to the Cache object
     */
    const SmartPointer<const GridTools::Cache<dim, spacedim>> cache;

    /**
     * A reference to the actual data vector.
     */
    const VectorType &data_vector;

    /**
     * Pointer to the kernel function.
     */
    const SmartPointer<Functions::CutOffFunctionBase<spacedim>> kernel;

    /**
     * ScratchData used to perform evaluations.
     */
    mutable MeshWorker::ScratchData<dim, spacedim> scratch;
  };
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
