// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_values_base_h
#define dealii_fe_values_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/std_cxx20/iota_view.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_values_views.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_related_data.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/read_vector.h>

#include <boost/signals2/connection.hpp>

#include <algorithm>
#include <memory>
#include <optional>
#include <type_traits>
#include <variant>

DEAL_II_NAMESPACE_OPEN

/**
 * FEValues, FEFaceValues and FESubfaceValues objects are interfaces to finite
 * element and mapping classes on the one hand side, to cells and quadrature
 * rules on the other side. They allow to evaluate values or derivatives of
 * shape functions at the quadrature points of a quadrature formula when
 * projected by a mapping from the unit cell onto a cell in real space. The
 * reason for this abstraction is possible optimization: Depending on the type
 * of finite element and mapping, some values can be computed once on the unit
 * cell. Others must be computed on each cell, but maybe computation of
 * several values at the same time offers ways for optimization. Since this
 * interplay may be complex and depends on the actual finite element, it
 * cannot be left to the applications programmer.
 *
 * FEValues, FEFaceValues and FESubfaceValues provide only data handling:
 * computations are left to objects of type Mapping and FiniteElement. These
 * provide functions <tt>get_*_data</tt> and <tt>fill_*_values</tt> which are
 * called by the constructor and <tt>reinit</tt> functions of
 * <tt>FEValues*</tt>, respectively.
 *
 * <h3>General usage</h3>
 *
 * Usually, an object of <tt>FEValues*</tt> is used in integration loops over
 * all cells of a triangulation (or faces of cells). To take full advantage of
 * the optimization features, it should be constructed before the loop so that
 * information that does not depend on the location and shape of cells can be
 * computed once and for all (this includes, for example, the values of shape
 * functions at quadrature points for the most common elements: we can
 * evaluate them on the unit cell and they will be the same when mapped to the
 * real cell). Then, in the loop over all cells, it must be re-initialized for
 * each grid cell to compute that part of the information that changes
 * depending on the actual cell (for example, the gradient of shape functions
 * equals the gradient on the unit cell -- which can be computed once and for
 * all -- times the Jacobian matrix of the mapping between unit and real cell,
 * which needs to be recomputed for each cell).
 *
 * A typical piece of code, adding up local contributions to the Laplace
 * matrix looks like this:
 *
 * @code
 * FEValues values (mapping, finite_element, quadrature, flags);
 * for (const auto &cell : dof_handler.active_cell_iterators())
 *   {
 *     values.reinit(cell);
 *     for (unsigned int q=0; q<quadrature.size(); ++q)
 *       for (unsigned int i=0; i<finite_element.n_dofs_per_cell(); ++i)
 *         for (unsigned int j=0; j<finite_element.n_dofs_per_cell(); ++j)
 *           A(i,j) += fe_values.shape_value(i,q) *
 *                     fe_values.shape_value(j,q) *
 *                     fe_values.JxW(q);
 *     ...
 *   }
 * @endcode
 *
 * The individual functions used here are described below. Note that by
 * design, the order of quadrature points used inside the FEValues object is
 * the same as defined by the quadrature formula passed to the constructor of
 * the FEValues object above.
 *
 * <h3>Member functions</h3>
 *
 * The functions of this class fall into different categories:
 * <ul>
 * <li> shape_value(), shape_grad(), etc: return one of the values of this
 * object at a time. These functions are inlined, so this is the suggested
 * access to all finite element values. There should be no loss in performance
 * with an optimizing compiler. If the finite element is vector valued, then
 * these functions return the only non-zero component of the requested shape
 * function. However, some finite elements have shape functions that have more
 * than one non-zero component (we call them non-"primitive"), and in this
 * case this set of functions will throw an exception since they cannot
 * generate a useful result. Rather, use the next set of functions.
 *
 * <li> shape_value_component(), shape_grad_component(), etc: This is the same
 * set of functions as above, except that for vector valued finite elements
 * they return only one vector component. This is useful for elements of which
 * shape functions have more than one non-zero component, since then the above
 * functions cannot be used, and you have to walk over all (or only the
 * non-zero) components of the shape function using this set of functions.
 *
 * <li> get_function_values(), get_function_gradients(), etc.: Compute a
 * finite element function or its derivative in quadrature points.
 *
 * <li> reinit: initialize the FEValues object for a certain cell. This
 * function is not in the present class but only in the derived classes and
 * has a variable call syntax. See the docs for the derived classes for more
 * information.
 * </ul>
 *
 *
 * <h3>Internals about the implementation</h3>
 *
 * The mechanisms by which this class work are discussed on the page on
 * @ref UpdateFlags "Update flags"
 * and about the
 * @ref FE_vs_Mapping_vs_FEValues "How Mapping, FiniteElement, and FEValues work together".
 *
 *
 * @ingroup feaccess
 */
template <int dim, int spacedim>
class FEValuesBase : public EnableObserverPointer
{
public:
  /**
   * Dimension in which this object operates.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Dimension of the space in which this object operates.
   */
  static constexpr unsigned int space_dimension = spacedim;

  /**
   * Number of quadrature points of the current object. Its value is
   * initialized by the value of max_n_quadrature_points and is updated,
   * e.g., if FEFaceValues::reinit() is called for a new cell/face.
   *
   * @note The default value equals to the value of max_n_quadrature_points.
   */
  const unsigned int n_quadrature_points;

  /**
   * Maximum number of quadrature points. This value might be different from
   * n_quadrature_points, e.g., if a QCollection with different face quadrature
   * rules has been passed to initialize FEFaceValues.
   *
   * This is mostly useful to initialize arrays to allocate the maximum amount
   * of memory that may be used when re-sizing later on to a the current
   * number of quadrature points given by n_quadrature_points.
   */
  const unsigned int max_n_quadrature_points;

  /**
   * Number of shape functions per cell. If we use this base class to evaluate
   * a finite element on faces of cells, this is still the number of degrees
   * of freedom per cell, not per face.
   */
  const unsigned int dofs_per_cell;


  /**
   * Constructor. Set up the array sizes with <tt>n_q_points</tt> quadrature
   * points, <tt>dofs_per_cell</tt> trial functions per cell and with the
   * given pattern to update the fields when the <tt>reinit</tt> function of
   * the derived classes is called. The fields themselves are not set up, this
   * must happen in the constructor of the derived class.
   */
  FEValuesBase(const unsigned int                  n_q_points,
               const unsigned int                  dofs_per_cell,
               const UpdateFlags                   update_flags,
               const Mapping<dim, spacedim>       &mapping,
               const FiniteElement<dim, spacedim> &fe);

  /**
   * The copy assignment is deleted since objects of this class are not
   * copyable.
   */
  FEValuesBase &
  operator=(const FEValuesBase &) = delete;

  /**
   * The copy constructor is deleted since objects of this class are not
   * copyable.
   */
  FEValuesBase(const FEValuesBase &) = delete;

  /**
   * Destructor.
   */
  virtual ~FEValuesBase() override;

  /**
   * Explicitly allow to check for cell similarity.
   * The detection of simple geometries with CellSimilarity is sensitive to the
   * first cell detected. When using multiple threads, each thread might get a
   * thread local copy of the FEValues object that is initialized to the first
   * cell the thread sees. As this cell might be different between runs and
   * number of threads used, this slight deviation leads to difference in
   * roundoff errors that propagate through the program. Therefore, deal.II
   * disables the CellSimilarity check by default in programs that use more than
   * one thread. This function can be used to disable this behavior: When
   * called, FEValues objects will always do the similarity check, even in cases
   * where the program uses multiple threads. This substantially accelerates the
   * operations of FEValues because many operations can be avoided if a cell is,
   * for example, just a translation of the previous cell. On the other hand,
   * you might get results that differ by an amount proportional to round-off
   * between the case of using or not using cell similarity information, and
   * because the order of cells assigned to individual threads may differ from
   * run to run, when you call this function you may end up with results that
   * differ by round-off between runs of the same program.
   */
  void
  always_allow_check_for_cell_similarity(const bool allow);

  /// @name Access to shape function values
  ///
  /// These fields are filled by the finite element.
  /** @{ */

  /**
   * Value of a shape function at a quadrature point on the cell, face or
   * subface selected the last time the <tt>reinit</tt> function of the
   * derived class was called.
   *
   * If the shape function is vector-valued, then this returns the only
   * non-zero component. If the shape function has more than one non-zero
   * component (i.e. it is not primitive), then throw an exception of type
   * ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_value_component() function.
   *
   * @param i Number of the shape function $\varphi_i$ to be evaluated. Note
   * that this number runs from zero to dofs_per_cell, even in the case of an
   * FEFaceValues or FESubfaceValues object.
   *
   * @param q_point Number of the quadrature point at which function is to be
   * evaluated
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  const double &
  shape_value(const unsigned int i, const unsigned int q_point) const;

  /**
   * Compute one vector component of the value of a shape function at a
   * quadrature point. If the finite element is scalar, then only component
   * zero is allowed and the return value equals that of the shape_value()
   * function. If the finite element is vector valued but all shape functions
   * are primitive (i.e. they are non-zero in only one component), then the
   * value returned by shape_value() equals that of this function for exactly
   * one component. This function is therefore only of greater interest if the
   * shape function is not primitive, but then it is necessary since the other
   * function cannot be used.
   *
   * @param i Number of the shape function $\varphi_i$ to be evaluated.
   *
   * @param q_point Number of the quadrature point at which function is to be
   * evaluated.
   *
   * @param component vector component to be evaluated.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  double
  shape_value_component(const unsigned int i,
                        const unsigned int q_point,
                        const unsigned int component) const;

  /**
   * Compute the gradient of the <tt>i</tt>th shape function at the
   * <tt>quadrature_point</tt>th quadrature point with respect to real cell
   * coordinates.  If you want to get the derivative in one of the coordinate
   * directions, use the appropriate function of the Tensor class to extract
   * one component of the Tensor returned by this function. Since only a
   * reference to the gradient's value is returned, there should be no major
   * performance drawback.
   *
   * If the shape function is vector-valued, then this returns the only
   * non-zero component. If the shape function has more than one non-zero
   * component (i.e. it is not primitive), then it will throw an exception of
   * type ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_grad_component() function.
   *
   * The same holds for the arguments of this function as for the
   * shape_value() function.
   *
   * @param i Number of the shape function $\varphi_i$ to be evaluated.
   *
   * @param q_point Number of the quadrature point at which function
   * is to be evaluated.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  const Tensor<1, spacedim> &
  shape_grad(const unsigned int i, const unsigned int q_point) const;

  /**
   * Return one vector component of the gradient of a shape function at a
   * quadrature point. If the finite element is scalar, then only component
   * zero is allowed and the return value equals that of the shape_grad()
   * function. If the finite element is vector valued but all shape functions
   * are primitive (i.e. they are non-zero in only one component), then the
   * value returned by shape_grad() equals that of this function for exactly
   * one component. This function is therefore only of greater interest if the
   * shape function is not primitive, but then it is necessary since the other
   * function cannot be used.
   *
   * The same holds for the arguments of this function as for the
   * shape_value_component() function.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  Tensor<1, spacedim>
  shape_grad_component(const unsigned int i,
                       const unsigned int q_point,
                       const unsigned int component) const;

  /**
   * Second derivatives of the <tt>i</tt>th shape function at the
   * <tt>q_point</tt>th quadrature point with respect to real cell
   * coordinates. If you want to get the derivatives in one of the coordinate
   * directions, use the appropriate function of the Tensor class to extract
   * one component. Since only a reference to the hessian values is returned,
   * there should be no major performance drawback.
   *
   * If the shape function is vector-valued, then this returns the only
   * non-zero component. If the shape function has more than one non-zero
   * component (i.e. it is not primitive), then throw an exception of type
   * ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_hessian_component() function.
   *
   * The same holds for the arguments of this function as for the
   * shape_value() function.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  const Tensor<2, spacedim> &
  shape_hessian(const unsigned int i, const unsigned int q_point) const;

  /**
   * Return one vector component of the hessian of a shape function at a
   * quadrature point. If the finite element is scalar, then only component
   * zero is allowed and the return value equals that of the shape_hessian()
   * function. If the finite element is vector valued but all shape functions
   * are primitive (i.e. they are non-zero in only one component), then the
   * value returned by shape_hessian() equals that of this function for
   * exactly one component. This function is therefore only of greater
   * interest if the shape function is not primitive, but then it is necessary
   * since the other function cannot be used.
   *
   * The same holds for the arguments of this function as for the
   * shape_value_component() function.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  Tensor<2, spacedim>
  shape_hessian_component(const unsigned int i,
                          const unsigned int q_point,
                          const unsigned int component) const;

  /**
   * Third derivatives of the <tt>i</tt>th shape function at the
   * <tt>q_point</tt>th quadrature point with respect to real cell
   * coordinates. If you want to get the 3rd derivatives in one of the
   * coordinate directions, use the appropriate function of the Tensor class
   * to extract one component. Since only a reference to the 3rd derivative
   * values is returned, there should be no major performance drawback.
   *
   * If the shape function is vector-valued, then this returns the only
   * non-zero component. If the shape function has more than one non-zero
   * component (i.e. it is not primitive), then throw an exception of type
   * ExcShapeFunctionNotPrimitive. In that case, use the
   * shape_3rdderivative_component() function.
   *
   * The same holds for the arguments of this function as for the
   * shape_value() function.
   *
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   */
  const Tensor<3, spacedim> &
  shape_3rd_derivative(const unsigned int i, const unsigned int q_point) const;

  /**
   * Return one vector component of the third derivative of a shape function
   * at a quadrature point. If the finite element is scalar, then only
   * component zero is allowed and the return value equals that of the
   * shape_3rdderivative() function. If the finite element is vector valued
   * but all shape functions are primitive (i.e. they are non-zero in only one
   * component), then the value returned by shape_3rdderivative() equals that
   * of this function for exactly one component. This function is therefore
   * only of greater interest if the shape function is not primitive, but then
   * it is necessary since the other function cannot be used.
   *
   * The same holds for the arguments of this function as for the
   * shape_value_component() function.
   *
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   */
  Tensor<3, spacedim>
  shape_3rd_derivative_component(const unsigned int i,
                                 const unsigned int q_point,
                                 const unsigned int component) const;

  /** @} */
  /// @name Access to values of global finite element fields
  /** @{ */

  /**
   * Return the values of a finite element function at the quadrature points
   * of the current cell, face, or subface (selected the last time the reinit()
   * function was called). That is, if the first argument @p fe_function is a
   * vector of nodal values of a finite element function $u_h(\mathbf x)$
   * defined on a DoFHandler object, then the output vector (the second
   * argument,
   * @p values) is the vector of values $u_h(\mathbf x_q^K)$ where $x_q^K$ are
   * the quadrature points on the current cell $K$.
   * This function is first discussed in the Results
   * section of step-4, and the related get_function_gradients() function
   * is also used in step-15 along with numerous other
   * tutorial programs.
   *
   * If the current cell is not active (i.e., it has children), then the finite
   * element function is, strictly speaking, defined by shape functions
   * that live on these child cells. Rather than evaluating the shape functions
   * on the child cells, with the quadrature points defined on the current
   * cell, this function first interpolates the finite element function to shape
   * functions defined on the current cell, and then evaluates this interpolated
   * function.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component.  To get values of multi-component
   * elements, there is another get_function_values() below,
   * returning a vector of vectors of results.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] values The values of the function specified by fe_function at
   * the quadrature points of the current cell.  The object is assume to
   * already have the correct size. The data type stored by this output vector
   * must be what you get when you multiply the values of shape function times
   * the type used to store the values of the unknowns $U_j$ of your finite
   * element vector $U$ (represented by the @p fe_function argument). This
   * happens to be equal to the type of the elements of the solution vector.
   *
   * @post <code>values[q]</code> will contain the value of the field
   * described by fe_function at the $q$th quadrature point.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <typename Number>
  void
  get_function_values(const ReadVector<Number> &fe_function,
                      std::vector<Number>      &values) const;

  /**
   * This function does the same as the other get_function_values(), but
   * applied to multi-component (vector-valued) elements. The meaning of the
   * arguments is as explained there.
   *
   * @post <code>values[q]</code> is a vector of values of the field described
   * by fe_function at the $q$th quadrature point. The size of the vector
   * accessed by <code>values[q]</code> equals the number of components of the
   * finite element, i.e. <code>values[q](c)</code> returns the value of the
   * $c$th vector component at the $q$th quadrature point.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <typename Number>
  void
  get_function_values(const ReadVector<Number>    &fe_function,
                      std::vector<Vector<Number>> &values) const;

  /**
   * Generate function values from an arbitrary vector. This function
   * does in essence the same as the first function of this name above,
   * except that it does not make the assumption that the input vector
   * corresponds to a DoFHandler that describes the unknowns of a finite
   * element field (and for which we would then assume that
   * `fe_function.size() == dof_handler.n_dofs()`). Rather, the nodal
   * values corresponding to the current cell are elements of an otherwise
   * arbitrary vector, and these elements are indexed by the second
   * argument to this function. What the rest of the `fe_function` input
   * argument corresponds to is of no consequence to this function.
   *
   * Given this, the function above corresponds to passing `fe_function`
   * as first argument to the current function, and using the
   * `local_dof_indices` array that results from the following call as
   * second argument to the current function:
   * @code
   *   cell->get_dof_indices (local_dof_indices);
   * @endcode
   * (See DoFCellAccessor::get_dof_indices() for more information.)
   *
   * Likewise, the function above is equivalent to calling
   * @code
   *   cell->get_dof_values (fe_function, local_dof_values);
   * @endcode
   * and then calling the current function with `local_dof_values` as
   * first argument, and an array with indices `{0,...,fe.dofs_per_cell-1}`
   * as second argument.
   *
   * The point of the current function is that one sometimes wants to
   * evaluate finite element functions at quadrature points with nodal
   * values that are not stored in a global vector -- for example, one could
   * modify these local values first, such as by applying a limiter
   * or by ensuring that all nodal values are positive, before evaluating
   * the finite element field that corresponds to these local values on the
   * current cell. Another application is where one wants to postprocess
   * the solution on a cell into a different finite element space on every
   * cell, without actually creating a corresponding DoFHandler -- in that
   * case, all one would compute is a local representation of that
   * postprocessed function, characterized by its nodal values; this function
   * then allows the evaluation of that representation at quadrature points.
   *
   * @param[in] fe_function A vector of nodal values. This vector can have
   *   an arbitrary size, as long as all elements index by `indices` can
   *   actually be accessed.
   *
   * @param[in] indices A vector of indices into `fe_function`. This vector
   *   must have length equal to the number of degrees of freedom on the
   *   current cell, and must identify elements in `fe_function` in the
   *   order in which degrees of freedom are indexed on the reference cell.
   *
   * @param[out] values A vector of values of the given finite element field,
   *   at the quadrature points on the current object.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <typename Number>
  void
  get_function_values(const ReadVector<Number> &fe_function,
                      const ArrayView<const types::global_dof_index> &indices,
                      std::vector<Number> &values) const;

  /**
   * Generate vector function values from an arbitrary vector.
   *
   * This function corresponds to the previous one, just for the vector-valued
   * case.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <typename Number>
  void
  get_function_values(const ReadVector<Number> &fe_function,
                      const ArrayView<const types::global_dof_index> &indices,
                      std::vector<Vector<Number>> &values) const;


  /**
   * Generate vector function values from an arbitrary vector. This
   * function is similar to the previous one, but the `indices`
   * vector may also be a multiple of the number of dofs per
   * cell. Then, the vectors in <tt>value</tt> should allow for the same
   * multiple of the components of the finite element.
   *
   * Depending on the value of the last argument, the outer vector of
   * <tt>values</tt> has either the length of the quadrature rule
   * (<tt>quadrature_points_fastest == false</tt>) or the length of components
   * to be filled <tt>quadrature_points_fastest == true</tt>. If <tt>p</tt> is
   * the current quadrature point number and <tt>i</tt> is the vector
   * component of the solution desired, the access to <tt>values</tt> is
   * <tt>values[p][i]</tt> if <tt>quadrature_points_fastest == false</tt>, and
   * <tt>values[i][p]</tt> otherwise.
   *
   * Since this function allows for fairly general combinations of argument
   * sizes, be aware that the checks on the arguments may not detect errors.
   *
   * @dealiiRequiresUpdateFlags{update_values}
   */
  template <typename Number>
  void
  get_function_values(const ReadVector<Number> &fe_function,
                      const ArrayView<const types::global_dof_index> &indices,
                      ArrayView<std::vector<Number>>                  values,
                      const bool quadrature_points_fastest) const;

  /** @} */
  /// @name Access to derivatives of global finite element fields
  /** @{ */

  /**
   * Return the gradients of a finite element function at the quadrature points
   * of the current cell, face, or subface (selected the last time the reinit()
   * function was called). That is, if the first argument @p fe_function is a
   * vector of nodal values of a finite element function $u_h(\mathbf x)$
   * defined on a DoFHandler object, then the output vector (the second
   * argument,
   * @p values) is the vector of values $\nabla u_h(\mathbf x_q^K)$ where
   * $x_q^K$ are the quadrature points on the current cell $K$. This function is
   * first discussed in the Results section of step-4, and it is also used in
   * step-15 along with numerous other tutorial programs.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component. There is a corresponding
   * function of the same name for vector-valued finite elements.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] gradients The gradients of the function specified by
   * fe_function at the quadrature points of the current cell.  The gradients
   * are computed in real space (as opposed to on the unit cell).  The object
   * is assume to already have the correct size. The data type stored by this
   * output vector must be what you get when you multiply the gradients of
   * shape function times the type used to store the values of the unknowns
   * $U_j$ of your finite element vector $U$ (represented by the @p
   * fe_function argument).
   *
   * @post <code>gradients[q]</code> will contain the gradient of the field
   * described by fe_function at the $q$th quadrature point.
   * <code>gradients[q][d]</code> represents the derivative in coordinate
   * direction $d$ at quadrature point $q$.
   *
   * @note The actual data type of the input vector may be either a
   * Vector&lt;T&gt;, BlockVector&lt;T&gt;, or one of the PETSc or Trilinos
   * vector wrapper classes. It represents a global vector of DoF values
   * associated with the DoFHandler object with which this FEValues object was
   * last initialized.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <typename Number>
  void
  get_function_gradients(
    const ReadVector<Number>                 &fe_function,
    std::vector<Tensor<1, spacedim, Number>> &gradients) const;

  /**
   * This function does the same as the other get_function_gradients(), but
   * applied to multi-component (vector-valued) elements. The meaning of the
   * arguments is as explained there.
   *
   * @post <code>gradients[q]</code> is a vector of gradients of the field
   * described by fe_function at the $q$th quadrature point. The size of the
   * vector accessed by <code>gradients[q]</code> equals the number of
   * components of the finite element, i.e. <code>gradients[q][c]</code>
   * returns the gradient of the $c$th vector component at the $q$th
   * quadrature point. Consequently, <code>gradients[q][c][d]</code> is the
   * derivative in coordinate direction $d$ of the $c$th vector component of
   * the vector field at quadrature point $q$ of the current cell.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <typename Number>
  void
  get_function_gradients(
    const ReadVector<Number>                              &fe_function,
    std::vector<std::vector<Tensor<1, spacedim, Number>>> &gradients) const;

  /**
   * This function relates to the first of the get_function_gradients() function
   * above in the same way as the get_function_values() with similar arguments
   * relates to the first of the get_function_values() functions. See there for
   * more information.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <typename Number>
  void
  get_function_gradients(
    const ReadVector<Number>                       &fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<Tensor<1, spacedim, Number>>       &gradients) const;

  /**
   * This function relates to the first of the get_function_gradients() function
   * above in the same way as the get_function_values() with similar arguments
   * relates to the first of the get_function_values() functions. See there for
   * more information.
   *
   * @dealiiRequiresUpdateFlags{update_gradients}
   */
  template <typename Number>
  void
  get_function_gradients(
    const ReadVector<Number>                           &fe_function,
    const ArrayView<const types::global_dof_index>     &indices,
    ArrayView<std::vector<Tensor<1, spacedim, Number>>> gradients,
    const bool quadrature_points_fastest = false) const;

  /** @} */
  /// @name Access to second derivatives
  ///
  /// Hessian matrices and Laplacians of global finite element fields
  /** @{ */

  /**
   * Compute the tensor of second derivatives of a finite element at the
   * quadrature points of a cell. This function is the equivalent of the
   * corresponding get_function_values() function (see there for more
   * information) but evaluates the finite element field's second derivatives
   * instead of its value.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component. There is a corresponding
   * function of the same name for vector-valued finite elements.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] hessians The Hessians of the function specified by
   * fe_function at the quadrature points of the current cell.  The Hessians
   * are computed in real space (as opposed to on the unit cell).  The object
   * is assume to already have the correct size. The data type stored by this
   * output vector must be what you get when you multiply the Hessians of
   * shape function times the type used to store the values of the unknowns
   * $U_j$ of your finite element vector $U$ (represented by the @p
   * fe_function argument).
   *
   * @post <code>hessians[q]</code> will contain the Hessian of the field
   * described by fe_function at the $q$th quadrature point.
   * <code>hessians[q][i][j]</code> represents the $(i,j)$th component of the
   * matrix of second derivatives at quadrature point $q$.
   *
   * @note The actual data type of the input vector may be either a
   * Vector&lt;T&gt;, BlockVector&lt;T&gt;, or one of the PETSc or Trilinos
   * vector wrapper classes. It represents a global vector of DoF values
   * associated with the DoFHandler object with which this FEValues object was
   * last initialized.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_hessians(
    const ReadVector<Number>                 &fe_function,
    std::vector<Tensor<2, spacedim, Number>> &hessians) const;

  /**
   * This function does the same as the other get_function_hessians(), but
   * applied to multi-component (vector-valued) elements. The meaning of the
   * arguments is as explained there.
   *
   * @post <code>hessians[q]</code> is a vector of Hessians of the field
   * described by fe_function at the $q$th quadrature point. The size of the
   * vector accessed by <code>hessians[q]</code> equals the number of
   * components of the finite element, i.e. <code>hessians[q][c]</code>
   * returns the Hessian of the $c$th vector component at the $q$th quadrature
   * point. Consequently, <code>hessians[q][c][i][j]</code> is the $(i,j)$th
   * component of the matrix of second derivatives of the $c$th vector
   * component of the vector field at quadrature point $q$ of the current
   * cell.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_hessians(
    const ReadVector<Number>                              &fe_function,
    std::vector<std::vector<Tensor<2, spacedim, Number>>> &hessians,
    const bool quadrature_points_fastest = false) const;

  /**
   * This function relates to the first of the get_function_hessians() function
   * above in the same way as the get_function_values() with similar arguments
   * relates to the first of the get_function_values() functions. See there for
   * more information.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_hessians(
    const ReadVector<Number>                       &fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<Tensor<2, spacedim, Number>>       &hessians) const;

  /**
   * This function relates to the first of the get_function_hessians() function
   * above in the same way as the get_function_values() with similar arguments
   * relates to the first of the get_function_values() functions. See there for
   * more information.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_hessians(
    const ReadVector<Number>                           &fe_function,
    const ArrayView<const types::global_dof_index>     &indices,
    ArrayView<std::vector<Tensor<2, spacedim, Number>>> hessians,
    const bool quadrature_points_fastest = false) const;

  /**
   * Compute the (scalar) Laplacian (i.e. the trace of the tensor of second
   * derivatives) of a finite element at the quadrature points of a cell. This
   * function is the equivalent of the corresponding get_function_values()
   * function (see there for more information) but evaluates the finite
   * element field's second derivatives instead of its value.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component. There is a corresponding
   * function of the same name for vector-valued finite elements.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] laplacians The Laplacians of the function specified by
   * fe_function at the quadrature points of the current cell.  The Laplacians
   * are computed in real space (as opposed to on the unit cell).  The object
   * is assume to already have the correct size. The data type stored by this
   * output vector must be what you get when you multiply the Laplacians of
   * shape function times the type used to store the values of the unknowns
   * $U_j$ of your finite element vector $U$ (represented by the @p
   * fe_function argument). This happens to be equal to the type of the
   * elements of the input vector.
   *
   * @post <code>laplacians[q]</code> will contain the Laplacian of the field
   * described by fe_function at the $q$th quadrature point.
   *
   * @post For each component of the output vector, there holds
   * <code>laplacians[q]=trace(hessians[q])</code>, where <tt>hessians</tt>
   * would be the output of the get_function_hessians() function.
   *
   * @note The actual data type of the input vector may be either a
   * Vector&lt;T&gt;, BlockVector&lt;T&gt;, or one of the PETSc or Trilinos
   * vector wrapper classes. It represents a global vector of DoF values
   * associated with the DoFHandler object with which this FEValues object was
   * last initialized.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_laplacians(const ReadVector<Number> &fe_function,
                          std::vector<Number>      &laplacians) const;

  /**
   * This function does the same as the other get_function_laplacians(), but
   * applied to multi-component (vector-valued) elements. The meaning of the
   * arguments is as explained there.
   *
   * @post <code>laplacians[q]</code> is a vector of Laplacians of the field
   * described by fe_function at the $q$th quadrature point. The size of the
   * vector accessed by <code>laplacians[q]</code> equals the number of
   * components of the finite element, i.e. <code>laplacians[q][c]</code>
   * returns the Laplacian of the $c$th vector component at the $q$th
   * quadrature point.
   *
   * @post For each component of the output vector, there holds
   * <code>laplacians[q][c]=trace(hessians[q][c])</code>, where
   * <tt>hessians</tt> would be the output of the get_function_hessians()
   * function.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_laplacians(const ReadVector<Number>    &fe_function,
                          std::vector<Vector<Number>> &laplacians) const;

  /**
   * This function relates to the first of the get_function_laplacians()
   * function above in the same way as the get_function_values() with similar
   * arguments relates to the first of the get_function_values() functions. See
   * there for more information.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_laplacians(
    const ReadVector<Number>                       &fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<Number>                            &laplacians) const;

  /**
   * This function relates to the first of the get_function_laplacians()
   * function above in the same way as the get_function_values() with similar
   * arguments relates to the first of the get_function_values() functions. See
   * there for more information.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_laplacians(
    const ReadVector<Number>                       &fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<Vector<Number>>                    &laplacians) const;

  /**
   * This function relates to the first of the get_function_laplacians()
   * function above in the same way as the get_function_values() with similar
   * arguments relates to the first of the get_function_values() functions. See
   * there for more information.
   *
   * @dealiiRequiresUpdateFlags{update_hessians}
   */
  template <typename Number>
  void
  get_function_laplacians(
    const ReadVector<Number>                       &fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<std::vector<Number>>               &laplacians,
    const bool quadrature_points_fastest = false) const;

  /** @} */
  /// @name Access to third derivatives of global finite element fields
  /** @{ */

  /**
   * Compute the tensor of third derivatives of a finite element at the
   * quadrature points of a cell. This function is the equivalent of the
   * corresponding get_function_values() function (see there for more
   * information) but evaluates the finite element field's third derivatives
   * instead of its value.
   *
   * This function may only be used if the finite element in use is a scalar
   * one, i.e. has only one vector component. There is a corresponding
   * function of the same name for vector-valued finite elements.
   *
   * @param[in] fe_function A vector of values that describes (globally) the
   * finite element function that this function should evaluate at the
   * quadrature points of the current cell.
   *
   * @param[out] third_derivatives The third derivatives of the function
   * specified by fe_function at the quadrature points of the current cell.
   * The third derivatives are computed in real space (as opposed to on the
   * unit cell).  The object is assumed to already have the correct size. The
   * data type stored by this output vector must be what you get when you
   * multiply the third derivatives of shape function times the type used to
   * store the values of the unknowns $U_j$ of your finite element vector $U$
   * (represented by the @p fe_function argument).
   *
   * @post <code>third_derivatives[q]</code> will contain the third
   * derivatives of the field described by fe_function at the $q$th quadrature
   * point. <code>third_derivatives[q][i][j][k]</code> represents the
   * $(i,j,k)$th component of the 3rd order tensor of third derivatives at
   * quadrature point $q$.
   *
   * @note The actual data type of the input vector may be either a
   * Vector&lt;T&gt;, BlockVector&lt;T&gt;, or one of the PETSc or Trilinos
   * vector wrapper classes. It represents a global vector of DoF values
   * associated with the DoFHandler object with which this FEValues object was
   * last initialized.
   *
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   */
  template <typename Number>
  void
  get_function_third_derivatives(
    const ReadVector<Number>                 &fe_function,
    std::vector<Tensor<3, spacedim, Number>> &third_derivatives) const;

  /**
   * This function does the same as the other
   * get_function_third_derivatives(), but applied to multi-component
   * (vector-valued) elements. The meaning of the arguments is as explained
   * there.
   *
   * @post <code>third_derivatives[q]</code> is a vector of third derivatives
   * of the field described by fe_function at the $q$th quadrature point. The
   * size of the vector accessed by <code>third_derivatives[q]</code> equals
   * the number of components of the finite element, i.e.
   * <code>third_derivatives[q][c]</code> returns the third derivative of the
   * $c$th vector component at the $q$th quadrature point. Consequently,
   * <code>third_derivatives[q][c][i][j][k]</code> is the $(i,j,k)$th
   * component of the tensor of third derivatives of the $c$th vector
   * component of the vector field at quadrature point $q$ of the current
   * cell.
   *
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   */
  template <typename Number>
  void
  get_function_third_derivatives(
    const ReadVector<Number>                              &fe_function,
    std::vector<std::vector<Tensor<3, spacedim, Number>>> &third_derivatives,
    const bool quadrature_points_fastest = false) const;

  /**
   * This function relates to the first of the get_function_third_derivatives()
   * function above in the same way as the get_function_values() with similar
   * arguments relates to the first of the get_function_values() functions. See
   * there for more information.
   *
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   */
  template <typename Number>
  void
  get_function_third_derivatives(
    const ReadVector<Number>                       &fe_function,
    const ArrayView<const types::global_dof_index> &indices,
    std::vector<Tensor<3, spacedim, Number>>       &third_derivatives) const;

  /**
   * This function relates to the first of the get_function_third_derivatives()
   * function above in the same way as the get_function_values() with similar
   * arguments relates to the first of the get_function_values() functions. See
   * there for more information.
   *
   * @dealiiRequiresUpdateFlags{update_3rd_derivatives}
   */
  template <typename Number>
  void
  get_function_third_derivatives(
    const ReadVector<Number>                           &fe_function,
    const ArrayView<const types::global_dof_index>     &indices,
    ArrayView<std::vector<Tensor<3, spacedim, Number>>> third_derivatives,
    const bool quadrature_points_fastest = false) const;
  /** @} */

  /// @name Cell degrees of freedom
  /** @{ */

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero (inclusive) to `dofs_per_cell` (exclusive). This allows
   * one to write code using range-based `for` loops of the following kind:
   * @code
   *   FEValues<dim>      fe_values (...);
   *   FullMatrix<double> cell_matrix (...);
   *
   *   for (auto &cell : dof_handler.active_cell_iterators())
   *     {
   *       cell_matrix = 0;
   *       fe_values.reinit(cell);
   *       for (const auto q : fe_values.quadrature_point_indices())
   *         for (const auto i : fe_values.dof_indices())
   *           for (const auto j : fe_values.dof_indices())
   *             cell_matrix(i,j) += ...; // Do something for DoF indices (i,j)
   *                                      // at quadrature point q
   *     }
   * @endcode
   * Here, we are looping over all degrees of freedom on all cells, with
   * `i` and `j` taking on all valid indices for cell degrees of freedom, as
   * defined by the finite element passed to `fe_values`.
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  dof_indices() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from @p start_dof_index (inclusive) to `dofs_per_cell` (exclusive).
   * This allows one to write code using range-based `for` loops of the
   * following kind:
   * @code
   *   FEValues<dim>      fe_values (...);
   *   FullMatrix<double> cell_matrix (...);
   *
   *   for (auto &cell : dof_handler.active_cell_iterators())
   *     {
   *       cell_matrix = 0;
   *       fe_values.reinit(cell);
   *       for (const auto q : fe_values.quadrature_point_indices())
   *         for (const auto i : fe_values.dof_indices())
   *           for (const auto j : fe_values.dof_indices_starting_at(i))
   *             cell_matrix(i,j) += ...; // Do something for DoF indices (i,j)
   *                                      // at quadrature point q
   *     }
   * @endcode
   * Here, we are looping over all local degrees of freedom on all cells, with
   * `i` taking on all valid indices for cell degrees of freedom, as
   * defined by the finite element passed to `fe_values`, and `j` taking
   * on a specified subset of `i`'s range, starting at `i` itself and ending at
   * the number of cell degrees of freedom. In this way, we can construct the
   * upper half and the diagonal of a @ref GlossStiffnessMatrix "stiffness matrix" contribution (assuming it
   * is symmetric, and that only one half of it needs to be computed), for
   * example.
   *
   * @note If the @p start_dof_index is equal to the number of DoFs in the cell,
   * then the returned index range is empty.
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  dof_indices_starting_at(const unsigned int start_dof_index) const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero (inclusive) to @p end_dof_index (inclusive). This allows
   * one to write code using range-based `for` loops of the following kind:
   * @code
   *   FEValues<dim>      fe_values (...);
   *   FullMatrix<double> cell_matrix (...);
   *
   *   for (auto &cell : dof_handler.active_cell_iterators())
   *     {
   *       cell_matrix = 0;
   *       fe_values.reinit(cell);
   *       for (const auto q : fe_values.quadrature_point_indices())
   *         for (const auto i : fe_values.dof_indices())
   *           for (const auto j : fe_values.dof_indices_ending_at(i))
   *             cell_matrix(i,j) += ...; // Do something for DoF indices (i,j)
   *                                      // at quadrature point q
   *     }
   * @endcode
   * Here, we are looping over all local degrees of freedom on all cells, with
   * `i` taking on all valid indices for cell degrees of freedom, as
   * defined by the finite element passed to `fe_values`, and `j` taking
   * on a specified subset of `i`'s range, starting at zero and ending at
   * `i` itself. In this way, we can construct the lower half and the
   * diagonal of a @ref GlossStiffnessMatrix "stiffness matrix" contribution (assuming it is symmetric, and
   * that only one half of it needs to be computed), for example.
   *
   * @note If the @p end_dof_index is equal to zero, then the returned index
   * range is empty.
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  dof_indices_ending_at(const unsigned int end_dof_index) const;

  /** @} */

  /// @name Geometry of the cell
  /** @{ */

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero to `n_quadrature_points`. This allows to write code
   * using range-based `for` loops of the following kind:
   * @code
   *   FEValues<dim> fe_values (...);
   *
   *   for (auto &cell : dof_handler.active_cell_iterators())
   *     {
   *       fe_values.reinit(cell);
   *       for (const auto q_point : fe_values.quadrature_point_indices())
   *         ... do something at the quadrature point ...
   *     }
   * @endcode
   * Here, we are looping over all quadrature points on all cells, with
   * `q_point` taking on all valid indices for quadrature points, as defined
   * by the quadrature rule passed to `fe_values`.
   *
   * @see CPP11
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  quadrature_point_indices() const;

  /**
   * Return the location of the <tt>q_point</tt>th quadrature point in
   * real space.
   *
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   */
  const Point<spacedim> &
  quadrature_point(const unsigned int q_point) const;

  /**
   * Return a reference to the vector of quadrature points in real space.
   *
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   */
  const std::vector<Point<spacedim>> &
  get_quadrature_points() const;

  /**
   * Mapped quadrature weight. If this object refers to a volume evaluation
   * (i.e. the derived class is of type FEValues), then this is the Jacobi
   * determinant times the weight of the <tt>q_point</tt>th unit quadrature
   * point.
   *
   * For surface evaluations (i.e. classes FEFaceValues or FESubfaceValues),
   * it is the mapped surface element times the weight of the quadrature
   * point.
   *
   * You can think of the quantity returned by this function as the volume or
   * surface element $dx, ds$ in the integral that we implement here by
   * quadrature.
   *
   * @dealiiRequiresUpdateFlags{update_JxW_values}
   */
  double
  JxW(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by JxW().
   */
  const std::vector<double> &
  get_JxW_values() const;

  /**
   * Return the Jacobian of the transformation at the specified quadrature
   * point, i.e.  $J_{ij}=dx_i/d\hat x_j$
   *
   * @dealiiRequiresUpdateFlags{update_jacobians}
   */
  const DerivativeForm<1, dim, spacedim> &
  jacobian(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by
   * jacobian().
   *
   * @dealiiRequiresUpdateFlags{update_jacobians}
   */
  const std::vector<DerivativeForm<1, dim, spacedim>> &
  get_jacobians() const;

  /**
   * Return the second derivative of the transformation from unit to real
   * cell, i.e. the first derivative of the Jacobian, at the specified
   * quadrature point, i.e. $G_{ijk}=dJ_{jk}/d\hat x_i$.
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_grads}
   */
  const DerivativeForm<2, dim, spacedim> &
  jacobian_grad(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by
   * jacobian_grads().
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_grads}
   */
  const std::vector<DerivativeForm<2, dim, spacedim>> &
  get_jacobian_grads() const;

  /**
   * Return the second derivative of the transformation from unit to real
   * cell, i.e. the first derivative of the Jacobian, at the specified
   * quadrature point, pushed forward to the real cell coordinates, i.e.
   * $G_{ijk}=dJ_{iJ}/d\hat x_K (J_{jJ})^{-1} (J_{kK})^{-1}$.
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_grads}
   */
  const Tensor<3, spacedim> &
  jacobian_pushed_forward_grad(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by
   * jacobian_pushed_forward_grads().
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_grads}
   */
  const std::vector<Tensor<3, spacedim>> &
  get_jacobian_pushed_forward_grads() const;

  /**
   * Return the third derivative of the transformation from unit to real cell,
   * i.e. the second derivative of the Jacobian, at the specified quadrature
   * point, i.e. $G_{ijkl}=\frac{d^2J_{ij}}{d\hat x_k d\hat x_l}$.
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_2nd_derivatives}
   */
  const DerivativeForm<3, dim, spacedim> &
  jacobian_2nd_derivative(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by
   * jacobian_2nd_derivatives().
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_2nd_derivatives}
   */
  const std::vector<DerivativeForm<3, dim, spacedim>> &
  get_jacobian_2nd_derivatives() const;

  /**
   * Return the third derivative of the transformation from unit to real cell,
   * i.e. the second derivative of the Jacobian, at the specified quadrature
   * point, pushed forward to the real cell coordinates, i.e.
   * $G_{ijkl}=\frac{d^2J_{iJ}}{d\hat x_K d\hat x_L} (J_{jJ})^{-1}
   * (J_{kK})^{-1}(J_{lL})^{-1}$.
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_2nd_derivatives}
   */
  const Tensor<4, spacedim> &
  jacobian_pushed_forward_2nd_derivative(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by
   * jacobian_pushed_forward_2nd_derivatives().
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_2nd_derivatives}
   */
  const std::vector<Tensor<4, spacedim>> &
  get_jacobian_pushed_forward_2nd_derivatives() const;

  /**
   * Return the fourth derivative of the transformation from unit to real
   * cell, i.e. the third derivative of the Jacobian, at the specified
   * quadrature point, i.e. $G_{ijklm}=\frac{d^2J_{ij}}{d\hat x_k d\hat x_l
   * d\hat x_m}$.
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_3rd_derivatives}
   */
  const DerivativeForm<4, dim, spacedim> &
  jacobian_3rd_derivative(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by
   * jacobian_3rd_derivatives().
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_3rd_derivatives}
   */
  const std::vector<DerivativeForm<4, dim, spacedim>> &
  get_jacobian_3rd_derivatives() const;

  /**
   * Return the fourth derivative of the transformation from unit to real
   * cell, i.e. the third derivative of the Jacobian, at the specified
   * quadrature point, pushed forward to the real cell coordinates, i.e.
   * $G_{ijklm}=\frac{d^3J_{iJ}}{d\hat x_K d\hat x_L d\hat x_M} (J_{jJ})^{-1}
   * (J_{kK})^{-1} (J_{lL})^{-1} (J_{mM})^{-1}$.
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_3rd_derivatives}
   */
  const Tensor<5, spacedim> &
  jacobian_pushed_forward_3rd_derivative(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by
   * jacobian_pushed_forward_3rd_derivatives().
   *
   * @dealiiRequiresUpdateFlags{update_jacobian_pushed_forward_2nd_derivatives}
   */
  const std::vector<Tensor<5, spacedim>> &
  get_jacobian_pushed_forward_3rd_derivatives() const;

  /**
   * Return the inverse Jacobian of the transformation at the specified
   * quadrature point, i.e.  $J_{ij}=d\hat x_i/dx_j$
   *
   * @dealiiRequiresUpdateFlags{update_inverse_jacobians}
   */
  const DerivativeForm<1, spacedim, dim> &
  inverse_jacobian(const unsigned int q_point) const;

  /**
   * Return a reference to the array holding the values returned by
   * inverse_jacobian().
   *
   * @dealiiRequiresUpdateFlags{update_inverse_jacobians}
   */
  const std::vector<DerivativeForm<1, spacedim, dim>> &
  get_inverse_jacobians() const;

  /**
   * Return the normal vector at a quadrature point. If you call this
   * function for a face (i.e., when using a FEFaceValues or FESubfaceValues
   * object), then this function returns the outward normal vector to
   * the cell at the <tt>q_point</tt>th quadrature point of the face.
   *
   * In contrast, if you call this function for a cell of codimension one
   * (i.e., when using a `FEValues<dim,spacedim>` object with
   * `spacedim>dim`), then this function returns the normal vector to the
   * cell -- in other words, an approximation to the normal vector to the
   * manifold in which the triangulation is embedded. There are of
   * course two normal directions to a manifold in that case, and this
   * function returns the "up" direction as induced by the numbering of the
   * vertices.
   *
   * The length of the vector is normalized to one.
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  const Tensor<1, spacedim> &
  normal_vector(const unsigned int q_point) const;

  /**
   * Return the normal vectors at all quadrature points represented by
   * this object. See the normal_vector() function for what the normal
   * vectors represent.
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  const std::vector<Tensor<1, spacedim>> &
  get_normal_vectors() const;

  /** @} */

  /// @name Extractors Methods to extract individual components
  /** @{ */

  /**
   * Create a view of the current FEValues object that represents a particular
   * scalar component of the possibly vector-valued finite element. The
   * concept of views is explained in the documentation of the namespace
   * FEValuesViews and in particular in the
   * @ref vector_valued
   * topic.
   */
  const FEValuesViews::Scalar<dim, spacedim> &
  operator[](const FEValuesExtractors::Scalar &scalar) const;

  /**
   * Create a view of the current FEValues object that represents a set of
   * <code>dim</code> scalar components (i.e. a vector) of the vector-valued
   * finite element. The concept of views is explained in the documentation of
   * the namespace FEValuesViews and in particular in the
   * @ref vector_valued
   * topic.
   */
  const FEValuesViews::Vector<dim, spacedim> &
  operator[](const FEValuesExtractors::Vector &vector) const;

  /**
   * Create a view of the current FEValues object that represents a set of
   * <code>(dim*dim + dim)/2</code> scalar components (i.e. a symmetric 2nd
   * order tensor) of the vector-valued finite element. The concept of views
   * is explained in the documentation of the namespace FEValuesViews and in
   * particular in the
   * @ref vector_valued
   * topic.
   */
  const FEValuesViews::SymmetricTensor<2, dim, spacedim> &
  operator[](const FEValuesExtractors::SymmetricTensor<2> &tensor) const;


  /**
   * Create a view of the current FEValues object that represents a set of
   * <code>(dim*dim)</code> scalar components (i.e. a 2nd order tensor) of the
   * vector-valued finite element. The concept of views is explained in the
   * documentation of the namespace FEValuesViews and in particular in the
   * @ref vector_valued
   * topic.
   */
  const FEValuesViews::Tensor<2, dim, spacedim> &
  operator[](const FEValuesExtractors::Tensor<2> &tensor) const;

  /** @} */

  /// @name Access to the raw data
  /** @{ */

  /**
   * Constant reference to the selected mapping object.
   */
  const Mapping<dim, spacedim> &
  get_mapping() const;

  /**
   * Constant reference to the selected finite element object.
   */
  const FiniteElement<dim, spacedim> &
  get_fe() const;

  /**
   * Return the update flags set for this object.
   */
  UpdateFlags
  get_update_flags() const;

  /**
   * Return a triangulation iterator to the current cell.
   */
  typename Triangulation<dim, spacedim>::cell_iterator
  get_cell() const;

  /**
   * Return the relation of the current cell to the previous cell. This allows
   * re-use of some cell data (like local matrices for equations with constant
   * coefficients) if the result is <tt>CellSimilarity::translation</tt>.
   */
  CellSimilarity::Similarity
  get_cell_similarity() const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t
  memory_consumption() const;
  /** @} */


  /**
   * This exception is thrown if FEValuesBase is asked to return the value of
   * a field which was not required by the UpdateFlags for this FEValuesBase.
   *
   * @ingroup Exceptions
   */
  DeclException1(
    ExcAccessToUninitializedField,
    std::string,
    << "You are requesting information from an FEValues/FEFaceValues/FESubfaceValues "
    << "object for which this kind of information has not been computed. What "
    << "information these objects compute is determined by the update_* flags you "
    << "pass to the constructor. Here, the operation you are attempting requires "
    << "the <" << arg1
    << "> flag to be set, but it was apparently not specified "
    << "upon construction.");

  /**
   * FEValues::reinit() has not been called for any cell.
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcNotReinited,
                   "FEValues object is not reinit'ed to any cell");

  /**
   * Mismatch between the FEValues FiniteElement and
   * cell->get_dof_handler().get_fe()
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(
    ExcFEDontMatch,
    "The FiniteElement you provided to FEValues and the FiniteElement that belongs "
    "to the DoFHandler that provided the cell iterator do not match.");
  /**
   * A given shape function is not primitive, but it needs to be.
   *
   * @ingroup Exceptions
   */
  DeclException1(ExcShapeFunctionNotPrimitive,
                 int,
                 << "The shape function with index " << arg1
                 << " is not primitive, i.e. it is vector-valued and "
                 << "has more than one non-zero vector component. This "
                 << "function cannot be called for these shape functions. "
                 << "Maybe you want to use the same function with the "
                 << "_component suffix?");

  /**
   * The given FiniteElement is not a primitive element, see
   * FiniteElement::is_primitive().
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(
    ExcFENotPrimitive,
    "The given FiniteElement is not a primitive element but the requested operation "
    "only works for those. See FiniteElement::is_primitive() for more information.");

protected:
  /**
   * Objects of the FEValues class need to store an iterator
   * to the present cell in order to be able to extract the values of the
   * degrees of freedom on this cell in the get_function_values() and assorted
   * functions.
   *
   * The problem is that the iterators given to the various reinit() functions
   * can either be Triangulation iterators, or DoFHandler cell or level
   * iterators. All three are valid, and provide different functionality that is
   * used in different contexts; as a consequence we need to be able to store
   * all three. This class provides the ability to store an object of any of
   * these types, via a member variable that is a std::variant that encapsulates
   * an object of any of the three types. Because a std::variant always stores
   * an object of *one* of these types, we wrap the std::variant object into a
   * std::optional that allows us to encode a "not yet initialized" state.
   */
  class CellIteratorWrapper
  {
  public:
    DeclExceptionMsg(
      ExcNeedsDoFHandler,
      "You have previously called the FEValues::reinit() function with a "
      "cell iterator of type Triangulation<dim,spacedim>::cell_iterator. However, "
      "when you do this, you cannot call some functions in the FEValues "
      "class, such as the get_function_values/gradients/hessians/third_derivatives "
      "functions. If you need these functions, then you need to call "
      "FEValues::reinit() with an iterator type that allows to extract "
      "degrees of freedom, such as DoFHandler<dim,spacedim>::cell_iterator.");

    /**
     * Constructor. Creates an unusable object that is not associated with
     * any cell at all.
     */
    CellIteratorWrapper() = default;

    /**
     * Constructor.
     */
    CellIteratorWrapper(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell);

    /**
     * Constructor.
     */
    CellIteratorWrapper(
      const typename DoFHandler<dim, spacedim>::cell_iterator &cell);

    /**
     * Constructor.
     */
    CellIteratorWrapper(
      const typename DoFHandler<dim, spacedim>::level_cell_iterator &cell);

    /**
     * Indicate whether FEValues::reinit() was called.
     */
    bool
    is_initialized() const;

    /**
     * Conversion operator to an iterator for triangulations. This
     * conversion is implicit for the original iterators, since they are derived
     * classes. However, since here we have kind of a parallel class hierarchy,
     * we have to have a conversion operator.
     */
    operator typename Triangulation<dim, spacedim>::cell_iterator() const;

    /**
     * Return the number of degrees of freedom the DoF
     * handler object has to which the iterator belongs to.
     */
    types::global_dof_index
    n_dofs_for_dof_handler() const;

    /**
     * Call @p get_interpolated_dof_values of the iterator with the
     * given arguments.
     */
    template <typename Number>
    void
    get_interpolated_dof_values(const ReadVector<Number> &in,
                                Vector<Number>           &out) const;

  private:
    /**
     * The cell in question, if one has been assigned to this object. The
     * concrete data type can either be a Triangulation cell iterator, a
     * DoFHandler cell iterator, or a DoFHandler level cell iterator.
     */
    std::optional<
      std::variant<typename Triangulation<dim, spacedim>::cell_iterator,
                   typename DoFHandler<dim, spacedim>::cell_iterator,
                   typename DoFHandler<dim, spacedim>::level_cell_iterator>>
      cell;
  };

  /**
   * Store the cell selected last time the reinit() function was called.  This
   * is necessary for the <tt>get_function_*</tt> functions as well as the
   * functions of same name in the extractor classes.
   */
  CellIteratorWrapper present_cell;

  /**
   * A signal connection we use to ensure we get informed whenever the
   * triangulation changes by refinement. We need to know about that because
   * it invalidates all cell iterators and, as part of that, the
   * 'present_cell' iterator we keep around between subsequent calls to
   * reinit() in order to compute the cell similarity.
   */
  boost::signals2::connection tria_listener_refinement;

  /**
   * A signal connection we use to ensure we get informed whenever the
   * triangulation changes by mesh transformations. We need to know about that
   * because it invalidates all cell iterators and, as part of that, the
   * 'present_cell' iterator we keep around between subsequent calls to
   * reinit() in order to compute the cell similarity.
   */
  boost::signals2::connection tria_listener_mesh_transform;

  /**
   * A function that is connected to the triangulation in order to reset the
   * stored 'present_cell' iterator to an invalid one whenever the
   * triangulation is changed and the iterator consequently becomes invalid.
   */
  void
  invalidate_present_cell();

  /**
   * This function is called by the various reinit() functions in derived
   * classes. Given the cell indicated by the argument, test whether we have
   * to throw away the previously stored present_cell argument because it
   * would require us to compare cells from different triangulations. In
   * checking all this, also make sure that we have tria_listener connected to
   * the triangulation to which we will set present_cell right after calling
   * this function.
   */
  void
  maybe_invalidate_previous_present_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell);

  /**
   * A pointer to the mapping object associated with this FEValues object.
   */
  const ObserverPointer<const Mapping<dim, spacedim>,
                        FEValuesBase<dim, spacedim>>
    mapping;

  /**
   * A pointer to the internal data object of mapping, obtained from
   * Mapping::get_data(), Mapping::get_face_data(), or
   * Mapping::get_subface_data().
   */
  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
    mapping_data;

  /**
   * An object into which the Mapping::fill_fe_values() and similar functions
   * place their output.
   */
  internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    mapping_output;

  /**
   * A pointer to the finite element object associated with this FEValues
   * object.
   */
  const ObserverPointer<const FiniteElement<dim, spacedim>,
                        FEValuesBase<dim, spacedim>>
    fe;

  /**
   * A pointer to the internal data object of finite element, obtained from
   * FiniteElement::get_data(), Mapping::get_face_data(), or
   * FiniteElement::get_subface_data().
   */
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
    fe_data;

  /**
   * An object into which the FiniteElement::fill_fe_values() and similar
   * functions place their output.
   */
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    finite_element_output;


  /**
   * Original update flags handed to the constructor of FEValues.
   */
  UpdateFlags update_flags;

  /**
   * Initialize some update flags. Called from the @p initialize functions of
   * derived classes, which are in turn called from their constructors.
   *
   * Basically, this function finds out using the finite element and mapping
   * object already stored which flags need to be set to compute everything
   * the user wants, as expressed through the flags passed as argument.
   */
  UpdateFlags
  compute_update_flags(const UpdateFlags update_flags) const;

  /**
   * An enum variable that can store different states of the current cell in
   * comparison to the previously visited cell. If wanted, additional states
   * can be checked here and used in one of the methods used during reinit.
   */
  CellSimilarity::Similarity cell_similarity;

  /**
   * A function that checks whether the new cell is similar to the one
   * previously used. Then, a significant amount of the data can be reused,
   * e.g. the derivatives of the basis functions in real space, shape_grad.
   */
  void
  check_cell_similarity(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell);

private:
  /**
   * A cache for all possible FEValuesViews objects.
   */
  dealii::internal::FEValuesViews::Cache<dim, spacedim> fe_values_views_cache;

  /**
   * Whether checking for cell similarity is allowed.
   */
  bool check_for_cell_similarity_allowed;

  // Make the view classes friends of this class, since they access internal
  // data.
  template <int, int>
  friend class FEValuesViews::Scalar;
  template <int, int>
  friend class FEValuesViews::Vector;
  template <int, int, int>
  friend class FEValuesViews::SymmetricTensor;
  template <int, int, int>
  friend class FEValuesViews::Tensor;
};

#ifndef DOXYGEN

/*---------------------- Inline functions: FEValuesBase ---------------------*/

template <int dim, int spacedim>
inline const FEValuesViews::Scalar<dim, spacedim> &
FEValuesBase<dim, spacedim>::operator[](
  const FEValuesExtractors::Scalar &scalar) const
{
  AssertIndexRange(scalar.component, fe_values_views_cache.scalars.size());

  return fe_values_views_cache.scalars[scalar.component].value_or_initialize(
    [scalar, this]() {
      return FEValuesViews::Scalar<dim, spacedim>(*this, scalar.component);
    });
}



template <int dim, int spacedim>
inline const FEValuesViews::Vector<dim, spacedim> &
FEValuesBase<dim, spacedim>::operator[](
  const FEValuesExtractors::Vector &vector) const
{
  AssertIndexRange(vector.first_vector_component,
                   fe_values_views_cache.vectors.size());

  return fe_values_views_cache.vectors[vector.first_vector_component]
    .value_or_initialize([vector, this]() {
      return FEValuesViews::Vector<dim, spacedim>(
        *this, vector.first_vector_component);
    });
}



template <int dim, int spacedim>
inline const FEValuesViews::SymmetricTensor<2, dim, spacedim> &
FEValuesBase<dim, spacedim>::operator[](
  const FEValuesExtractors::SymmetricTensor<2> &tensor) const
{
  Assert(
    tensor.first_tensor_component <
      fe_values_views_cache.symmetric_second_order_tensors.size(),
    ExcIndexRange(tensor.first_tensor_component,
                  0,
                  fe_values_views_cache.symmetric_second_order_tensors.size()));

  return fe_values_views_cache
    .symmetric_second_order_tensors[tensor.first_tensor_component]
    .value_or_initialize([tensor, this]() {
      return FEValuesViews::SymmetricTensor<2, dim, spacedim>(
        *this, tensor.first_tensor_component);
    });
}



template <int dim, int spacedim>
inline const FEValuesViews::Tensor<2, dim, spacedim> &
FEValuesBase<dim, spacedim>::operator[](
  const FEValuesExtractors::Tensor<2> &tensor) const
{
  AssertIndexRange(tensor.first_tensor_component,
                   fe_values_views_cache.second_order_tensors.size());

  return fe_values_views_cache
    .second_order_tensors[tensor.first_tensor_component]
    .value_or_initialize([tensor, this]() {
      return FEValuesViews::Tensor<2, dim, spacedim>(
        *this, tensor.first_tensor_component);
    });
}



template <int dim, int spacedim>
inline const double &
FEValuesBase<dim, spacedim>::shape_value(const unsigned int i,
                                         const unsigned int q_point) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));
  Assert(fe->is_primitive(i), ExcShapeFunctionNotPrimitive(i));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->finite_element_output.shape_values(i, q_point);
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int row =
        this->finite_element_output
          .shape_function_to_row_table[i * fe->n_components() +
                                       fe->system_to_component_index(i).first];
      return this->finite_element_output.shape_values(row, q_point);
    }
}



template <int dim, int spacedim>
inline double
FEValuesBase<dim, spacedim>::shape_value_component(
  const unsigned int i,
  const unsigned int q_point,
  const unsigned int component) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_values,
         ExcAccessToUninitializedField("update_values"));
  AssertIndexRange(component, fe->n_components());
  Assert(present_cell.is_initialized(), ExcNotReinited());

  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return 0;

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int row =
    this->finite_element_output
      .shape_function_to_row_table[i * fe->n_components() + component];
  return this->finite_element_output.shape_values(row, q_point);
}



template <int dim, int spacedim>
inline const Tensor<1, spacedim> &
FEValuesBase<dim, spacedim>::shape_grad(const unsigned int i,
                                        const unsigned int q_point) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_gradients,
         ExcAccessToUninitializedField("update_gradients"));
  Assert(fe->is_primitive(i), ExcShapeFunctionNotPrimitive(i));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->finite_element_output.shape_gradients[i][q_point];
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int row =
        this->finite_element_output
          .shape_function_to_row_table[i * fe->n_components() +
                                       fe->system_to_component_index(i).first];
      return this->finite_element_output.shape_gradients[row][q_point];
    }
}



template <int dim, int spacedim>
inline Tensor<1, spacedim>
FEValuesBase<dim, spacedim>::shape_grad_component(
  const unsigned int i,
  const unsigned int q_point,
  const unsigned int component) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_gradients,
         ExcAccessToUninitializedField("update_gradients"));
  AssertIndexRange(component, fe->n_components());
  Assert(present_cell.is_initialized(), ExcNotReinited());
  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return Tensor<1, spacedim>();

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int row =
    this->finite_element_output
      .shape_function_to_row_table[i * fe->n_components() + component];
  return this->finite_element_output.shape_gradients[row][q_point];
}



template <int dim, int spacedim>
inline const Tensor<2, spacedim> &
FEValuesBase<dim, spacedim>::shape_hessian(const unsigned int i,
                                           const unsigned int q_point) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  Assert(fe->is_primitive(i), ExcShapeFunctionNotPrimitive(i));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->finite_element_output.shape_hessians[i][q_point];
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int row =
        this->finite_element_output
          .shape_function_to_row_table[i * fe->n_components() +
                                       fe->system_to_component_index(i).first];
      return this->finite_element_output.shape_hessians[row][q_point];
    }
}



template <int dim, int spacedim>
inline Tensor<2, spacedim>
FEValuesBase<dim, spacedim>::shape_hessian_component(
  const unsigned int i,
  const unsigned int q_point,
  const unsigned int component) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_hessians,
         ExcAccessToUninitializedField("update_hessians"));
  AssertIndexRange(component, fe->n_components());
  Assert(present_cell.is_initialized(), ExcNotReinited());
  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return Tensor<2, spacedim>();

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int row =
    this->finite_element_output
      .shape_function_to_row_table[i * fe->n_components() + component];
  return this->finite_element_output.shape_hessians[row][q_point];
}



template <int dim, int spacedim>
inline const Tensor<3, spacedim> &
FEValuesBase<dim, spacedim>::shape_3rd_derivative(
  const unsigned int i,
  const unsigned int q_point) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_3rd_derivatives,
         ExcAccessToUninitializedField("update_3rd_derivatives"));
  Assert(fe->is_primitive(i), ExcShapeFunctionNotPrimitive(i));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  // if the entire FE is primitive,
  // then we can take a short-cut:
  if (fe->is_primitive())
    return this->finite_element_output.shape_3rd_derivatives[i][q_point];
  else
    {
      // otherwise, use the mapping
      // between shape function
      // numbers and rows. note that
      // by the assertions above, we
      // know that this particular
      // shape function is primitive,
      // so we can call
      // system_to_component_index
      const unsigned int row =
        this->finite_element_output
          .shape_function_to_row_table[i * fe->n_components() +
                                       fe->system_to_component_index(i).first];
      return this->finite_element_output.shape_3rd_derivatives[row][q_point];
    }
}



template <int dim, int spacedim>
inline Tensor<3, spacedim>
FEValuesBase<dim, spacedim>::shape_3rd_derivative_component(
  const unsigned int i,
  const unsigned int q_point,
  const unsigned int component) const
{
  AssertIndexRange(i, fe->n_dofs_per_cell());
  Assert(this->update_flags & update_3rd_derivatives,
         ExcAccessToUninitializedField("update_3rd_derivatives"));
  AssertIndexRange(component, fe->n_components());
  Assert(present_cell.is_initialized(), ExcNotReinited());
  // check whether the shape function
  // is non-zero at all within
  // this component:
  if (fe->get_nonzero_components(i)[component] == false)
    return Tensor<3, spacedim>();

  // look up the right row in the
  // table and take the data from
  // there
  const unsigned int row =
    this->finite_element_output
      .shape_function_to_row_table[i * fe->n_components() + component];
  return this->finite_element_output.shape_3rd_derivatives[row][q_point];
}



template <int dim, int spacedim>
inline const FiniteElement<dim, spacedim> &
FEValuesBase<dim, spacedim>::get_fe() const
{
  return *fe;
}



template <int dim, int spacedim>
inline const Mapping<dim, spacedim> &
FEValuesBase<dim, spacedim>::get_mapping() const
{
  return *mapping;
}



template <int dim, int spacedim>
inline UpdateFlags
FEValuesBase<dim, spacedim>::get_update_flags() const
{
  return this->update_flags;
}



template <int dim, int spacedim>
inline const std::vector<Point<spacedim>> &
FEValuesBase<dim, spacedim>::get_quadrature_points() const
{
  Assert(this->update_flags & update_quadrature_points,
         ExcAccessToUninitializedField("update_quadrature_points"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.quadrature_points;
}



template <int dim, int spacedim>
inline const std::vector<double> &
FEValuesBase<dim, spacedim>::get_JxW_values() const
{
  Assert(this->update_flags & update_JxW_values,
         ExcAccessToUninitializedField("update_JxW_values"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.JxW_values;
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<1, dim, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobians() const
{
  Assert(this->update_flags & update_jacobians,
         ExcAccessToUninitializedField("update_jacobians"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobians;
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<2, dim, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_grads() const
{
  Assert(this->update_flags & update_jacobian_grads,
         ExcAccessToUninitializedField("update_jacobians_grads"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_grads;
}



template <int dim, int spacedim>
inline const Tensor<3, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_pushed_forward_grad(
  const unsigned int q_point) const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_grads,
         ExcAccessToUninitializedField("update_jacobian_pushed_forward_grads"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_pushed_forward_grads[q_point];
}



template <int dim, int spacedim>
inline const std::vector<Tensor<3, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_pushed_forward_grads() const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_grads,
         ExcAccessToUninitializedField("update_jacobian_pushed_forward_grads"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_pushed_forward_grads;
}



template <int dim, int spacedim>
inline const DerivativeForm<3, dim, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_2nd_derivative(
  const unsigned int q_point) const
{
  Assert(this->update_flags & update_jacobian_2nd_derivatives,
         ExcAccessToUninitializedField("update_jacobian_2nd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_2nd_derivatives[q_point];
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<3, dim, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_2nd_derivatives() const
{
  Assert(this->update_flags & update_jacobian_2nd_derivatives,
         ExcAccessToUninitializedField("update_jacobian_2nd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_2nd_derivatives;
}



template <int dim, int spacedim>
inline const Tensor<4, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_pushed_forward_2nd_derivative(
  const unsigned int q_point) const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_2nd_derivatives,
         ExcAccessToUninitializedField(
           "update_jacobian_pushed_forward_2nd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_pushed_forward_2nd_derivatives[q_point];
}



template <int dim, int spacedim>
inline const std::vector<Tensor<4, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_pushed_forward_2nd_derivatives() const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_2nd_derivatives,
         ExcAccessToUninitializedField(
           "update_jacobian_pushed_forward_2nd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_pushed_forward_2nd_derivatives;
}



template <int dim, int spacedim>
inline const DerivativeForm<4, dim, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_3rd_derivative(
  const unsigned int q_point) const
{
  Assert(this->update_flags & update_jacobian_3rd_derivatives,
         ExcAccessToUninitializedField("update_jacobian_3rd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_3rd_derivatives[q_point];
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<4, dim, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_3rd_derivatives() const
{
  Assert(this->update_flags & update_jacobian_3rd_derivatives,
         ExcAccessToUninitializedField("update_jacobian_3rd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_3rd_derivatives;
}



template <int dim, int spacedim>
inline const Tensor<5, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_pushed_forward_3rd_derivative(
  const unsigned int q_point) const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_3rd_derivatives,
         ExcAccessToUninitializedField(
           "update_jacobian_pushed_forward_3rd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_pushed_forward_3rd_derivatives[q_point];
}



template <int dim, int spacedim>
inline const std::vector<Tensor<5, spacedim>> &
FEValuesBase<dim, spacedim>::get_jacobian_pushed_forward_3rd_derivatives() const
{
  Assert(this->update_flags & update_jacobian_pushed_forward_3rd_derivatives,
         ExcAccessToUninitializedField(
           "update_jacobian_pushed_forward_3rd_derivatives"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.jacobian_pushed_forward_3rd_derivatives;
}



template <int dim, int spacedim>
inline const std::vector<DerivativeForm<1, spacedim, dim>> &
FEValuesBase<dim, spacedim>::get_inverse_jacobians() const
{
  Assert(this->update_flags & update_inverse_jacobians,
         ExcAccessToUninitializedField("update_inverse_jacobians"));
  Assert(present_cell.is_initialized(), ExcNotReinited());
  return this->mapping_output.inverse_jacobians;
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEValuesBase<dim, spacedim>::dof_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    0U, dofs_per_cell);
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEValuesBase<dim, spacedim>::dof_indices_starting_at(
  const unsigned int start_dof_index) const
{
  Assert(start_dof_index <= dofs_per_cell,
         ExcIndexRange(start_dof_index, 0, dofs_per_cell + 1));
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    start_dof_index, dofs_per_cell);
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEValuesBase<dim, spacedim>::dof_indices_ending_at(
  const unsigned int end_dof_index) const
{
  Assert(end_dof_index < dofs_per_cell,
         ExcIndexRange(end_dof_index, 0, dofs_per_cell));
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    0U, end_dof_index + 1);
}



template <int dim, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FEValuesBase<dim, spacedim>::quadrature_point_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    0U, n_quadrature_points);
}



template <int dim, int spacedim>
inline const Point<spacedim> &
FEValuesBase<dim, spacedim>::quadrature_point(const unsigned int q_point) const
{
  Assert(this->update_flags & update_quadrature_points,
         ExcAccessToUninitializedField("update_quadrature_points"));
  AssertIndexRange(q_point, this->mapping_output.quadrature_points.size());
  Assert(present_cell.is_initialized(), ExcNotReinited());

  return this->mapping_output.quadrature_points[q_point];
}



template <int dim, int spacedim>
inline double
FEValuesBase<dim, spacedim>::JxW(const unsigned int q_point) const
{
  Assert(this->update_flags & update_JxW_values,
         ExcAccessToUninitializedField("update_JxW_values"));
  AssertIndexRange(q_point, this->mapping_output.JxW_values.size());
  Assert(present_cell.is_initialized(), ExcNotReinited());

  return this->mapping_output.JxW_values[q_point];
}



template <int dim, int spacedim>
inline const DerivativeForm<1, dim, spacedim> &
FEValuesBase<dim, spacedim>::jacobian(const unsigned int q_point) const
{
  Assert(this->update_flags & update_jacobians,
         ExcAccessToUninitializedField("update_jacobians"));
  AssertIndexRange(q_point, this->mapping_output.jacobians.size());
  Assert(present_cell.is_initialized(), ExcNotReinited());

  return this->mapping_output.jacobians[q_point];
}



template <int dim, int spacedim>
inline const DerivativeForm<2, dim, spacedim> &
FEValuesBase<dim, spacedim>::jacobian_grad(const unsigned int q_point) const
{
  Assert(this->update_flags & update_jacobian_grads,
         ExcAccessToUninitializedField("update_jacobians_grads"));
  AssertIndexRange(q_point, this->mapping_output.jacobian_grads.size());
  Assert(present_cell.is_initialized(), ExcNotReinited());

  return this->mapping_output.jacobian_grads[q_point];
}



template <int dim, int spacedim>
inline const DerivativeForm<1, spacedim, dim> &
FEValuesBase<dim, spacedim>::inverse_jacobian(const unsigned int q_point) const
{
  Assert(this->update_flags & update_inverse_jacobians,
         ExcAccessToUninitializedField("update_inverse_jacobians"));
  AssertIndexRange(q_point, this->mapping_output.inverse_jacobians.size());
  Assert(present_cell.is_initialized(), ExcNotReinited());

  return this->mapping_output.inverse_jacobians[q_point];
}



template <int dim, int spacedim>
inline const Tensor<1, spacedim> &
FEValuesBase<dim, spacedim>::normal_vector(const unsigned int q_point) const
{
  Assert(this->update_flags & update_normal_vectors,
         (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
           "update_normal_vectors")));
  AssertIndexRange(q_point, this->mapping_output.normal_vectors.size());
  Assert(present_cell.is_initialized(), ExcNotReinited());

  return this->mapping_output.normal_vectors[q_point];
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
