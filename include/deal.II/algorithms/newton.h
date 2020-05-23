// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2019 by the deal.II authors
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


#ifndef dealii_newton_h
#define dealii_newton_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>
#include <deal.II/algorithms/operator.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class ParameterHandler;
#endif

namespace Algorithms
{
  /**
   * Operator class performing Newton's iteration with standard step size
   * control and adaptive matrix generation.
   *
   * This class performs a Newton iteration up to convergence determined by
   * #control. If after an update the norm of the residual has become larger,
   * then step size control is activated and the update is subsequently
   * divided by two until the residual actually becomes smaller (or the
   * minimal scaling factor determined by #n_stepsize_iterations is reached).
   *
   * Since assembling matrices, depending on the implementation, tends to be
   * costly, this method applies an adaptive reassembling strategy. Only if
   * the reduction factor for the residual is more than #threshold, the event
   * Algorithms::bad_derivative is submitted to #inverse_derivative. It is up
   * to this object to implement reassembling accordingly.
   *
   * <h3>Contents of the AnyData objects</h3>
   *
   * The only value used by the Newton method is the first vector in the
   * parameter <tt>out</tt> of operator()(). It serves as the start vector of
   * Newton's method and in the end contains the solution. All other vectors
   * of <tt>out</tt> are ignored by Newton's method and its inner Operator
   * objects. All vectors of <tt>in</tt> are forwarded to the inner Operator
   * objects, with additional information added as follows.
   *
   * When calling (*#residual)(), the AnyData <tt>in</tt> given to the Newton
   * iteration is prepended by a vector <tt>"Newton iterate"</tt>, the current
   * value of the Newton iterate, which can be used to evaluate the residual
   * at this point.
   *
   * For the call to (*#inverse_derivative), the vector <tt>"Newton
   * residual"</tt> is inserted before <tt>"Newton iterate"</tt>.
   *
   * @author Guido Kanschat, 2006, 2010
   */
  template <typename VectorType>
  class Newton : public OperatorBase
  {
  public:
    /**
     * Constructor, receiving the applications computing the residual and
     * solving the linear problem, respectively.
     */
    Newton(OperatorBase &residual, OperatorBase &inverse_derivative);

    /**
     * Declare the parameters applicable to Newton's method.
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * Read the parameters in the ParameterHandler.
     */
    void
    parse_parameters(ParameterHandler &param);

    /**
     * Initialize the pointer data_out for debugging.
     */
    void
    initialize(OutputOperator<VectorType> &output);

    /**
     * The actual Newton iteration. The initial value is in <tt>out(0)</tt>,
     * which also contains the result after convergence. Values in <tt>in</tt>
     * are not used by Newton, but will be handed down to the objects
     * #residual and #inverse_derivative.
     */
    virtual void
    operator()(AnyData &out, const AnyData &in) override;

    virtual void
    notify(const Event &) override;

    /**
     * Set the maximal residual reduction allowed without triggering
     * assembling in the next step. Return the previous value.
     */
    double
    threshold(double new_value);

    /**
     * Control object for the Newton iteration.
     */
    ReductionControl control;

  private:
    /**
     * The operator computing the residual.
     */
    SmartPointer<OperatorBase, Newton<VectorType>> residual;

    /**
     * The operator applying the inverse derivative to the residual.
     */
    SmartPointer<OperatorBase, Newton<VectorType>> inverse_derivative;

    /**
     * The operator handling the output in case the debug_vectors is true.
     * Call the initialize function first.
     */
    SmartPointer<OutputOperator<VectorType>, Newton<VectorType>> data_out;

    /**
     * This flag is set by the function assemble(), indicating that the matrix
     * must be assembled anew upon start.
     */
    bool assemble_now;

    /**
     * A flag used to decide how many stepsize iteration should be made.
     * Default is the original value of 21.
     *
     * Enter zero here to turn off stepsize control.
     *
     * @note Controlled by <tt>Stepsize iterations</tt> in parameter file
     */
    unsigned int n_stepsize_iterations;

    /**
     * Threshold for re-assembling matrix.
     *
     * If the quotient of two consecutive residuals is smaller than this
     * threshold, the system matrix is not assembled in this step.
     *
     * @note This parameter should be adjusted to the residual gain of the
     * inner solver.
     *
     * The default values is zero, resulting in reassembling in every Newton
     * step.
     */
    double assemble_threshold;

  public:
    /**
     * Print residual, update and updated solution after each step into file
     * <tt>Newton_NNN</tt>?
     */
    bool debug_vectors;
    /**
     * Write debug output to @p deallog; the higher the number, the more
     * output.
     */
    unsigned int debug;
  };
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif
