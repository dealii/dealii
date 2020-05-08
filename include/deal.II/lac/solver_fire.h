// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_solver_fire_h
#define dealii_solver_fire_h


#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>

#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/solver.h>

#include <functional>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup Solvers */
/*@{*/

/**
 * FIRE (Fast Inertial Relaxation Engine) for minimization of (potentially
 * non-linear) objective function $E(\mathbf x)$, $\mathbf x$ is a vector of
 * $n$ variables ($n$ is the number of variables of the objective function).
 * Like all other solver classes, it can work on any kind of vector and matrix
 * as long as they satisfy certain requirements (for the requirements on
 * matrices and vectors in order to work with this class, see the documentation
 * of the Solver base class). The type of the solution vector must be passed as
 * template argument, and defaults to dealii::Vector<double>.
 *
 * FIRE is a damped dynamics method described in
 * <a href="https://doi.org/10.1103/PhysRevLett.97.170201">Structural
 * Relaxation Made Simple</a> by Bitzek et al. 2006, typically used to find
 * stable equilibrium configurations of atomistic systems in computational
 * material science. Starting from a given initial configuration of the
 * atomistic system, the algorithm relies on inertia to obtain (nearest)
 * configuration with least potential energy.
 *
 * Notation:
 *  - The global vector of unknown variables: $\mathbf x$.
 *  - Objective function:                     $E(\mathbf x)$.
 *  - Rate of change of unknowns:             $\mathbf v$.
 *  - Gradient of the objective
 *    function w.r.t unknowns:                $\mathbf g = \nabla E(\mathbf x)$.
 *  - Mass matrix:                            $\mathbf M$.
 *  - Initial guess of unknowns:              $\mathbf x_0$.
 *  - Time step:                              $\Delta t$.
 *
 * Given initial values for $\Delta t$, $\alpha = \alpha_0$, $\epsilon$,
 * $\mathbf x = \mathbf x_0$ and $\mathbf v= \mathbf 0$ along with a given mass
 * matrix $\mathbf M$, FIRE algorithm is as follows,
 * 1. Calculate $\mathbf g = \nabla E(\mathbf x)$ and check for convergence
 *    ($\mathbf g \cdot \mathbf g < \epsilon^2 $).
 * 2. Update $\mathbf x$ and $V$ using simple (forward) Euler integration step,
 *    <BR>
 *        $\mathbf x = \mathbf x + \Delta t \mathbf v$,                 <BR>
 *        $\mathbf v = \mathbf v + \Delta t \mathbf M^{-1} \cdot \mathbf g$.
 * 3. Calculate $p = \mathbf g \cdot \mathbf v$.
 * 4. Set $\mathbf v = (1-\alpha) \mathbf v
 *                   + \alpha \frac{|\mathbf v|}{|\mathbf g|} \mathbf g$.
 * 5. If $p<0$ and number of steps since $p$ was last negative is larger
 *    than certain value, then increase time step $\Delta t$ and decrease
 *    $\alpha$.
 * 6. If $p>0$, then decrease the time step, freeze the system i.e.,
 *    $\mathbf v = \mathbf 0$ and reset $\alpha = \alpha_0$.
 * 7. Return to 1.
 *
 * Also see
 * <a href="http://onlinelibrary.wiley.com/doi/10.1002/pamm.201110246/full">
 * Energy-Minimization in Atomic-to-Continuum Scale-Bridging Methods </a> by
 * Eidel et al. 2011.
 *
 * @author Vishal Boddu, Denis Davydov, 2017
 */
template <typename VectorType = Vector<double>>
class SolverFIRE : public SolverBase<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the initial time step for the (forward)
     * Euler integration step to 0.1, the maximum time step to 1 and the
     * maximum change allowed in any variable (per iteration) to 1.
     */
    explicit AdditionalData(const double initial_timestep    = 0.1,
                            const double maximum_timestep    = 1,
                            const double maximum_linfty_norm = 1);

    /**
     * Initial time step for the (forward) Euler integration step.
     */
    const double initial_timestep;

    /**
     * Maximum time step for the (forward) Euler integration step.
     */
    const double maximum_timestep;

    /**
     * Maximum change allowed in any variable of the objective function.
     */
    const double maximum_linfty_norm;
  };

  /**
   * Constructor.
   */
  SolverFIRE(SolverControl &           solver_control,
             VectorMemory<VectorType> &vector_memory,
             const AdditionalData &    data = AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverFIRE(SolverControl &       solver_control,
             const AdditionalData &data = AdditionalData());

  /**
   * Obtain a set of variables @p x that minimize an objective function
   * described by the polymorphic function wrapper @p compute, with a given
   * preconditioner @p inverse_mass_matrix and initial @p x values.
   * The function @p compute returns the objective function's value and updates
   * the objective function's gradient (with respect to the variables) when
   * passed in as first argument based on the second argument-- the state of
   * variables.
   */
  template <typename PreconditionerType = DiagonalMatrix<VectorType>>
  void
  solve(const std::function<double(VectorType &, const VectorType &)> &compute,
        VectorType &                                                   x,
        const PreconditionerType &inverse_mass_matrix);

  /**
   * Solve for x that minimizes $E(\mathbf x)$ for the <EM>special case</EM>
   * when $E(\mathbf x)
   * = \frac{1}{2} \mathbf x^{T} \mathbf A \mathbf x - \mathbf x^{T} \mathbf b$.
   */
  template <typename MatrixType, typename PreconditionerType>
  void
  solve(const MatrixType &        A,
        VectorType &              x,
        const VectorType &        b,
        const PreconditionerType &preconditioner);

protected:
  /**
   * Interface for derived class. This function gets the current iteration
   * @p x (variables), @p v (x's time derivative) and @p g (the gradient) in
   * each step.
   * It can be used for graphical output of the convergence history.
   */
  virtual void
  print_vectors(const unsigned int,
                const VectorType &x,
                const VectorType &v,
                const VectorType &g) const;

  /**
   * Additional data to the solver.
   */
  const AdditionalData additional_data;
};

/*@}*/

/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN

template <typename VectorType>
SolverFIRE<VectorType>::AdditionalData::AdditionalData(
  const double initial_timestep,
  const double maximum_timestep,
  const double maximum_linfty_norm)
  : initial_timestep(initial_timestep)
  , maximum_timestep(maximum_timestep)
  , maximum_linfty_norm(maximum_linfty_norm)
{
  AssertThrow(initial_timestep > 0. && maximum_timestep > 0. &&
                maximum_linfty_norm > 0.,
              ExcMessage("Expected positive values for initial_timestep, "
                         "maximum_timestep and maximum_linfty_norm but one "
                         "or more of the these values are not positive."));
}



template <typename VectorType>
SolverFIRE<VectorType>::SolverFIRE(SolverControl &           solver_control,
                                   VectorMemory<VectorType> &vector_memory,
                                   const AdditionalData &    data)
  : SolverBase<VectorType>(solver_control, vector_memory)
  , additional_data(data)
{}



template <typename VectorType>
SolverFIRE<VectorType>::SolverFIRE(SolverControl &       solver_control,
                                   const AdditionalData &data)
  : SolverBase<VectorType>(solver_control)
  , additional_data(data)
{}



template <typename VectorType>
template <typename PreconditionerType>
void
SolverFIRE<VectorType>::solve(
  const std::function<double(VectorType &, const VectorType &)> &compute,
  VectorType &                                                   x,
  const PreconditionerType &inverse_mass_matrix)
{
  LogStream::Prefix prefix("FIRE");

  // FIRE algorithm constants
  const double DELAYSTEP       = 5;
  const double TIMESTEP_GROW   = 1.1;
  const double TIMESTEP_SHRINK = 0.5;
  const double ALPHA_0         = 0.1;
  const double ALPHA_SHRINK    = 0.99;

  using real_type = typename VectorType::real_type;

  typename VectorMemory<VectorType>::Pointer v(this->memory);
  typename VectorMemory<VectorType>::Pointer g(this->memory);

  // Set velocities to zero but not gradients
  // as we are going to compute them soon.
  v->reinit(x, false);
  g->reinit(x, true);

  // Refer to v and g with some readable names.
  VectorType &velocities = *v;
  VectorType &gradients  = *g;

  // Update gradients for the new x.
  compute(gradients, x);

  unsigned int iter = 0;

  SolverControl::State conv = SolverControl::iterate;
  conv = this->iteration_status(iter, gradients * gradients, x);
  if (conv != SolverControl::iterate)
    return;

  // Refer to additional data members with some readable names.
  const auto &maximum_timestep = additional_data.maximum_timestep;
  double      timestep         = additional_data.initial_timestep;

  // First scaling factor.
  double alpha = ALPHA_0;

  unsigned int previous_iter_with_positive_v_dot_g = 0;

  while (conv == SolverControl::iterate)
    {
      ++iter;
      // Euler integration step.
      x.add(timestep, velocities);                     // x += dt     * v
      inverse_mass_matrix.vmult(gradients, gradients); // g  = M^{-1} * g
      velocities.add(-timestep, gradients);            // v -= dt     * h

      // Compute gradients for the new x.
      compute(gradients, x);

      const real_type gradient_norm_squared = gradients * gradients;
      conv = this->iteration_status(iter, gradient_norm_squared, x);
      if (conv != SolverControl::iterate)
        break;

      // v_dot_g = V * G
      const real_type v_dot_g = velocities * gradients;

      if (v_dot_g < 0.)
        {
          const real_type velocities_norm_squared = velocities * velocities;

          // Check if we divide by zero in DEBUG mode.
          Assert(gradient_norm_squared > 0., ExcInternalError());

          // beta = - alpha |V|/|G|
          const real_type beta =
            -alpha * std::sqrt(velocities_norm_squared / gradient_norm_squared);

          // V = (1-alpha) V + beta G.
          velocities.sadd(1. - alpha, beta, gradients);

          if (iter - previous_iter_with_positive_v_dot_g > DELAYSTEP)
            {
              // Increase timestep and decrease alpha.
              timestep = std::min(timestep * TIMESTEP_GROW, maximum_timestep);
              alpha *= ALPHA_SHRINK;
            }
        }
      else
        {
          // Decrease timestep, reset alpha and set V = 0.
          previous_iter_with_positive_v_dot_g = iter;
          timestep *= TIMESTEP_SHRINK;
          alpha      = ALPHA_0;
          velocities = 0.;
        }

      real_type vmax = velocities.linfty_norm();

      // Change timestep if any dof would move more than maximum_linfty_norm.
      if (vmax > 0.)
        {
          const double minimal_timestep =
            additional_data.maximum_linfty_norm / vmax;
          if (minimal_timestep < timestep)
            timestep = minimal_timestep;
        }

      print_vectors(iter, x, velocities, gradients);

    } // While we need to iterate.

  // In the case of failure: throw exception.
  if (conv != SolverControl::success)
    AssertThrow(false,
                SolverControl::NoConvergence(iter, gradients * gradients));
}



template <typename VectorType>
template <typename MatrixType, typename PreconditionerType>
void
SolverFIRE<VectorType>::solve(const MatrixType &        A,
                              VectorType &              x,
                              const VectorType &        b,
                              const PreconditionerType &preconditioner)
{
  std::function<double(VectorType &, const VectorType &)> compute_func =
    [&](VectorType &g, const VectorType &x) -> double {
    // Residual of the quadratic form $ \frac{1}{2} xAx - xb $.
    // G = b - Ax
    A.residual(g, x, b);

    // Gradient G = Ax -b.
    g *= -1.;

    // The quadratic form $\frac{1}{2} xAx - xb $.
    return 0.5 * A.matrix_norm_square(x) - x * b;
  };

  this->solve(compute_func, x, preconditioner);
}



template <typename VectorType>
void
SolverFIRE<VectorType>::print_vectors(const unsigned int,
                                      const VectorType &,
                                      const VectorType &,
                                      const VectorType &) const
{}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
