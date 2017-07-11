#ifndef dealii__solver_fire_h
#define dealii__solver_fire_h


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
 * non-linear) objective function $E(U)$, $U$ is a vector of $n$ variables
 * ($n$ is the number of variables of the objective function).
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
 * The problem of solving the system of linear equations $AU = F$ for $U$
 * can be re-casted into the problem of finding $U$ for which the quadratic
 * function $\frac{1}{2} U^{T}AU - U^{T}F$ is minimized.
 *
 * Notation:
 * The global vector of unknowns:   $U$.        <BR>
 * Objective function:              $E(U)$.     <BR>
 * Rate of change of unknowns:      $V$.        <BR>
 * Gradient of the objective
 * function w.r.t unknowns:         $G$.        <BR>
 * Mass matrix:                     $M$.        <BR>
 * Initial guess of unknowns:       $U_0$.      <BR>
 * Time step:                       $\Delta t$. <BR>
 *
 * Given initial values for $\Delta t$, $\alpha = \alpha_0$, $U = U_0$ and $V=0$
 * along with a given mass matrix $M$,
 * FIRE algorithm is as follows,
 * 1. Calculate $G = \nabla E(U)$ and check for convergence.
 *    Calculate $U$ and $V$ using any common Molecular Dynamics time integrator.
 *    Here we use simple Euler integration step, <BR>
 *        $U = U + \Delta t V$,                  <BR>
 *        $V = V + \Delta t M^{-1} \cdot G$.
 * 2. Calculate $P = G \cdot V$.
 * 3. Set $V = (1-\alpha) V + \alpha \frac{|V|}{|G|} G$.
 * 4. If $P \leq 0$ and number of steps since P was non-negative is larger than
 *    certain value, then increase time step $\Delta t$ and decrease $\alpha$.
 * 5. If $P>0$, then decrease the time step, freeze the system i.e., V = 0 and
 *    reset $\alpha = \alpha_0$.
 * 6. Return to 1.
 *
 * Alse see
 * <a href="http://onlinelibrary.wiley.com/doi/10.1002/pamm.201110246/full">
 * Energy-Minimization in Atomic-to-Continuum Scale-Bridging Methods </a> by
 * Eidel et al. 2011.
 *
 * @author Vishal Boddu 2017
 */
template<typename VectorType = Vector<double> >
class SolverFIRE : public Solver<VectorType>
{

public:

  struct AdditionalData
  {
    explicit
    AdditionalData (const double  initial_timestep    = 0.1,
                    const double  maximum_timestep    = 1,
                    const double  maximum_linfty_norm = 1);

    /**
     * Initial time step.
     */
    const double initial_timestep;

    /**
     * Maximum time step.
     */
    const double maximum_timestep;

    /**
     * Maximum change allowed in any degree of freedom.
     */
    const double maximum_linfty_norm;

  };

  /**
   * Constructor.
   */
  SolverFIRE (SolverControl            &solver_control,
              VectorMemory<VectorType> &vector_memory,
              const AdditionalData     &data          );

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverFIRE (SolverControl         &solver_control,
              const AdditionalData  &data          );

  /**
   * Virtual destructor.
   */
  virtual ~SolverFIRE();

  /**
   * Obtain a set of #p u (variables) that minimize an objective function
   * described by the polymorphic function wrapper @p compute, with a given
   * preconditioner @p inverse_masses and initial @p u values.
   * The function @p compute takes in the state of the (u) variables as
   * argument and returns a pair of objective function's value and
   * objective function's gradient (with respect to the variables).
   */
  template<typename PreconditionerType = DiagonalMatrix<VectorType>>
  void solve
  (std::function<double(VectorType &, const VectorType &)>  compute,
   VectorType                                              &u,
   const PreconditionerType                                &inverse_masses);

  /**
   * Solve for x that minimizes the quadratic function
   * $E(x) = \frac{1}{2} x^{T}Ax - x^{T}b$.
   */
  template<typename MatrixType, typename PreconditionerType>
  void solve (const MatrixType         &A,
              VectorType               &x,
              const VectorType         &b,
              const PreconditionerType &precondition);

protected:

  /**
   * Interface for derived class. This function gets the current iteration
   * u, u's time derivative and the gradient in each step. It can be used
   * for a graphical output of the convergence history.
   */
  void print_vectors (const unsigned int,
                      const VectorType &,
                      const VectorType &,
                      const VectorType &) const;

  /**
   * Additional parameters.
   */
  const AdditionalData additional_data;

};

/*@}*/

/*------------------------- Implementation ----------------------------*/

#ifndef DOXYGEN

template<typename VectorType>
SolverFIRE<VectorType>::AdditionalData::
AdditionalData (const double  initial_timestep,
                const double  maximum_timestep,
                const double  maximum_linfty_norm)
  :
  initial_timestep(initial_timestep),
  maximum_timestep(maximum_timestep),
  maximum_linfty_norm(maximum_linfty_norm)
{
  AssertThrow (initial_timestep    > 0. &&
               maximum_timestep    > 0. &&
               maximum_linfty_norm > 0.,
               ExcMessage("Expected positive values for initial_timestep, "
                          "maximum_timestep and maximum_linfty_norm but one "
                          "or more of the these values are not positive."));
}



template<typename VectorType>
SolverFIRE<VectorType>::
SolverFIRE (SolverControl            &solver_control,
            VectorMemory<VectorType> &vector_memory,
            const AdditionalData     &data          )
  :
  Solver<VectorType>(solver_control, vector_memory),
  additional_data(data)
{}



template<typename VectorType>
SolverFIRE<VectorType>::
SolverFIRE (SolverControl         &solver_control,
            const AdditionalData  &data          )
  :
  Solver<VectorType>(solver_control),
  additional_data(data)
{}



template<typename VectorType>
SolverFIRE<VectorType>::~SolverFIRE()
{}



template<typename VectorType>
template<typename PreconditionerType>
void
SolverFIRE<VectorType>::solve
(std::function<double(VectorType &, const VectorType &)>  compute,
 VectorType                                              &u,
 const PreconditionerType                                &inverse_masses)
{
  deallog.push("FIRE");

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
  v->reinit(u,false);
  g->reinit(u,true);

  // Refer to v and g with some readable names.
  VectorType &velocities = *v;
  VectorType &gradients  = *g;

  // Update gradients for the new u.
  compute(gradients, u);

  unsigned int iter = 0;

  SolverControl::State conv = SolverControl::iterate;
  conv = this->iteration_status (iter, gradients * gradients, u);
  if (conv != SolverControl::iterate)
    return;

  // Refer to additional data members with some readable names.
  const auto &maximum_timestep   = additional_data.maximum_timestep;
  double timestep                = additional_data.initial_timestep;

  // First scaling factor.
  double alpha = ALPHA_0;

  unsigned int previous_iter_with_positive_v_dot_g = 0;

  while (conv == SolverControl::iterate)
    {
      ++iter;
      // Euler integration step.
      u.add (timestep, velocities);                // U += dt     * V
      inverse_masses.vmult(gradients, gradients);  // G  = M^{-1} * G
      velocities.add (-timestep, gradients);       // V -= dt     * G

      // Compute gradients for the new u.
      compute(gradients, u);

      const real_type gradient_norm_squared = gradients * gradients;
      conv = this->iteration_status(iter, gradient_norm_squared, u);
      if (conv != SolverControl::iterate)
        break;

      // v_dot_g = V * G
      const real_type v_dot_g = velocities * gradients;

      if (v_dot_g < 0.)
        {
          const real_type velocities_norm_squared =
            velocities * velocities;

          // Check if we divide by zero in DEBUG mode.
          Assert (gradient_norm_squared > 0., ExcInternalError());

          // beta = - alpha |V|/|G|
          const real_type beta = -alpha *
                                 std::sqrt (velocities_norm_squared
                                            /
                                            gradient_norm_squared);

          // V = (1-alpha) V + beta G.
          velocities.sadd (1. - alpha, beta, gradients);

          if (iter - previous_iter_with_positive_v_dot_g > DELAYSTEP)
            {
              // Increase timestep and decrease alpha.
              timestep = std::min (timestep*TIMESTEP_GROW, maximum_timestep);
              alpha *= ALPHA_SHRINK;
            }
        }
      else
        {
          // Decrease timestep, reset alpha and set V = 0.
          previous_iter_with_positive_v_dot_g = iter;
          timestep *= TIMESTEP_SHRINK;
          alpha = ALPHA_0;
          velocities = 0.;
        }

      real_type vmax = velocities.linfty_norm();

      // Change timestep if any dof would move more than maximum_linfty_norm.
      if (vmax > 0.)
        {
          const double minimal_timestep = additional_data.maximum_linfty_norm
                                          /
                                          vmax;
          if (minimal_timestep < timestep)
            timestep = minimal_timestep;
        }

      print_vectors(iter, u, velocities, gradients);

    } // While we need to iterate.

  deallog.pop();

  // In the case of failure: throw exception.
  if (conv != SolverControl::success)
    AssertThrow (false,
                 SolverControl::NoConvergence (iter, gradients * gradients));

}



template <typename VectorType>
template<typename MatrixType, typename PreconditionerType>
void SolverFIRE<VectorType>::solve (const MatrixType         &A,
                                    VectorType               &x,
                                    const VectorType         &b,
                                    const PreconditionerType &precondition)
{

  std::function<double(VectorType &,  const VectorType &)> compute_func =
    [&]               (decltype(x) &G, decltype(b)       &x) -> double
  {
    // Residual of the quadratic form $ \frac{1}{2} xAx - xb $.
    // G = b - Ax
    A.residual(G, x, b);

    // Gradient G = Ax -b.
    G *= -1.;

    // The quadratic form $ \frac{1}{2} x^t A x - x^{t} b $.
    return 0.5*A.matrix_norm_square(x) - x*b;
  };

  this->solve (compute_func, x, precondition);
}



template <typename VectorType>
void
SolverFIRE<VectorType>::
print_vectors (const unsigned int,
               const VectorType &,
               const VectorType &,
               const VectorType &) const
{}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
