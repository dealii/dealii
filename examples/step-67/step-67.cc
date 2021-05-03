/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Martin Kronbichler, 2020
 */

// The include files are similar to the previous matrix-free tutorial programs
// step-37, step-48, and step-59
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/time_stepping.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <iostream>

// The following file includes the CellwiseInverseMassMatrix data structure
// that we will use for the mass matrix inversion, the only new include
// file for this tutorial program:
#include <deal.II/matrix_free/operators.h>



namespace Euler_DG
{
  using namespace dealii;

  // Similarly to the other matrix-free tutorial programs, we collect all
  // parameters that control the execution of the program at the top of the
  // file. Besides the dimension and polynomial degree we want to run with, we
  // also specify a number of points in the Gaussian quadrature formula we
  // want to use for the nonlinear terms in the Euler equations. Furthermore,
  // we specify the time interval for the time-dependent problem, and
  // implement two different test cases. The first one is an analytical
  // solution in 2D, whereas the second is a channel flow around a cylinder as
  // described in the introduction. Depending on the test case, we also change
  // the final time up to which we run the simulation, and a variable
  // `output_tick` that specifies in which intervals we want to write output
  // (assuming that the tick is larger than the time step size).
  constexpr unsigned int testcase             = 0;
  constexpr unsigned int dimension            = 2;
  constexpr unsigned int n_global_refinements = 3;
  constexpr unsigned int fe_degree            = 5;
  constexpr unsigned int n_q_points_1d        = fe_degree + 2;

  using Number = double;

  constexpr double gamma       = 1.4;
  constexpr double final_time  = testcase == 0 ? 10 : 2.0;
  constexpr double output_tick = testcase == 0 ? 1 : 0.05;

  // Next off are some details of the time integrator, namely a Courant number
  // that scales the time step size in terms of the formula $\Delta t =
  // \text{Cr} n_\text{stages} \frac{h}{(p+1)^{1.5} (\|\mathbf{u} +
  // c)_\text{max}}$, as well as a selection of a few low-storage Runge--Kutta
  // methods. We specify the Courant number per stage of the Runge--Kutta
  // scheme, as this gives a more realistic expression of the numerical cost
  // for schemes of various numbers of stages.
  const double courant_number = 0.15 / std::pow(fe_degree, 1.5);
  enum LowStorageRungeKuttaScheme
  {
    stage_3_order_3, /* Kennedy, Carpenter, Lewis, 2000 */
    stage_5_order_4, /* Kennedy, Carpenter, Lewis, 2000 */
    stage_7_order_4, /* Tselios, Simos, 2007 */
    stage_9_order_5, /* Kennedy, Carpenter, Lewis, 2000 */
  };
  constexpr LowStorageRungeKuttaScheme lsrk_scheme = stage_5_order_4;

  // Eventually, we select a detail of the spatial discretization, namely the
  // numerical flux (Riemann solver) at the faces between cells. For this
  // program, we have implemented a modified variant of the Lax--Friedrichs
  // flux and the Harten--Lax--van Leer (HLL) flux.
  enum EulerNumericalFlux
  {
    lax_friedrichs_modified,
    harten_lax_vanleer,
  };
  constexpr EulerNumericalFlux numerical_flux_type = lax_friedrichs_modified;



  // @sect3{Equation data}

  // We now define a class with the exact solution for the test case 0 and one
  // with a background flow field for test case 1 of the channel. Given that
  // the Euler equations are a problem with $d+2$ equations in $d$ dimensions,
  // we need to tell the Function base class about the correct number of
  // components.
  template <int dim>
  class ExactSolution : public Function<dim>
  {
  public:
    ExactSolution(const double time)
      : Function<dim>(dim + 2, time)
    {}

    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };



  // As far as the actual function implemented is concerned, the analytical
  // test case is an isentropic vortex case (see e.g. the book by Hesthaven
  // and Warburton, Example 6.1 in Section 6.6 on page 209) which fulfills the
  // Euler equations with zero force term on the right hand side. Given that
  // definition, we return either the density, the momentum, or the energy
  // depending on which component is requested. Note that the original
  // definition of the density involves the $\frac{1}{\gamma -1}$-th power of
  // some expression. Since `std::pow()` has pretty slow implementations on
  // some systems, we replace it by logarithm followed by exponentiation (of
  // base 2), which is mathematically equivalent but usually much better
  // optimized. This formula might lose accuracy in the last digits
  // for very small numbers compared to `std::pow()`, but we are happy with
  // it anyway, since small numbers map to data close to 1.
  //
  // For the channel test case, we simply select a density of 1, a velocity of
  // 0.4 in $x$ direction and zero in the other directions, and an energy that
  // corresponds to a speed of sound of 1.3 measured against the background
  // velocity field, computed from the relation $E = \frac{c^2}{\gamma (\gamma
  // -1)} + \frac 12 \rho \|u\|^2$.
  template <int dim>
  double ExactSolution<dim>::value(const Point<dim> & x,
                                   const unsigned int component) const
  {
    const double t = this->get_time();

    switch (testcase)
      {
        case 0:
          {
            Assert(dim == 2, ExcNotImplemented());
            const double beta = 5;

            Point<dim> x0;
            x0[0] = 5.;
            const double radius_sqr =
              (x - x0).norm_square() - 2. * (x[0] - x0[0]) * t + t * t;
            const double factor =
              beta / (numbers::PI * 2) * std::exp(1. - radius_sqr);
            const double density_log = std::log2(
              std::abs(1. - (gamma - 1.) / gamma * 0.25 * factor * factor));
            const double density = std::exp2(density_log * (1. / (gamma - 1.)));
            const double u       = 1. - factor * (x[1] - x0[1]);
            const double v       = factor * (x[0] - t - x0[0]);

            if (component == 0)
              return density;
            else if (component == 1)
              return density * u;
            else if (component == 2)
              return density * v;
            else
              {
                const double pressure =
                  std::exp2(density_log * (gamma / (gamma - 1.)));
                return pressure / (gamma - 1.) +
                       0.5 * (density * u * u + density * v * v);
              }
          }

        case 1:
          {
            if (component == 0)
              return 1.;
            else if (component == 1)
              return 0.4;
            else if (component == dim + 1)
              return 3.097857142857143;
            else
              return 0.;
          }

        default:
          Assert(false, ExcNotImplemented());
          return 0.;
      }
  }



  // @sect3{Low-storage explicit Runge--Kutta time integrators}

  // The next few lines implement a few low-storage variants of Runge--Kutta
  // methods. These methods have specific Butcher tableaux with coefficients
  // $b_i$ and $a_i$ as shown in the introduction. As usual in Runge--Kutta
  // method, we can deduce time steps, $c_i = \sum_{j=1}^{i-2} b_i + a_{i-1}$
  // from those coefficients. The main advantage of this kind of scheme is the
  // fact that only two vectors are needed per stage, namely the accumulated
  // part of the solution $\mathbf{w}$ (that will hold the solution
  // $\mathbf{w}^{n+1}$ at the new time $t^{n+1}$ after the last stage), the
  // update vector $\mathbf{r}_i$ that gets evaluated during the stages, plus
  // one vector $\mathbf{k}_i$ to hold the evaluation of the operator. Such a
  // Runge--Kutta setup reduces the memory storage and memory access. As the
  // memory bandwidth is often the performance-limiting factor on modern
  // hardware when the evaluation of the differential operator is
  // well-optimized, performance can be improved over standard time
  // integrators. This is true also when taking into account that a
  // conventional Runge--Kutta scheme might allow for slightly larger time
  // steps as more free parameters allow for better stability properties.
  //
  // In this tutorial programs, we concentrate on a few variants of
  // low-storage schemes defined in the article by Kennedy, Carpenter, and
  // Lewis (2000), as well as one variant described by Tselios and Simos
  // (2007). There is a large series of other schemes available, which could
  // be addressed by additional sets of coefficients or slightly different
  // update formulas.
  //
  // We define a single class for the four integrators, distinguished by the
  // enum described above. To each scheme, we then fill the vectors for the
  // $b_i$ and $a_i$ to the given variables in the class.
  class LowStorageRungeKuttaIntegrator
  {
  public:
    LowStorageRungeKuttaIntegrator(const LowStorageRungeKuttaScheme scheme)
    {
      TimeStepping::runge_kutta_method lsrk;
      // First comes the three-stage scheme of order three by Kennedy et al.
      // (2000). While its stability region is significantly smaller than for
      // the other schemes, it only involves three stages, so it is very
      // competitive in terms of the work per stage.
      switch (scheme)
        {
          case stage_3_order_3:
            {
              lsrk = TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3;
              break;
            }

            // The next scheme is a five-stage scheme of order four, again
            // defined in the paper by Kennedy et al. (2000).
          case stage_5_order_4:
            {
              lsrk = TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4;
              break;
            }

            // The following scheme of seven stages and order four has been
            // explicitly derived for acoustics problems. It is a balance of
            // accuracy for imaginary eigenvalues among fourth order schemes,
            // combined with a large stability region. Since DG schemes are
            // dissipative among the highest frequencies, this does not
            // necessarily translate to the highest possible time step per
            // stage. In the context of the present tutorial program, the
            // numerical flux plays a crucial role in the dissipation and thus
            // also the maximal stable time step size. For the modified
            // Lax--Friedrichs flux, this scheme is similar to the
            // `stage_5_order_4` scheme in terms of step size per stage if only
            // stability is considered, but somewhat less efficient for the HLL
            // flux.
          case stage_7_order_4:
            {
              lsrk = TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4;
              break;
            }

            // The last scheme included here is the nine-stage scheme of order
            // five from Kennedy et al. (2000). It is the most accurate among
            // the schemes used here, but the higher order of accuracy
            // sacrifices some stability, so the step length normalized per
            // stage is less than for the fourth order schemes.
          case stage_9_order_5:
            {
              lsrk = TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5;
              break;
            }

          default:
            AssertThrow(false, ExcNotImplemented());
        }
      TimeStepping::LowStorageRungeKutta<
        LinearAlgebra::distributed::Vector<Number>>
        rk_integrator(lsrk);
      rk_integrator.get_coefficients(ai, bi, ci);
    }

    unsigned int n_stages() const
    {
      return bi.size();
    }

    // The main function of the time integrator is to go through the stages,
    // evaluate the operator, prepare the $\mathbf{r}_i$ vector for the next
    // evaluation, and update the solution vector $\mathbf{w}$. We hand off
    // the work to the `pde_operator` involved in order to be able to merge
    // the vector operations of the Runge--Kutta setup with the evaluation of
    // the differential operator for better performance, so all we do here is
    // to delegate the vectors and coefficients.
    //
    // We separately call the operator for the first stage because we need
    // slightly modified arguments there: We evaluate the solution from
    // the old solution $\mathbf{w}^n$ rather than a $\mathbf r_i$ vector, so
    // the first argument is `solution`. We here let the stage vector
    // $\mathbf{r}_i$ also hold the temporary result of the evaluation, as it
    // is not used otherwise. For all subsequent stages, we use the vector
    // `vec_ki` as the second vector argument to store the result of the
    // operator evaluation. Finally, when we are at the last stage, we must
    // skip the computation of the vector $\mathbf{r}_{s+1}$ as there is no
    // coefficient $a_s$ available (nor will it be used).
    template <typename VectorType, typename Operator>
    void perform_time_step(const Operator &pde_operator,
                           const double    current_time,
                           const double    time_step,
                           VectorType &    solution,
                           VectorType &    vec_ri,
                           VectorType &    vec_ki) const
    {
      AssertDimension(ai.size() + 1, bi.size());

      pde_operator.perform_stage(current_time,
                                 bi[0] * time_step,
                                 ai[0] * time_step,
                                 solution,
                                 vec_ri,
                                 solution,
                                 vec_ri);

      for (unsigned int stage = 1; stage < bi.size(); ++stage)
        {
          const double c_i = ci[stage];
          pde_operator.perform_stage(current_time + c_i * time_step,
                                     bi[stage] * time_step,
                                     (stage == bi.size() - 1 ?
                                        0 :
                                        ai[stage] * time_step),
                                     vec_ri,
                                     vec_ki,
                                     solution,
                                     vec_ri);
        }
    }

  private:
    std::vector<double> bi;
    std::vector<double> ai;
    std::vector<double> ci;
  };



  // @sect3{Implementation of point-wise operations of the Euler equations}

  // In the following functions, we implement the various problem-specific
  // operators pertaining to the Euler equations. Each function acts on the
  // vector of conserved variables $[\rho, \rho\mathbf{u}, E]$ that we hold in
  // the solution vectors, and computes various derived quantities.
  //
  // First out is the computation of the velocity, that we derive from the
  // momentum variable $\rho \mathbf{u}$ by division by $\rho$. One thing to
  // note here is that we decorate all those functions with the keyword
  // `DEAL_II_ALWAYS_INLINE`. This is a special macro that maps to a
  // compiler-specific keyword that tells the compiler to never create a
  // function call for any of those functions, and instead move the
  // implementation <a
  // href="https://en.wikipedia.org/wiki/Inline_function">inline</a> to where
  // they are called. This is critical for performance because we call into some
  // of those functions millions or billions of times: For example, we both use
  // the velocity for the computation of the flux further down, but also for the
  // computation of the pressure, and both of these places are evaluated at
  // every quadrature point of every cell. Making sure these functions are
  // inlined ensures not only that the processor does not have to execute a jump
  // instruction into the function (and the corresponding return jump), but also
  // that the compiler can re-use intermediate information from one function's
  // context in code that comes after the place where the function was called.
  // (We note that compilers are generally quite good at figuring out which
  // functions to inline by themselves. Here is a place where compilers may or
  // may not have figured it out by themselves but where we know for sure that
  // inlining is a win.)
  //
  // Another trick we apply is a separate variable for the inverse density
  // $\frac{1}{\rho}$. This enables the compiler to only perform a single
  // division for the flux, despite the division being used at several
  // places. As divisions are around ten to twenty times as expensive as
  // multiplications or additions, avoiding redundant divisions is crucial for
  // performance. We note that taking the inverse first and later multiplying
  // with it is not equivalent to a division in floating point arithmetic due
  // to roundoff effects, so the compiler is not allowed to exchange one way by
  // the other with standard optimization flags. However, it is also not
  // particularly difficult to write the code in the right way.
  //
  // To summarize, the chosen strategy of always inlining and careful
  // definition of expensive arithmetic operations allows us to write compact
  // code without passing all intermediate results around, despite making sure
  // that the code maps to excellent machine code.
  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim, Number>
    euler_velocity(const Tensor<1, dim + 2, Number> &conserved_variables)
  {
    const Number inverse_density = Number(1.) / conserved_variables[0];

    Tensor<1, dim, Number> velocity;
    for (unsigned int d = 0; d < dim; ++d)
      velocity[d] = conserved_variables[1 + d] * inverse_density;

    return velocity;
  }

  // The next function computes the pressure from the vector of conserved
  // variables, using the formula $p = (\gamma - 1) \left(E - \frac 12 \rho
  // \mathbf{u}\cdot \mathbf{u}\right)$. As explained above, we use the
  // velocity from the `euler_velocity()` function. Note that we need to
  // specify the first template argument `dim` here because the compiler is
  // not able to deduce it from the arguments of the tensor, whereas the
  // second argument (number type) can be automatically deduced.
  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Number
    euler_pressure(const Tensor<1, dim + 2, Number> &conserved_variables)
  {
    const Tensor<1, dim, Number> velocity =
      euler_velocity<dim>(conserved_variables);

    Number rho_u_dot_u = conserved_variables[1] * velocity[0];
    for (unsigned int d = 1; d < dim; ++d)
      rho_u_dot_u += conserved_variables[1 + d] * velocity[d];

    return (gamma - 1.) * (conserved_variables[dim + 1] - 0.5 * rho_u_dot_u);
  }

  // Here is the definition of the Euler flux function, i.e., the definition
  // of the actual equation. Given the velocity and pressure (that the
  // compiler optimization will make sure are done only once), this is
  // straight-forward given the equation stated in the introduction.
  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim + 2, Tensor<1, dim, Number>>
    euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables)
  {
    const Tensor<1, dim, Number> velocity =
      euler_velocity<dim>(conserved_variables);
    const Number pressure = euler_pressure<dim>(conserved_variables);

    Tensor<1, dim + 2, Tensor<1, dim, Number>> flux;
    for (unsigned int d = 0; d < dim; ++d)
      {
        flux[0][d] = conserved_variables[1 + d];
        for (unsigned int e = 0; e < dim; ++e)
          flux[e + 1][d] = conserved_variables[e + 1] * velocity[d];
        flux[d + 1][d] += pressure;
        flux[dim + 1][d] =
          velocity[d] * (conserved_variables[dim + 1] + pressure);
      }

    return flux;
  }

  // This next function is a helper to simplify the implementation of the
  // numerical flux, implementing the action of a tensor of tensors (with
  // non-standard outer dimension of size `dim + 2`, so the standard overloads
  // provided by deal.II's tensor classes do not apply here) with another
  // tensor of the same inner dimension, i.e., a matrix-vector product.
  template <int n_components, int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, n_components, Number>
    operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix,
              const Tensor<1, dim, Number> &                         vector)
  {
    Tensor<1, n_components, Number> result;
    for (unsigned int d = 0; d < n_components; ++d)
      result[d] = matrix[d] * vector;
    return result;
  }

  // This function implements the numerical flux (Riemann solver). It gets the
  // state from the two sides of an interface and the normal vector, oriented
  // from the side of the solution $\mathbf{w}^-$ towards the solution
  // $\mathbf{w}^+$. In finite volume methods which rely on piece-wise
  // constant data, the numerical flux is the central ingredient as it is the
  // only place where the physical information is entered. In DG methods, the
  // numerical flux is less central due to the polynomials within the elements
  // and the physical flux used there. As a result of higher-degree
  // interpolation with consistent values from both sides in the limit of a
  // continuous solution, the numerical flux can be seen as a control of the
  // jump of the solution from both sides to weakly impose continuity. It is
  // important to realize that a numerical flux alone cannot stabilize a
  // high-order DG method in the presence of shocks, and thus any DG method
  // must be combined with further shock-capturing techniques to handle those
  // cases. In this tutorial, we focus on wave-like solutions of the Euler
  // equations in the subsonic regime without strong discontinuities where our
  // basic scheme is sufficient.
  //
  // Nonetheless, the numerical flux is decisive in terms of the numerical
  // dissipation of the overall scheme and influences the admissible time step
  // size with explicit Runge--Kutta methods. We consider two choices, a
  // modified Lax--Friedrichs scheme and the widely used Harten--Lax--van Leer
  // (HLL) flux. For both variants, we first need to get the velocities and
  // pressures from both sides of the interface and evaluate the physical
  // Euler flux.
  //
  // For the local Lax--Friedrichs flux, the definition is $\hat{\mathbf{F}}
  // =\frac{\mathbf{F}(\mathbf{w}^-)+\mathbf{F}(\mathbf{w}^+)}{2} +
  // \frac{\lambda}{2}\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes
  // \mathbf{n^-}$, where the factor $\lambda =
  // \max\left(\|\mathbf{u}^-\|+c^-, \|\mathbf{u}^+\|+c^+\right)$ gives the
  // maximal wave speed and $c = \sqrt{\gamma p / \rho}$ is the speed of
  // sound. Here, we choose two modifications of that expression for reasons
  // of computational efficiency, given the small impact of the flux on the
  // solution. For the above definition of the factor $\lambda$, we would need
  // to take four square roots, two for the two velocity norms and two for the
  // speed of sound on either side. The first modification is hence to rather
  // use $\sqrt{\|\mathbf{u}\|^2+c^2}$ as an estimate of the maximal speed
  // (which is at most a factor of 2 away from the actual maximum, as shown in
  // the introduction). This allows us to pull the square root out of the
  // maximum and get away with a single square root computation. The second
  // modification is to further relax on the parameter $\lambda$---the smaller
  // it is, the smaller the dissipation factor (which is multiplied by the
  // jump in $\mathbf{w}$, which might result in a smaller or bigger
  // dissipation in the end). This allows us to fit the spectrum into the
  // stability region of the explicit Runge--Kutta integrator with bigger time
  // steps. However, we cannot make dissipation too small because otherwise
  // imaginary eigenvalues grow larger. Finally, the current conservative
  // formulation is not energy-stable in the limit of $\lambda\to 0$ as it is
  // not skew-symmetric, and would need additional measures such as split-form
  // DG schemes in that case.
  //
  // For the HLL flux, we follow the formula from literature, introducing an
  // additional weighting of the two states from Lax--Friedrichs by a
  // parameter $s$. It is derived from the physical transport directions of
  // the Euler equations in terms of the current direction of velocity and
  // sound speed. For the velocity, we here choose a simple arithmetic average
  // which is sufficient for DG scenarios and moderate jumps in material
  // parameters.
  //
  // Since the numerical flux is multiplied by the normal vector in the weak
  // form, we multiply by the result by the normal vector for all terms in the
  // equation. In these multiplications, the `operator*` defined above enables
  // a compact notation similar to the mathematical definition.
  //
  // In this and the following functions, we use variable suffixes `_m` and
  // `_p` to indicate quantities derived from $\mathbf{w}^-$ and $\mathbf{w}^+$,
  // i.e., values "here" and "there" relative to the current cell when looking
  // at a neighbor cell.
  template <int dim, typename Number>
  inline DEAL_II_ALWAYS_INLINE //
    Tensor<1, dim + 2, Number>
    euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m,
                         const Tensor<1, dim + 2, Number> &u_p,
                         const Tensor<1, dim, Number> &    normal)
  {
    const auto velocity_m = euler_velocity<dim>(u_m);
    const auto velocity_p = euler_velocity<dim>(u_p);

    const auto pressure_m = euler_pressure<dim>(u_m);
    const auto pressure_p = euler_pressure<dim>(u_p);

    const auto flux_m = euler_flux<dim>(u_m);
    const auto flux_p = euler_flux<dim>(u_p);

    switch (numerical_flux_type)
      {
        case lax_friedrichs_modified:
          {
            const auto lambda =
              0.5 * std::sqrt(std::max(velocity_p.norm_square() +
                                         gamma * pressure_p * (1. / u_p[0]),
                                       velocity_m.norm_square() +
                                         gamma * pressure_m * (1. / u_m[0])));

            return 0.5 * (flux_m * normal + flux_p * normal) +
                   0.5 * lambda * (u_m - u_p);
          }

        case harten_lax_vanleer:
          {
            const auto avg_velocity_normal =
              0.5 * ((velocity_m + velocity_p) * normal);
            const auto   avg_c = std::sqrt(std::abs(
              0.5 * gamma *
              (pressure_p * (1. / u_p[0]) + pressure_m * (1. / u_m[0]))));
            const Number s_pos =
              std::max(Number(), avg_velocity_normal + avg_c);
            const Number s_neg =
              std::min(Number(), avg_velocity_normal - avg_c);
            const Number inverse_s = Number(1.) / (s_pos - s_neg);

            return inverse_s *
                   ((s_pos * (flux_m * normal) - s_neg * (flux_p * normal)) -
                    s_pos * s_neg * (u_m - u_p));
          }

        default:
          {
            Assert(false, ExcNotImplemented());
            return {};
          }
      }
  }



  // This and the next function are helper functions to provide compact
  // evaluation calls as multiple points get batched together via a
  // VectorizedArray argument (see the step-37 tutorial for details). This
  // function is used for the subsonic outflow boundary conditions where we
  // need to set the energy component to a prescribed value. The next one
  // requests the solution on all components and is used for inflow boundaries
  // where all components of the solution are set.
  template <int dim, typename Number>
  VectorizedArray<Number>
  evaluate_function(const Function<dim> &                      function,
                    const Point<dim, VectorizedArray<Number>> &p_vectorized,
                    const unsigned int                         component)
  {
    VectorizedArray<Number> result;
    for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
      {
        Point<dim> p;
        for (unsigned int d = 0; d < dim; ++d)
          p[d] = p_vectorized[d][v];
        result[v] = function.value(p, component);
      }
    return result;
  }


  template <int dim, typename Number, int n_components = dim + 2>
  Tensor<1, n_components, VectorizedArray<Number>>
  evaluate_function(const Function<dim> &                      function,
                    const Point<dim, VectorizedArray<Number>> &p_vectorized)
  {
    AssertDimension(function.n_components, n_components);
    Tensor<1, n_components, VectorizedArray<Number>> result;
    for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
      {
        Point<dim> p;
        for (unsigned int d = 0; d < dim; ++d)
          p[d] = p_vectorized[d][v];
        for (unsigned int d = 0; d < n_components; ++d)
          result[d][v] = function.value(p, d);
      }
    return result;
  }



  // @sect3{The EulerOperation class}

  // This class implements the evaluators for the Euler problem, in analogy to
  // the `LaplaceOperator` class of step-37 or step-59. Since the present
  // operator is non-linear and does not require a matrix interface (to be
  // handed over to preconditioners), we skip the various `vmult` functions
  // otherwise present in matrix-free operators and only implement an `apply`
  // function as well as the combination of `apply` with the required vector
  // updates for the low-storage Runge--Kutta time integrator mentioned above
  // (called `perform_stage`). Furthermore, we have added three additional
  // functions involving matrix-free routines, namely one to compute an
  // estimate of the time step scaling (that is combined with the Courant
  // number for the actual time step size) based on the velocity and speed of
  // sound in the elements, one for the projection of solutions (specializing
  // VectorTools::project() for the DG case), and one to compute the errors
  // against a possible analytical solution or norms against some background
  // state.
  //
  // The rest of the class is similar to other matrix-free tutorials. As
  // discussed in the introduction, we provide a few functions to allow a user
  // to pass in various forms of boundary conditions on different parts of the
  // domain boundary marked by types::boundary_id variables, as well as
  // possible body forces.
  template <int dim, int degree, int n_points_1d>
  class EulerOperator
  {
  public:
    static constexpr unsigned int n_quadrature_points_1d = n_points_1d;

    EulerOperator(TimerOutput &timer_output);

    void reinit(const Mapping<dim> &   mapping,
                const DoFHandler<dim> &dof_handler);

    void set_inflow_boundary(const types::boundary_id       boundary_id,
                             std::unique_ptr<Function<dim>> inflow_function);

    void set_subsonic_outflow_boundary(
      const types::boundary_id       boundary_id,
      std::unique_ptr<Function<dim>> outflow_energy);

    void set_wall_boundary(const types::boundary_id boundary_id);

    void set_body_force(std::unique_ptr<Function<dim>> body_force);

    void apply(const double                                      current_time,
               const LinearAlgebra::distributed::Vector<Number> &src,
               LinearAlgebra::distributed::Vector<Number> &      dst) const;

    void
    perform_stage(const Number cur_time,
                  const Number factor_solution,
                  const Number factor_ai,
                  const LinearAlgebra::distributed::Vector<Number> &current_ri,
                  LinearAlgebra::distributed::Vector<Number> &      vec_ki,
                  LinearAlgebra::distributed::Vector<Number> &      solution,
                  LinearAlgebra::distributed::Vector<Number> &next_ri) const;

    void project(const Function<dim> &                       function,
                 LinearAlgebra::distributed::Vector<Number> &solution) const;

    std::array<double, 3> compute_errors(
      const Function<dim> &                             function,
      const LinearAlgebra::distributed::Vector<Number> &solution) const;

    double compute_cell_transport_speed(
      const LinearAlgebra::distributed::Vector<Number> &solution) const;

    void
    initialize_vector(LinearAlgebra::distributed::Vector<Number> &vector) const;

  private:
    MatrixFree<dim, Number> data;

    TimerOutput &timer;

    std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
      inflow_boundaries;
    std::map<types::boundary_id, std::unique_ptr<Function<dim>>>
                                   subsonic_outflow_boundaries;
    std::set<types::boundary_id>   wall_boundaries;
    std::unique_ptr<Function<dim>> body_force;

    void local_apply_inverse_mass_matrix(
      const MatrixFree<dim, Number> &                   data,
      LinearAlgebra::distributed::Vector<Number> &      dst,
      const LinearAlgebra::distributed::Vector<Number> &src,
      const std::pair<unsigned int, unsigned int> &     cell_range) const;

    void local_apply_cell(
      const MatrixFree<dim, Number> &                   data,
      LinearAlgebra::distributed::Vector<Number> &      dst,
      const LinearAlgebra::distributed::Vector<Number> &src,
      const std::pair<unsigned int, unsigned int> &     cell_range) const;

    void local_apply_face(
      const MatrixFree<dim, Number> &                   data,
      LinearAlgebra::distributed::Vector<Number> &      dst,
      const LinearAlgebra::distributed::Vector<Number> &src,
      const std::pair<unsigned int, unsigned int> &     face_range) const;

    void local_apply_boundary_face(
      const MatrixFree<dim, Number> &                   data,
      LinearAlgebra::distributed::Vector<Number> &      dst,
      const LinearAlgebra::distributed::Vector<Number> &src,
      const std::pair<unsigned int, unsigned int> &     face_range) const;
  };



  template <int dim, int degree, int n_points_1d>
  EulerOperator<dim, degree, n_points_1d>::EulerOperator(TimerOutput &timer)
    : timer(timer)
  {}



  // For the initialization of the Euler operator, we set up the MatrixFree
  // variable contained in the class. This can be done given a mapping to
  // describe possible curved boundaries as well as a DoFHandler object
  // describing the degrees of freedom. Since we use a discontinuous Galerkin
  // discretization in this tutorial program where no constraints are imposed
  // strongly on the solution field, we do not need to pass in an
  // AffineConstraints object and rather use a dummy for the
  // construction. With respect to quadrature, we want to select two different
  // ways of computing the underlying integrals: The first is a flexible one,
  // based on a template parameter `n_points_1d` (that will be assigned the
  // `n_q_points_1d` value specified at the top of this file). More accurate
  // integration is necessary to avoid the aliasing problem due to the
  // variable coefficients in the Euler operator. The second less accurate
  // quadrature formula is a tight one based on `fe_degree+1` and needed for
  // the inverse mass matrix. While that formula provides an exact inverse
  // only on affine element shapes and not on deformed elements, it enables
  // the fast inversion of the mass matrix by tensor product techniques,
  // necessary to ensure optimal computational efficiency overall.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::reinit(
    const Mapping<dim> &   mapping,
    const DoFHandler<dim> &dof_handler)
  {
    const std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_handler};
    const AffineConstraints<double>            dummy;
    const std::vector<const AffineConstraints<double> *> constraints = {&dummy};
    const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d),
                                                    QGauss<1>(fe_degree + 1)};

    typename MatrixFree<dim, Number>::AdditionalData additional_data;
    additional_data.mapping_update_flags =
      (update_gradients | update_JxW_values | update_quadrature_points |
       update_values);
    additional_data.mapping_update_flags_inner_faces =
      (update_JxW_values | update_quadrature_points | update_normal_vectors |
       update_values);
    additional_data.mapping_update_flags_boundary_faces =
      (update_JxW_values | update_quadrature_points | update_normal_vectors |
       update_values);
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim, Number>::AdditionalData::none;

    data.reinit(
      mapping, dof_handlers, constraints, quadratures, additional_data);
  }



  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::initialize_vector(
    LinearAlgebra::distributed::Vector<Number> &vector) const
  {
    data.initialize_dof_vector(vector);
  }



  // The subsequent four member functions are the ones that must be called from
  // outside to specify the various types of boundaries. For an inflow boundary,
  // we must specify all components in terms of density $\rho$, momentum $\rho
  // \mathbf{u}$ and energy $E$. Given this information, we then store the
  // function alongside the respective boundary id in a map member variable of
  // this class. Likewise, we proceed for the subsonic outflow boundaries (where
  // we request a function as well, which we use to retrieve the energy) and for
  // wall (no-penetration) boundaries where we impose zero normal velocity (no
  // function necessary, so we only request the boundary id). For the present
  // DG code where boundary conditions are solely applied as part of the weak
  // form (during time integration), the call to set the boundary conditions
  // can appear both before or after the `reinit()` call to this class. This
  // is different from continuous finite element codes where the boundary
  // conditions determine the content of the AffineConstraints object that is
  // sent into MatrixFree for initialization, thus requiring to be set before
  // the initialization of the matrix-free data structures.
  //
  // The checks added in each of the four function are used to
  // ensure that boundary conditions are mutually exclusive on the various
  // parts of the boundary, i.e., that a user does not accidentally designate a
  // boundary as both an inflow and say a subsonic outflow boundary.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::set_inflow_boundary(
    const types::boundary_id       boundary_id,
    std::unique_ptr<Function<dim>> inflow_function)
  {
    AssertThrow(subsonic_outflow_boundaries.find(boundary_id) ==
                    subsonic_outflow_boundaries.end() &&
                  wall_boundaries.find(boundary_id) == wall_boundaries.end(),
                ExcMessage("You already set the boundary with id " +
                           std::to_string(static_cast<int>(boundary_id)) +
                           " to another type of boundary before now setting " +
                           "it as inflow"));
    AssertThrow(inflow_function->n_components == dim + 2,
                ExcMessage("Expected function with dim+2 components"));

    inflow_boundaries[boundary_id] = std::move(inflow_function);
  }


  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::set_subsonic_outflow_boundary(
    const types::boundary_id       boundary_id,
    std::unique_ptr<Function<dim>> outflow_function)
  {
    AssertThrow(inflow_boundaries.find(boundary_id) ==
                    inflow_boundaries.end() &&
                  wall_boundaries.find(boundary_id) == wall_boundaries.end(),
                ExcMessage("You already set the boundary with id " +
                           std::to_string(static_cast<int>(boundary_id)) +
                           " to another type of boundary before now setting " +
                           "it as subsonic outflow"));
    AssertThrow(outflow_function->n_components == dim + 2,
                ExcMessage("Expected function with dim+2 components"));

    subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function);
  }


  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::set_wall_boundary(
    const types::boundary_id boundary_id)
  {
    AssertThrow(inflow_boundaries.find(boundary_id) ==
                    inflow_boundaries.end() &&
                  subsonic_outflow_boundaries.find(boundary_id) ==
                    subsonic_outflow_boundaries.end(),
                ExcMessage("You already set the boundary with id " +
                           std::to_string(static_cast<int>(boundary_id)) +
                           " to another type of boundary before now setting " +
                           "it as wall boundary"));

    wall_boundaries.insert(boundary_id);
  }


  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::set_body_force(
    std::unique_ptr<Function<dim>> body_force)
  {
    AssertDimension(body_force->n_components, dim);

    this->body_force = std::move(body_force);
  }



  // @sect4{Local evaluators}

  // Now we proceed to the local evaluators for the Euler problem. The
  // evaluators are relatively simple and follow what has been presented in
  // step-37, step-48, or step-59. The first notable difference is the fact
  // that we use an FEEvaluation with a non-standard number of quadrature
  // points. Whereas we previously always set the number of quadrature points
  // to equal the polynomial degree plus one (ensuring exact integration on
  // affine element shapes), we now set the number quadrature points as a
  // separate variable (e.g. the polynomial degree plus two or three halves of
  // the polynomial degree) to more accurately handle nonlinear terms. Since
  // the evaluator is fed with the appropriate loop lengths via the template
  // argument and keeps the number of quadrature points in the whole cell in
  // the variable FEEvaluation::n_q_points, we now automatically operate on
  // the more accurate formula without further changes.
  //
  // The second difference is due to the fact that we are now evaluating a
  // multi-component system, as opposed to the scalar systems considered
  // previously. The matrix-free framework provides several ways to handle the
  // multi-component case. The variant shown here utilizes an FEEvaluation
  // object with multiple components embedded into it, specified by the fourth
  // template argument `dim + 2` for the components in the Euler system. As a
  // consequence, the return type of FEEvaluation::get_value() is not a scalar
  // any more (that would return a VectorizedArray type, collecting data from
  // several elements), but a Tensor of `dim+2` components. The functionality
  // is otherwise similar to the scalar case; it is handled by a template
  // specialization of a base class, called FEEvaluationAccess. An alternative
  // variant would have been to use several FEEvaluation objects, a scalar one
  // for the density, a vector-valued one with `dim` components for the
  // momentum, and another scalar evaluator for the energy. To ensure that
  // those components point to the correct part of the solution, the
  // constructor of FEEvaluation takes three optional integer arguments after
  // the required MatrixFree field, namely the number of the DoFHandler for
  // multi-DoFHandler systems (taking the first by default), the number of the
  // quadrature point in case there are multiple Quadrature objects (see more
  // below), and as a third argument the component within a vector system. As
  // we have a single vector for all components, we would go with the third
  // argument, and set it to `0` for the density, `1` for the vector-valued
  // momentum, and `dim+1` for the energy slot. FEEvaluation then picks the
  // appropriate subrange of the solution vector during
  // FEEvaluationBase::read_dof_values() and
  // FEEvaluation::distributed_local_to_global() or the more compact
  // FEEvaluation::gather_evaluate() and FEEvaluation::integrate_scatter()
  // calls.
  //
  // When it comes to the evaluation of the body force vector, we distinguish
  // between two cases for efficiency reasons: In case we have a constant
  // function (derived from Functions::ConstantFunction), we can precompute
  // the value outside the loop over quadrature points and simply use the
  // value everywhere. For a more general function, we instead need to call
  // the `evaluate_function()` method we provided above; this path is more
  // expensive because we need to access the memory associated with the
  // quadrature point data.
  //
  // The rest follows the other tutorial programs. Since we have implemented
  // all physics for the Euler equations in the separate `euler_flux()`
  // function, all we have to do here is to call this function
  // given the current solution evaluated at quadrature points, returned by
  // `phi.get_value(q)`, and tell the FEEvaluation object to queue the flux
  // for testing it by the gradients of the shape functions (which is a Tensor
  // of outer `dim+2` components, each holding a tensor of `dim` components
  // for the $x,y,z$ component of the Euler flux). One final thing worth
  // mentioning is the order in which we queue the data for testing by the
  // value of the test function, `phi.submit_value()`, in case we are given an
  // external function: We must do this after calling `phi.get_value(q)`,
  // because `get_value()` (reading the solution) and `submit_value()`
  // (queuing the value for multiplication by the test function and summation
  // over quadrature points) access the same underlying data field. Here it
  // would be easy to achieve also without temporary variable `w_q` since
  // there is no mixing between values and gradients. For more complicated
  // setups, one has to first copy out e.g. both the value and gradient at a
  // quadrature point and then queue results again by
  // FEEvaluationBase::submit_value() and FEEvaluationBase::submit_gradient().
  //
  // As a final note, we mention that we do not use the first MatrixFree
  // argument of this function, which is a call-back from MatrixFree::loop().
  // The interfaces imposes the present list of arguments, but since we are in
  // a member function where the MatrixFree object is already available as the
  // `data` variable, we stick with that to avoid confusion.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::local_apply_cell(
    const MatrixFree<dim, Number> &,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src,
    const std::pair<unsigned int, unsigned int> &     cell_range) const
  {
    FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data);

    Tensor<1, dim, VectorizedArray<Number>> constant_body_force;
    const Functions::ConstantFunction<dim> *constant_function =
      dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get());

    if (constant_function)
      constant_body_force = evaluate_function<dim, Number, dim>(
        *constant_function, Point<dim, VectorizedArray<Number>>());

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, EvaluationFlags::values);

        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            const auto w_q = phi.get_value(q);
            phi.submit_gradient(euler_flux<dim>(w_q), q);
            if (body_force.get() != nullptr)
              {
                const Tensor<1, dim, VectorizedArray<Number>> force =
                  constant_function ? constant_body_force :
                                      evaluate_function<dim, Number, dim>(
                                        *body_force, phi.quadrature_point(q));

                Tensor<1, dim + 2, VectorizedArray<Number>> forcing;
                for (unsigned int d = 0; d < dim; ++d)
                  forcing[d + 1] = w_q[0] * force[d];
                for (unsigned int d = 0; d < dim; ++d)
                  forcing[dim + 1] += force[d] * w_q[d + 1];

                phi.submit_value(forcing, q);
              }
          }

        phi.integrate_scatter(((body_force.get() != nullptr) ?
                                 EvaluationFlags::values :
                                 EvaluationFlags::nothing) |
                                EvaluationFlags::gradients,
                              dst);
      }
  }



  // The next function concerns the computation of integrals on interior
  // faces, where we need evaluators from both cells adjacent to the face. We
  // associate the variable `phi_m` with the solution component $\mathbf{w}^-$
  // and the variable `phi_p` with the solution component $\mathbf{w}^+$. We
  // distinguish the two sides in the constructor of FEFaceEvaluation by the
  // second argument, with `true` for the interior side and `false` for the
  // exterior side, with interior and exterior denoting the orientation with
  // respect to the normal vector.
  //
  // Note that the calls FEFaceEvaluation::gather_evaluate() and
  // FEFaceEvaluation::integrate_scatter() combine the access to the vectors
  // and the sum factorization parts. This combined operation not only saves a
  // line of code, but also contains an important optimization: Given that we
  // use a nodal basis in terms of the Lagrange polynomials in the points of
  // the Gauss-Lobatto quadrature formula, only $(p+1)^{d-1}$ out of the
  // $(p+1)^d$ basis functions evaluate to non-zero on each face. Thus, the
  // evaluator only accesses the necessary data in the vector and skips the
  // parts which are multiplied by zero. If we had first read the vector, we
  // would have needed to load all data from the vector, as the call in
  // isolation would not know what data is required in subsequent
  // operations. If the subsequent FEFaceEvaluation::evaluate() call requests
  // values and derivatives, indeed all $(p+1)^d$ vector entries for each
  // component are needed, as the normal derivative is nonzero for all basis
  // functions.
  //
  // The arguments to the evaluators as well as the procedure is similar to
  // the cell evaluation. We again use the more accurate (over-)integration
  // scheme due to the nonlinear terms, specified as the third template
  // argument in the list. At the quadrature points, we then go to our
  // free-standing function for the numerical flux. It receives the solution
  // evaluated at quadrature points from both sides (i.e., $\mathbf{w}^-$ and
  // $\mathbf{w}^+$), as well as the normal vector onto the minus side. As
  // explained above, the numerical flux is already multiplied by the normal
  // vector from the minus side. We need to switch the sign because the
  // boundary term comes with a minus sign in the weak form derived in the
  // introduction. The flux is then queued for testing both on the minus sign
  // and on the plus sign, with switched sign as the normal vector from the
  // plus side is exactly opposed to the one from the minus side.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::local_apply_face(
    const MatrixFree<dim, Number> &,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src,
    const std::pair<unsigned int, unsigned int> &     face_range) const
  {
    FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_m(data,
                                                                      true);
    FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_p(data,
                                                                      false);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_p.reinit(face);
        phi_p.gather_evaluate(src, EvaluationFlags::values);

        phi_m.reinit(face);
        phi_m.gather_evaluate(src, EvaluationFlags::values);

        for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
          {
            const auto numerical_flux =
              euler_numerical_flux<dim>(phi_m.get_value(q),
                                        phi_p.get_value(q),
                                        phi_m.get_normal_vector(q));
            phi_m.submit_value(-numerical_flux, q);
            phi_p.submit_value(numerical_flux, q);
          }

        phi_p.integrate_scatter(EvaluationFlags::values, dst);
        phi_m.integrate_scatter(EvaluationFlags::values, dst);
      }
  }



  // For faces located at the boundary, we need to impose the appropriate
  // boundary conditions. In this tutorial program, we implement four cases as
  // mentioned above. (A fifth case, for supersonic outflow conditions is
  // discussed in the "Results" section below.) The discontinuous Galerkin
  // method imposes boundary conditions not as constraints, but only
  // weakly. Thus, the various conditions are imposed by finding an appropriate
  // <i>exterior</i> quantity $\mathbf{w}^+$ that is then handed to the
  // numerical flux function also used for the interior faces. In essence,
  // we "pretend" a state on the outside of the domain in such a way that
  // if that were reality, the solution of the PDE would satisfy the boundary
  // conditions we want.
  //
  // For wall boundaries, we need to impose a no-normal-flux condition on the
  // momentum variable, whereas we use a Neumann condition for the density and
  // energy with $\rho^+ = \rho^-$ and $E^+ = E^-$. To achieve the no-normal
  // flux condition, we set the exterior values to the interior values and
  // subtract two times the velocity in wall-normal direction, i.e., in the
  // direction of the normal vector.
  //
  // For inflow boundaries, we simply set the given Dirichlet data
  // $\mathbf{w}_\mathrm{D}$ as a boundary value. An alternative would have been
  // to use $\mathbf{w}^+ = -\mathbf{w}^- + 2 \mathbf{w}_\mathrm{D}$, the
  // so-called mirror principle.
  //
  // The imposition of outflow is essentially a Neumann condition, i.e.,
  // setting $\mathbf{w}^+ = \mathbf{w}^-$. For the case of subsonic outflow,
  // we still need to impose a value for the energy, which we derive from the
  // respective function. A special step is needed for the case of
  // <i>backflow</i>, i.e., the case where there is a momentum flux into the
  // domain on the Neumann portion. According to the literature (a fact that can
  // be derived by appropriate energy arguments), we must switch to another
  // variant of the flux on inflow parts, see Gravemeier, Comerford,
  // Yoshihara, Ismail, Wall, "A novel formulation for Neumann inflow
  // conditions in biomechanics", Int. J. Numer. Meth. Biomed. Eng., vol. 28
  // (2012). Here, the momentum term needs to be added once again, which
  // corresponds to removing the flux contribution on the momentum
  // variables. We do this in a post-processing step, and only for the case
  // when we both are at an outflow boundary and the dot product between the
  // normal vector and the momentum (or, equivalently, velocity) is
  // negative. As we work on data of several quadrature points at once for
  // SIMD vectorizations, we here need to explicitly loop over the array
  // entries of the SIMD array.
  //
  // In the implementation below, we check for the various types
  // of boundaries at the level of quadrature points. Of course, we could also
  // have moved the decision out of the quadrature point loop and treat entire
  // faces as of the same kind, which avoids some map/set lookups in the inner
  // loop over quadrature points. However, the loss of efficiency is hardly
  // noticeable, so we opt for the simpler code here. Also note that the final
  // `else` clause will catch the case when some part of the boundary was not
  // assigned any boundary condition via `EulerOperator::set_..._boundary(...)`.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::local_apply_boundary_face(
    const MatrixFree<dim, Number> &,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src,
    const std::pair<unsigned int, unsigned int> &     face_range) const
  {
    FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, true);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi.reinit(face);
        phi.gather_evaluate(src, EvaluationFlags::values);

        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            const auto w_m    = phi.get_value(q);
            const auto normal = phi.get_normal_vector(q);

            auto rho_u_dot_n = w_m[1] * normal[0];
            for (unsigned int d = 1; d < dim; ++d)
              rho_u_dot_n += w_m[1 + d] * normal[d];

            bool at_outflow = false;

            Tensor<1, dim + 2, VectorizedArray<Number>> w_p;
            const auto boundary_id = data.get_boundary_id(face);
            if (wall_boundaries.find(boundary_id) != wall_boundaries.end())
              {
                w_p[0] = w_m[0];
                for (unsigned int d = 0; d < dim; ++d)
                  w_p[d + 1] = w_m[d + 1] - 2. * rho_u_dot_n * normal[d];
                w_p[dim + 1] = w_m[dim + 1];
              }
            else if (inflow_boundaries.find(boundary_id) !=
                     inflow_boundaries.end())
              w_p =
                evaluate_function(*inflow_boundaries.find(boundary_id)->second,
                                  phi.quadrature_point(q));
            else if (subsonic_outflow_boundaries.find(boundary_id) !=
                     subsonic_outflow_boundaries.end())
              {
                w_p          = w_m;
                w_p[dim + 1] = evaluate_function(
                  *subsonic_outflow_boundaries.find(boundary_id)->second,
                  phi.quadrature_point(q),
                  dim + 1);
                at_outflow = true;
              }
            else
              AssertThrow(false,
                          ExcMessage("Unknown boundary id, did "
                                     "you set a boundary condition for "
                                     "this part of the domain boundary?"));

            auto flux = euler_numerical_flux<dim>(w_m, w_p, normal);

            if (at_outflow)
              for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
                {
                  if (rho_u_dot_n[v] < -1e-12)
                    for (unsigned int d = 0; d < dim; ++d)
                      flux[d + 1][v] = 0.;
                }

            phi.submit_value(-flux, q);
          }

        phi.integrate_scatter(EvaluationFlags::values, dst);
      }
  }



  // The next function implements the inverse mass matrix operation. The
  // algorithms and rationale have been discussed extensively in the
  // introduction, so we here limit ourselves to the technicalities of the
  // MatrixFreeOperators::CellwiseInverseMassMatrix class. It does similar
  // operations as the forward evaluation of the mass matrix, except with a
  // different interpolation matrix, representing the inverse $S^{-1}$
  // factors. These represent a change of basis from the specified basis (in
  // this case, the Lagrange basis in the points of the Gauss--Lobatto
  // quadrature formula) to the Lagrange basis in the points of the Gauss
  // quadrature formula. In the latter basis, we can apply the inverse of the
  // point-wise `JxW` factor, i.e., the quadrature weight times the
  // determinant of the Jacobian of the mapping from reference to real
  // coordinates. Once this is done, the basis is changed back to the nodal
  // Gauss-Lobatto basis again. All of these operations are done by the
  // `apply()` function below. What we need to provide is the local fields to
  // operate on (which we extract from the global vector by an FEEvaluation
  // object) and write the results back to the destination vector of the mass
  // matrix operation.
  //
  // One thing to note is that we added two integer arguments (that are
  // optional) to the constructor of FEEvaluation, the first being 0
  // (selecting among the DoFHandler in multi-DoFHandler systems; here, we
  // only have one) and the second being 1 to make the quadrature formula
  // selection. As we use the quadrature formula 0 for the over-integration of
  // nonlinear terms, we use the formula 1 with the default $p+1$ (or
  // `fe_degree+1` in terms of the variable name) points for the mass
  // matrix. This leads to square contributions to the mass matrix and ensures
  // exact integration, as explained in the introduction.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::local_apply_inverse_mass_matrix(
    const MatrixFree<dim, Number> &,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src,
    const std::pair<unsigned int, unsigned int> &     cell_range) const
  {
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
      inverse(phi);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);

        inverse.apply(phi.begin_dof_values(), phi.begin_dof_values());

        phi.set_dof_values(dst);
      }
  }



  // @sect4{The apply() and related functions}

  // We now come to the function which implements the evaluation of the Euler
  // operator as a whole, i.e., $\mathcal M^{-1} \mathcal L(t, \mathbf{w})$,
  // calling into the local evaluators presented above. The steps should be
  // clear from the previous code. One thing to note is that we need to adjust
  // the time in the functions we have associated with the various parts of
  // the boundary, in order to be consistent with the equation in case the
  // boundary data is time-dependent. Then, we call MatrixFree::loop() to
  // perform the cell and face integrals, including the necessary ghost data
  // exchange in the `src` vector. The seventh argument to the function,
  // `true`, specifies that we want to zero the `dst` vector as part of the
  // loop, before we start accumulating integrals into it. This variant is
  // preferred over explicitly calling `dst = 0.;` before the loop as the
  // zeroing operation is done on a subrange of the vector in parts that are
  // written by the integrals nearby. This enhances data locality and allows
  // for caching, saving one roundtrip of vector data to main memory and
  // enhancing performance. The last two arguments to the loop determine which
  // data is exchanged: Since we only access the values of the shape functions
  // one faces, typical of first-order hyperbolic problems, and since we have
  // a nodal basis with nodes at the reference element surface, we only need
  // to exchange those parts. This again saves precious memory bandwidth.
  //
  // Once the spatial operator $\mathcal L$ is applied, we need to make a
  // second round and apply the inverse mass matrix. Here, we call
  // MatrixFree::cell_loop() since only cell integrals appear. The cell loop
  // is cheaper than the full loop as access only goes to the degrees of
  // freedom associated with the locally owned cells, which is simply the
  // locally owned degrees of freedom for DG discretizations. Thus, no ghost
  // exchange is needed here.
  //
  // Around all these functions, we put timer scopes to record the
  // computational time for statistics about the contributions of the various
  // parts.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::apply(
    const double                                      current_time,
    const LinearAlgebra::distributed::Vector<Number> &src,
    LinearAlgebra::distributed::Vector<Number> &      dst) const
  {
    {
      TimerOutput::Scope t(timer, "apply - integrals");

      for (auto &i : inflow_boundaries)
        i.second->set_time(current_time);
      for (auto &i : subsonic_outflow_boundaries)
        i.second->set_time(current_time);

      data.loop(&EulerOperator::local_apply_cell,
                &EulerOperator::local_apply_face,
                &EulerOperator::local_apply_boundary_face,
                this,
                dst,
                src,
                true,
                MatrixFree<dim, Number>::DataAccessOnFaces::values,
                MatrixFree<dim, Number>::DataAccessOnFaces::values);
    }

    {
      TimerOutput::Scope t(timer, "apply - inverse mass");

      data.cell_loop(&EulerOperator::local_apply_inverse_mass_matrix,
                     this,
                     dst,
                     dst);
    }
  }



  // Let us move to the function that does an entire stage of a Runge--Kutta
  // update. It calls EulerOperator::apply() followed by some updates
  // to the vectors, namely `next_ri = solution + factor_ai * k_i` and
  // `solution += factor_solution * k_i`. Rather than performing these
  // steps through the vector interfaces, we here present an alternative
  // strategy that is faster on cache-based architectures. As the memory
  // consumed by the vectors is often much larger than what fits into caches,
  // the data has to effectively come from the slow RAM memory. The situation
  // can be improved by loop fusion, i.e., performing both the updates to
  // `next_ki` and `solution` within a single sweep. In that case, we would
  // read the two vectors `rhs` and `solution` and write into `next_ki` and
  // `solution`, compared to at least 4 reads and two writes in the baseline
  // case. Here, we go one step further and perform the loop immediately when
  // the mass matrix inversion has finished on a part of the
  // vector. MatrixFree::cell_loop() provides a mechanism to attach an
  // `std::function` both before the loop over cells first touches a vector
  // entry (which we do not use here, but is e.g. used for zeroing the vector)
  // and a second `std::function` to be called after the loop last touches
  // an entry. The callback is in form of a range over the given vector (in
  // terms of the local index numbering in the MPI universe) that can be
  // addressed by `local_element()` functions.
  //
  // For this second callback, we create a lambda that works on a range and
  // write the respective update on this range. Ideally, we would add the
  // `DEAL_II_OPENMP_SIMD_PRAGMA` before the local loop to suggest to the
  // compiler to SIMD parallelize this loop (which means in practice that we
  // ensure that there is no overlap, also called aliasing, between the index
  // ranges of the pointers we use inside the loops). It turns out that at the
  // time of this writing, GCC 7.2 fails to compile an OpenMP pragma inside a
  // lambda function, so we comment this pragma out below. If your compiler is
  // newer, you should be able to uncomment these lines again.
  //
  // Note that we select a different code path for the last
  // Runge--Kutta stage when we do not need to update the `next_ri`
  // vector. This strategy gives a considerable speedup. Whereas the inverse
  // mass matrix and vector updates take more than 60% of the computational
  // time with default vector updates on a 40-core machine, the percentage is
  // around 35% with the more optimized variant. In other words, this is a
  // speedup of around a third.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::perform_stage(
    const Number                                      current_time,
    const Number                                      factor_solution,
    const Number                                      factor_ai,
    const LinearAlgebra::distributed::Vector<Number> &current_ri,
    LinearAlgebra::distributed::Vector<Number> &      vec_ki,
    LinearAlgebra::distributed::Vector<Number> &      solution,
    LinearAlgebra::distributed::Vector<Number> &      next_ri) const
  {
    {
      TimerOutput::Scope t(timer, "rk_stage - integrals L_h");

      for (auto &i : inflow_boundaries)
        i.second->set_time(current_time);
      for (auto &i : subsonic_outflow_boundaries)
        i.second->set_time(current_time);

      data.loop(&EulerOperator::local_apply_cell,
                &EulerOperator::local_apply_face,
                &EulerOperator::local_apply_boundary_face,
                this,
                vec_ki,
                current_ri,
                true,
                MatrixFree<dim, Number>::DataAccessOnFaces::values,
                MatrixFree<dim, Number>::DataAccessOnFaces::values);
    }


    {
      TimerOutput::Scope t(timer, "rk_stage - inv mass + vec upd");
      data.cell_loop(
        &EulerOperator::local_apply_inverse_mass_matrix,
        this,
        next_ri,
        vec_ki,
        std::function<void(const unsigned int, const unsigned int)>(),
        [&](const unsigned int start_range, const unsigned int end_range) {
          const Number ai = factor_ai;
          const Number bi = factor_solution;
          if (ai == Number())
            {
              /* DEAL_II_OPENMP_SIMD_PRAGMA */
              for (unsigned int i = start_range; i < end_range; ++i)
                {
                  const Number k_i          = next_ri.local_element(i);
                  const Number sol_i        = solution.local_element(i);
                  solution.local_element(i) = sol_i + bi * k_i;
                }
            }
          else
            {
              /* DEAL_II_OPENMP_SIMD_PRAGMA */
              for (unsigned int i = start_range; i < end_range; ++i)
                {
                  const Number k_i          = next_ri.local_element(i);
                  const Number sol_i        = solution.local_element(i);
                  solution.local_element(i) = sol_i + bi * k_i;
                  next_ri.local_element(i)  = sol_i + ai * k_i;
                }
            }
        });
    }
  }



  // Having discussed the implementation of the functions that deal with
  // advancing the solution by one time step, let us now move to functions
  // that implement other, ancillary operations. Specifically, these are
  // functions that compute projections, evaluate errors, and compute the speed
  // of information transport on a cell.
  //
  // The first of these functions is essentially equivalent to
  // VectorTools::project(), just much faster because it is specialized for DG
  // elements where there is no need to set up and solve a linear system, as
  // each element has independent basis functions. The reason why we show the
  // code here, besides a small speedup of this non-critical operation, is that
  // it shows additional functionality provided by
  // MatrixFreeOperators::CellwiseInverseMassMatrix.
  //
  // The projection operation works as follows: If we denote the matrix of
  // shape functions evaluated at quadrature points by $S$, the projection on
  // cell $K$ is an operation of the form $\underbrace{S J^K S^\mathrm
  // T}_{\mathcal M^K} \mathbf{w}^K = S J^K
  // \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$, where $J^K$ is the diagonal
  // matrix containing the determinant of the Jacobian times the quadrature
  // weight (JxW), $\mathcal M^K$ is the cell-wise mass matrix, and
  // $\tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$ is the evaluation of the
  // field to be projected onto quadrature points. (In reality the matrix $S$
  // has additional structure through the tensor product, as explained in the
  // introduction.) This system can now equivalently be written as
  // $\mathbf{w}^K = \left(S J^K S^\mathrm T\right)^{-1} S J^K
  // \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q} = S^{-\mathrm T}
  // \left(J^K\right)^{-1} S^{-1} S J^K
  // \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$. Now, the term $S^{-1} S$ and
  // then $\left(J^K\right)^{-1} J^K$ cancel, resulting in the final
  // expression $\mathbf{w}^K = S^{-\mathrm T}
  // \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$. This operation is
  // implemented by
  // MatrixFreeOperators::CellwiseInverseMassMatrix::transform_from_q_points_to_basis().
  // The name is derived from the fact that this projection is simply
  // the multiplication by $S^{-\mathrm T}$, a basis change from the
  // nodal basis in the points of the Gaussian quadrature to the given finite
  // element basis. Note that we call FEEvaluation::set_dof_values() to write
  // the result into the vector, overwriting previous content, rather than
  // accumulating the results as typical in integration tasks -- we can do
  // this because every vector entry has contributions from only a single
  // cell for discontinuous Galerkin discretizations.
  template <int dim, int degree, int n_points_1d>
  void EulerOperator<dim, degree, n_points_1d>::project(
    const Function<dim> &                       function,
    LinearAlgebra::distributed::Vector<Number> &solution) const
  {
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number>
      inverse(phi);
    solution.zero_out_ghost_values();
    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        phi.reinit(cell);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_dof_value(evaluate_function(function,
                                                 phi.quadrature_point(q)),
                               q);
        inverse.transform_from_q_points_to_basis(dim + 2,
                                                 phi.begin_dof_values(),
                                                 phi.begin_dof_values());
        phi.set_dof_values(solution);
      }
  }



  // The next function again repeats functionality also provided by the
  // deal.II library, namely VectorTools::integrate_difference(). We here show
  // the explicit code to highlight how the vectorization across several cells
  // works and how to accumulate results via that interface: Recall that each
  // <i>lane</i> of the vectorized array holds data from a different cell. By
  // the loop over all cell batches that are owned by the current MPI process,
  // we could then fill a VectorizedArray of results; to obtain a global sum,
  // we would need to further go on and sum across the entries in the SIMD
  // array. However, such a procedure is not stable as the SIMD array could in
  // fact not hold valid data for all its lanes. This happens when the number
  // of locally owned cells is not a multiple of the SIMD width. To avoid
  // invalid data, we must explicitly skip those invalid lanes when accessing
  // the data. While one could imagine that we could make it work by simply
  // setting the empty lanes to zero (and thus, not contribute to a sum), the
  // situation is more complicated than that: What if we were to compute a
  // velocity out of the momentum? Then, we would need to divide by the
  // density, which is zero -- the result would consequently be NaN and
  // contaminate the result. This trap is avoided by accumulating the results
  // from the valid SIMD range as we loop through the cell batches, using the
  // function MatrixFree::n_active_entries_per_cell_batch() to give us the
  // number of lanes with valid data. It equals VectorizedArray::size() on
  // most cells, but can be less on the last cell batch if the number of cells
  // has a remainder compared to the SIMD width.
  template <int dim, int degree, int n_points_1d>
  std::array<double, 3> EulerOperator<dim, degree, n_points_1d>::compute_errors(
    const Function<dim> &                             function,
    const LinearAlgebra::distributed::Vector<Number> &solution) const
  {
    TimerOutput::Scope t(timer, "compute errors");
    double             errors_squared[3] = {};
    FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, 0, 0);

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(solution, EvaluationFlags::values);
        VectorizedArray<Number> local_errors_squared[3] = {};
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            const auto error =
              evaluate_function(function, phi.quadrature_point(q)) -
              phi.get_value(q);
            const auto JxW = phi.JxW(q);

            local_errors_squared[0] += error[0] * error[0] * JxW;
            for (unsigned int d = 0; d < dim; ++d)
              local_errors_squared[1] += (error[d + 1] * error[d + 1]) * JxW;
            local_errors_squared[2] += (error[dim + 1] * error[dim + 1]) * JxW;
          }
        for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
             ++v)
          for (unsigned int d = 0; d < 3; ++d)
            errors_squared[d] += local_errors_squared[d][v];
      }

    Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared);

    std::array<double, 3> errors;
    for (unsigned int d = 0; d < 3; ++d)
      errors[d] = std::sqrt(errors_squared[d]);

    return errors;
  }



  // This final function of the EulerOperator class is used to estimate the
  // transport speed, scaled by the mesh size, that is relevant for setting
  // the time step size in the explicit time integrator. In the Euler
  // equations, there are two speeds of transport, namely the convective
  // velocity $\mathbf{u}$ and the propagation of sound waves with sound
  // speed $c = \sqrt{\gamma p/\rho}$ relative to the medium moving at
  // velocity $\mathbf u$.
  //
  // In the formula for the time step size, we are interested not by
  // these absolute speeds, but by the amount of time it takes for
  // information to cross a single cell. For information transported along with
  // the medium, $\mathbf u$ is scaled by the mesh size,
  // so an estimate of the maximal velocity can be obtained by computing
  // $\|J^{-\mathrm T} \mathbf{u}\|_\infty$, where $J$ is the Jacobian of the
  // transformation from real to the reference domain. Note that
  // FEEvaluationBase::inverse_jacobian() returns the inverse and transpose
  // Jacobian, representing the metric term from real to reference
  // coordinates, so we do not need to transpose it again. We store this limit
  // in the variable `convective_limit` in the code below.
  //
  // The sound propagation is isotropic, so we need to take mesh sizes in any
  // direction into account. The appropriate mesh size scaling is then given
  // by the minimal singular value of $J$ or, equivalently, the maximal
  // singular value of $J^{-1}$. Note that one could approximate this quantity
  // by the minimal distance between vertices of a cell when ignoring curved
  // cells. To get the maximal singular value of the Jacobian, the general
  // strategy would be some LAPACK function. Since all we need here is an
  // estimate, we can avoid the hassle of decomposing a tensor of
  // VectorizedArray numbers into several matrices and go into an (expensive)
  // eigenvalue function without vectorization, and instead use a few
  // iterations (five in the code below) of the power method applied to
  // $J^{-1}J^{-\mathrm T}$. The speed of convergence of this method depends
  // on the ratio of the largest to the next largest eigenvalue and the
  // initial guess, which is the vector of all ones. This might suggest that
  // we get slow convergence on cells close to a cube shape where all
  // lengths are almost the same. However, this slow convergence means that
  // the result will sit between the two largest singular values, which both
  // are close to the maximal value anyway. In all other cases, convergence
  // will be quick. Thus, we can merely hardcode 5 iterations here and be
  // confident that the result is good.
  template <int dim, int degree, int n_points_1d>
  double EulerOperator<dim, degree, n_points_1d>::compute_cell_transport_speed(
    const LinearAlgebra::distributed::Vector<Number> &solution) const
  {
    TimerOutput::Scope t(timer, "compute transport speed");
    Number             max_transport = 0;
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1);

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(solution, EvaluationFlags::values);
        VectorizedArray<Number> local_max = 0.;
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            const auto solution = phi.get_value(q);
            const auto velocity = euler_velocity<dim>(solution);
            const auto pressure = euler_pressure<dim>(solution);

            const auto inverse_jacobian = phi.inverse_jacobian(q);
            const auto convective_speed = inverse_jacobian * velocity;
            VectorizedArray<Number> convective_limit = 0.;
            for (unsigned int d = 0; d < dim; ++d)
              convective_limit =
                std::max(convective_limit, std::abs(convective_speed[d]));

            const auto speed_of_sound =
              std::sqrt(gamma * pressure * (1. / solution[0]));

            Tensor<1, dim, VectorizedArray<Number>> eigenvector;
            for (unsigned int d = 0; d < dim; ++d)
              eigenvector[d] = 1.;
            for (unsigned int i = 0; i < 5; ++i)
              {
                eigenvector = transpose(inverse_jacobian) *
                              (inverse_jacobian * eigenvector);
                VectorizedArray<Number> eigenvector_norm = 0.;
                for (unsigned int d = 0; d < dim; ++d)
                  eigenvector_norm =
                    std::max(eigenvector_norm, std::abs(eigenvector[d]));
                eigenvector /= eigenvector_norm;
              }
            const auto jac_times_ev   = inverse_jacobian * eigenvector;
            const auto max_eigenvalue = std::sqrt(
              (jac_times_ev * jac_times_ev) / (eigenvector * eigenvector));
            local_max =
              std::max(local_max,
                       max_eigenvalue * speed_of_sound + convective_limit);
          }

        // Similarly to the previous function, we must make sure to accumulate
        // speed only on the valid cells of a cell batch.
        for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
             ++v)
          for (unsigned int d = 0; d < 3; ++d)
            max_transport = std::max(max_transport, local_max[v]);
      }

    max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD);

    return max_transport;
  }



  // @sect3{The EulerProblem class}

  // This class combines the EulerOperator class with the time integrator and
  // the usual global data structures such as FiniteElement and DoFHandler, to
  // actually run the simulations of the Euler problem.
  //
  // The member variables are a triangulation, a finite element, a mapping (to
  // create high-order curved surfaces, see e.g. step-10), and a DoFHandler to
  // describe the degrees of freedom. In addition, we keep an instance of the
  // EulerOperator described above around, which will do all heavy lifting in
  // terms of integrals, and some parameters for time integration like the
  // current time or the time step size.
  //
  // Furthermore, we use a PostProcessor instance to write some additional
  // information to the output file, in similarity to what was done in
  // step-33. The interface of the DataPostprocessor class is intuitive,
  // requiring us to provide information about what needs to be evaluated
  // (typically only the values of the solution, except for the Schlieren plot
  // that we only enable in 2D where it makes sense), and the names of what
  // gets evaluated. Note that it would also be possible to extract most
  // information by calculator tools within visualization programs such as
  // ParaView, but it is so much more convenient to do it already when writing
  // the output.
  template <int dim>
  class EulerProblem
  {
  public:
    EulerProblem();

    void run();

  private:
    void make_grid_and_dofs();

    void output_results(const unsigned int result_number);

    LinearAlgebra::distributed::Vector<Number> solution;

    ConditionalOStream pcout;

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim> triangulation;
#else
    Triangulation<dim> triangulation;
#endif

    FESystem<dim>        fe;
    MappingQGeneric<dim> mapping;
    DoFHandler<dim>      dof_handler;

    TimerOutput timer;

    EulerOperator<dim, fe_degree, n_q_points_1d> euler_operator;

    double time, time_step;

    class Postprocessor : public DataPostprocessor<dim>
    {
    public:
      Postprocessor();

      virtual void evaluate_vector_field(
        const DataPostprocessorInputs::Vector<dim> &inputs,
        std::vector<Vector<double>> &computed_quantities) const override;

      virtual std::vector<std::string> get_names() const override;

      virtual std::vector<
        DataComponentInterpretation::DataComponentInterpretation>
      get_data_component_interpretation() const override;

      virtual UpdateFlags get_needed_update_flags() const override;

    private:
      const bool do_schlieren_plot;
    };
  };



  template <int dim>
  EulerProblem<dim>::Postprocessor::Postprocessor()
    : do_schlieren_plot(dim == 2)
  {}



  // For the main evaluation of the field variables, we first check that the
  // lengths of the arrays equal the expected values (the lengths `2*dim+4` or
  // `2*dim+5` are derived from the sizes of the names we specify in the
  // get_names() function below). Then we loop over all evaluation points and
  // fill the respective information: First we fill the primal solution
  // variables of density $\rho$, momentum $\rho \mathbf{u}$ and energy $E$,
  // then we compute the derived velocity $\mathbf u$, the pressure $p$, the
  // speed of sound $c=\sqrt{\gamma p / \rho}$, as well as the Schlieren plot
  // showing $s = |\nabla \rho|^2$ in case it is enabled. (See step-69 for
  // another example where we create a Schlieren plot.)
  template <int dim>
  void EulerProblem<dim>::Postprocessor::evaluate_vector_field(
    const DataPostprocessorInputs::Vector<dim> &inputs,
    std::vector<Vector<double>> &               computed_quantities) const
  {
    const unsigned int n_evaluation_points = inputs.solution_values.size();

    if (do_schlieren_plot == true)
      Assert(inputs.solution_gradients.size() == n_evaluation_points,
             ExcInternalError());

    Assert(computed_quantities.size() == n_evaluation_points,
           ExcInternalError());
    Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError());
    Assert(computed_quantities[0].size() ==
             dim + 2 + (do_schlieren_plot == true ? 1 : 0),
           ExcInternalError());

    for (unsigned int q = 0; q < n_evaluation_points; ++q)
      {
        Tensor<1, dim + 2> solution;
        for (unsigned int d = 0; d < dim + 2; ++d)
          solution[d] = inputs.solution_values[q](d);

        const double         density  = solution[0];
        const Tensor<1, dim> velocity = euler_velocity<dim>(solution);
        const double         pressure = euler_pressure<dim>(solution);

        for (unsigned int d = 0; d < dim; ++d)
          computed_quantities[q](d) = velocity[d];
        computed_quantities[q](dim)     = pressure;
        computed_quantities[q](dim + 1) = std::sqrt(gamma * pressure / density);

        if (do_schlieren_plot == true)
          computed_quantities[q](dim + 2) =
            inputs.solution_gradients[q][0] * inputs.solution_gradients[q][0];
      }
  }



  template <int dim>
  std::vector<std::string> EulerProblem<dim>::Postprocessor::get_names() const
  {
    std::vector<std::string> names;
    for (unsigned int d = 0; d < dim; ++d)
      names.emplace_back("velocity");
    names.emplace_back("pressure");
    names.emplace_back("speed_of_sound");

    if (do_schlieren_plot == true)
      names.emplace_back("schlieren_plot");

    return names;
  }



  // For the interpretation of quantities, we have scalar density, energy,
  // pressure, speed of sound, and the Schlieren plot, and vectors for the
  // momentum and the velocity.
  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  EulerProblem<dim>::Postprocessor::get_data_component_interpretation() const
  {
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation;
    for (unsigned int d = 0; d < dim; ++d)
      interpretation.push_back(
        DataComponentInterpretation::component_is_part_of_vector);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    if (do_schlieren_plot == true)
      interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);

    return interpretation;
  }



  // With respect to the necessary update flags, we only need the values for
  // all quantities but the Schlieren plot, which is based on the density
  // gradient.
  template <int dim>
  UpdateFlags EulerProblem<dim>::Postprocessor::get_needed_update_flags() const
  {
    if (do_schlieren_plot == true)
      return update_values | update_gradients;
    else
      return update_values;
  }



  // The constructor for this class is unsurprising: We set up a parallel
  // triangulation based on the `MPI_COMM_WORLD` communicator, a vector finite
  // element with `dim+2` components for density, momentum, and energy, a
  // high-order mapping of the same degree as the underlying finite element,
  // and initialize the time and time step to zero.
  template <int dim>
  EulerProblem<dim>::EulerProblem()
    : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
#ifdef DEAL_II_WITH_P4EST
    , triangulation(MPI_COMM_WORLD)
#endif
    , fe(FE_DGQ<dim>(fe_degree), dim + 2)
    , mapping(fe_degree)
    , dof_handler(triangulation)
    , timer(pcout, TimerOutput::never, TimerOutput::wall_times)
    , euler_operator(timer)
    , time(0)
    , time_step(0)
  {}



  // As a mesh, this tutorial program implements two options, depending on the
  // global variable `testcase`: For the analytical variant (`testcase==0`),
  // the domain is $(0, 10) \times (-5, 5)$, with Dirichlet boundary
  // conditions (inflow) all around the domain. For `testcase==1`, we set the
  // domain to a cylinder in a rectangular box, derived from the flow past
  // cylinder testcase for incompressible viscous flow by Sch&auml;fer and
  // Turek (1996). Here, we have a larger variety of boundaries. The inflow
  // part at the left of the channel is given the inflow type, for which we
  // choose a constant inflow profile, whereas we set a subsonic outflow at
  // the right. For the boundary around the cylinder (boundary id equal to 2)
  // as well as the channel walls (boundary id equal to 3) we use the wall
  // boundary type, which is no-normal flow. Furthermore, for the 3D cylinder
  // we also add a gravity force in vertical direction. Having the base mesh
  // in place (including the manifolds set by
  // GridGenerator::channel_with_cylinder()), we can then perform the
  // specified number of global refinements, create the unknown numbering from
  // the DoFHandler, and hand the DoFHandler and Mapping objects to the
  // initialization of the EulerOperator.
  template <int dim>
  void EulerProblem<dim>::make_grid_and_dofs()
  {
    switch (testcase)
      {
        case 0:
          {
            Point<dim> lower_left;
            for (unsigned int d = 1; d < dim; ++d)
              lower_left[d] = -5;

            Point<dim> upper_right;
            upper_right[0] = 10;
            for (unsigned int d = 1; d < dim; ++d)
              upper_right[d] = 5;

            GridGenerator::hyper_rectangle(triangulation,
                                           lower_left,
                                           upper_right);
            triangulation.refine_global(2);

            euler_operator.set_inflow_boundary(
              0, std::make_unique<ExactSolution<dim>>(0));

            break;
          }

        case 1:
          {
            GridGenerator::channel_with_cylinder(
              triangulation, 0.03, 1, 0, true);

            euler_operator.set_inflow_boundary(
              0, std::make_unique<ExactSolution<dim>>(0));
            euler_operator.set_subsonic_outflow_boundary(
              1, std::make_unique<ExactSolution<dim>>(0));

            euler_operator.set_wall_boundary(2);
            euler_operator.set_wall_boundary(3);

            if (dim == 3)
              euler_operator.set_body_force(
                std::make_unique<Functions::ConstantFunction<dim>>(
                  std::vector<double>({0., 0., -0.2})));

            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }

    triangulation.refine_global(n_global_refinements);

    dof_handler.distribute_dofs(fe);

    euler_operator.reinit(mapping, dof_handler);
    euler_operator.initialize_vector(solution);

    // In the following, we output some statistics about the problem. Because we
    // often end up with quite large numbers of cells or degrees of freedom, we
    // would like to print them with a comma to separate each set of three
    // digits. This can be done via "locales", although the way this works is
    // not particularly intuitive. step-32 explains this in slightly more
    // detail.
    std::locale s = pcout.get_stream().getloc();
    pcout.get_stream().imbue(std::locale(""));
    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << " ( = " << (dim + 2) << " [vars] x "
          << triangulation.n_global_active_cells() << " [cells] x "
          << Utilities::pow(fe_degree + 1, dim) << " [dofs/cell/var] )"
          << std::endl;
    pcout.get_stream().imbue(s);
  }



  // For output, we first let the Euler operator compute the errors of the
  // numerical results. More precisely, we compute the error against the
  // analytical result for the analytical solution case, whereas we compute
  // the deviation against the background field with constant density and
  // energy and constant velocity in $x$ direction for the second test case.
  //
  // The next step is to create output. This is similar to what is done in
  // step-33: We let the postprocessor defined above control most of the
  // output, except for the primal field that we write directly. For the
  // analytical solution test case, we also perform another projection of the
  // analytical solution and print the difference between that field and the
  // numerical solution. Once we have defined all quantities to be written, we
  // build the patches for output. Similarly to step-65, we create a
  // high-order VTK output by setting the appropriate flag, which enables us
  // to visualize fields of high polynomial degrees. Finally, we call the
  // `DataOutInterface::write_vtu_in_parallel()` function to write the result
  // to the given file name. This function uses special MPI parallel write
  // facilities, which are typically more optimized for parallel file systems
  // than the standard library's `std::ofstream` variants used in most other
  // tutorial programs. A particularly nice feature of the
  // `write_vtu_in_parallel()` function is the fact that it can combine output
  // from all MPI ranks into a single file, making it unnecessary to have a
  // central record of all such files (namely, the "pvtu" file).
  //
  // For parallel programs, it is often instructive to look at the partitioning
  // of cells among processors. To this end, one can pass a vector of numbers
  // to DataOut::add_data_vector() that contains as many entries as the
  // current processor has active cells; these numbers should then be the
  // rank of the processor that owns each of these cells. Such a vector
  // could, for example, be obtained from
  // GridTools::get_subdomain_association(). On the other hand, on each MPI
  // process, DataOut will only read those entries that correspond to locally
  // owned cells, and these of course all have the same value: namely, the rank
  // of the current process. What is in the remaining entries of the vector
  // doesn't actually matter, and so we can just get away with a cheap trick: We
  // just fill *all* values of the vector we give to DataOut::add_data_vector()
  // with the rank of the current MPI process. The key is that on each process,
  // only the entries corresponding to the locally owned cells will be read,
  // ignoring the (wrong) values in other entries. The fact that every process
  // submits a vector in which the correct subset of entries is correct is all
  // that is necessary.
  template <int dim>
  void EulerProblem<dim>::output_results(const unsigned int result_number)
  {
    const std::array<double, 3> errors =
      euler_operator.compute_errors(ExactSolution<dim>(time), solution);
    const std::string quantity_name = testcase == 0 ? "error" : "norm";

    pcout << "Time:" << std::setw(8) << std::setprecision(3) << time
          << ", dt: " << std::setw(8) << std::setprecision(2) << time_step
          << ", " << quantity_name << " rho: " << std::setprecision(4)
          << std::setw(10) << errors[0] << ", rho * u: " << std::setprecision(4)
          << std::setw(10) << errors[1] << ", energy:" << std::setprecision(4)
          << std::setw(10) << errors[2] << std::endl;

    {
      TimerOutput::Scope t(timer, "output");

      Postprocessor postprocessor;
      DataOut<dim>  data_out;

      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);

      data_out.attach_dof_handler(dof_handler);
      {
        std::vector<std::string> names;
        names.emplace_back("density");
        for (unsigned int d = 0; d < dim; ++d)
          names.emplace_back("momentum");
        names.emplace_back("energy");

        std::vector<DataComponentInterpretation::DataComponentInterpretation>
          interpretation;
        interpretation.push_back(
          DataComponentInterpretation::component_is_scalar);
        for (unsigned int d = 0; d < dim; ++d)
          interpretation.push_back(
            DataComponentInterpretation::component_is_part_of_vector);
        interpretation.push_back(
          DataComponentInterpretation::component_is_scalar);

        data_out.add_data_vector(dof_handler, solution, names, interpretation);
      }
      data_out.add_data_vector(solution, postprocessor);

      LinearAlgebra::distributed::Vector<Number> reference;
      if (testcase == 0 && dim == 2)
        {
          reference.reinit(solution);
          euler_operator.project(ExactSolution<dim>(time), reference);
          reference.sadd(-1., 1, solution);
          std::vector<std::string> names;
          names.emplace_back("error_density");
          for (unsigned int d = 0; d < dim; ++d)
            names.emplace_back("error_momentum");
          names.emplace_back("error_energy");

          std::vector<DataComponentInterpretation::DataComponentInterpretation>
            interpretation;
          interpretation.push_back(
            DataComponentInterpretation::component_is_scalar);
          for (unsigned int d = 0; d < dim; ++d)
            interpretation.push_back(
              DataComponentInterpretation::component_is_part_of_vector);
          interpretation.push_back(
            DataComponentInterpretation::component_is_scalar);

          data_out.add_data_vector(dof_handler,
                                   reference,
                                   names,
                                   interpretation);
        }

      Vector<double> mpi_owner(triangulation.n_active_cells());
      mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      data_out.add_data_vector(mpi_owner, "owner");

      data_out.build_patches(mapping,
                             fe.degree,
                             DataOut<dim>::curved_inner_cells);

      const std::string filename =
        "solution_" + Utilities::int_to_string(result_number, 3) + ".vtu";
      data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
    }
  }



  // The EulerProblem::run() function puts all pieces together. It starts off
  // by calling the function that creates the mesh and sets up data structures,
  // and then initializing the time integrator and the two temporary vectors of
  // the low-storage integrator. We call these vectors `rk_register_1` and
  // `rk_register_2`, and use the first vector to represent the quantity
  // $\mathbf{r}_i$ and the second one for $\mathbf{k}_i$ in the formulas for
  // the Runge--Kutta scheme outlined in the introduction. Before we start the
  // time loop, we compute the time step size by the
  // `EulerOperator::compute_cell_transport_speed()` function. For reasons of
  // comparison, we compare the result obtained there with the minimal mesh
  // size and print them to screen. For velocities and speeds of sound close
  // to unity as in this tutorial program, the predicted effective mesh size
  // will be close, but they could vary if scaling were different.
  template <int dim>
  void EulerProblem<dim>::run()
  {
    {
      const unsigned int n_vect_number = VectorizedArray<Number>::size();
      const unsigned int n_vect_bits   = 8 * sizeof(Number) * n_vect_number;

      pcout << "Running with "
            << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
            << " MPI processes" << std::endl;
      pcout << "Vectorization over " << n_vect_number << " "
            << (std::is_same<Number, double>::value ? "doubles" : "floats")
            << " = " << n_vect_bits << " bits ("
            << Utilities::System::get_current_vectorization_level() << ")"
            << std::endl;
    }

    make_grid_and_dofs();

    const LowStorageRungeKuttaIntegrator integrator(lsrk_scheme);

    LinearAlgebra::distributed::Vector<Number> rk_register_1;
    LinearAlgebra::distributed::Vector<Number> rk_register_2;
    rk_register_1.reinit(solution);
    rk_register_2.reinit(solution);

    euler_operator.project(ExactSolution<dim>(time), solution);

    double min_vertex_distance = std::numeric_limits<double>::max();
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        min_vertex_distance =
          std::min(min_vertex_distance, cell->minimum_vertex_distance());
    min_vertex_distance =
      Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD);

    time_step = courant_number * integrator.n_stages() /
                euler_operator.compute_cell_transport_speed(solution);
    pcout << "Time step size: " << time_step
          << ", minimal h: " << min_vertex_distance
          << ", initial transport scaling: "
          << 1. / euler_operator.compute_cell_transport_speed(solution)
          << std::endl
          << std::endl;

    output_results(0);

    // Now we are ready to start the time loop, which we run until the time
    // has reached the desired end time. Every 5 time steps, we compute a new
    // estimate for the time step -- since the solution is nonlinear, it is
    // most effective to adapt the value during the course of the
    // simulation. In case the Courant number was chosen too aggressively, the
    // simulation will typically blow up with time step NaN, so that is easy
    // to detect here. One thing to note is that roundoff errors might
    // propagate to the leading digits due to an interaction of slightly
    // different time step selections that in turn lead to slightly different
    // solutions. To decrease this sensitivity, it is common practice to round
    // or truncate the time step size to a few digits, e.g. 3 in this case. In
    // case the current time is near the prescribed 'tick' value for output
    // (e.g. 0.02), we also write the output. After the end of the time loop,
    // we summarize the computation by printing some statistics, which is
    // mostly done by the TimerOutput::print_wall_time_statistics() function.
    unsigned int timestep_number = 0;

    while (time < final_time - 1e-12)
      {
        ++timestep_number;
        if (timestep_number % 5 == 0)
          time_step =
            courant_number * integrator.n_stages() /
            Utilities::truncate_to_n_digits(
              euler_operator.compute_cell_transport_speed(solution), 3);

        {
          TimerOutput::Scope t(timer, "rk time stepping total");
          integrator.perform_time_step(euler_operator,
                                       time,
                                       time_step,
                                       solution,
                                       rk_register_1,
                                       rk_register_2);
        }

        time += time_step;

        if (static_cast<int>(time / output_tick) !=
              static_cast<int>((time - time_step) / output_tick) ||
            time >= final_time - 1e-12)
          output_results(
            static_cast<unsigned int>(std::round(time / output_tick)));
      }

    timer.print_wall_time_statistics(MPI_COMM_WORLD);
    pcout << std::endl;
  }

} // namespace Euler_DG



// The main() function is not surprising and follows what was done in all
// previous MPI programs: As we run an MPI program, we need to call `MPI_Init()`
// and `MPI_Finalize()`, which we do through the
// Utilities::MPI::MPI_InitFinalize data structure. Note that we run the program
// only with MPI, and set the thread count to 1.
int main(int argc, char **argv)
{
  using namespace Euler_DG;
  using namespace dealii;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  try
    {
      deallog.depth_console(0);

      EulerProblem<dimension> euler_problem;
      euler_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
