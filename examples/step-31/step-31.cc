/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2007 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Martin Kronbichler, Uppsala University,
 *          Wolfgang Bangerth, Texas A&M University 2007, 2008
 */


// @sect3{Include files}

// The first step, as always, is to include the functionality of these
// well-known deal.II library files and some C++ header files.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

// Then we need to include some header files that provide vector, matrix, and
// preconditioner classes that implement interfaces to the respective Trilinos
// classes. In particular, we will need interfaces to the matrix and vector
// classes based on Trilinos as well as Trilinos preconditioners:
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>

// Finally, here are two C++ headers that haven't been included yet by one of
// the aforelisted header files:
#include <fstream>
#include <sstream>
#include <limits>


// At the end of this top-matter, we import all deal.II names into the global
// namespace:
namespace Step31
{
  using namespace dealii;


  // @sect3{Equation data}

  // Again, the next stage in the program is the definition of the equation
  // data, that is, the various boundary conditions, the right hand sides and
  // the initial condition (remember that we're about to solve a
  // time-dependent system). The basic strategy for this definition is the
  // same as in step-22. Regarding the details, though, there are some
  // differences.

  // The first thing is that we don't set any nonhomogeneous boundary
  // conditions on the velocity, since as is explained in the introduction we
  // will use no-flux conditions $\mathbf{n}\cdot\mathbf{u}=0$. So what is
  // left are <code>dim-1</code> conditions for the tangential part of the
  // normal component of the stress tensor, $\textbf{n} \cdot [p \textbf{1} -
  // \eta\varepsilon(\textbf{u})]$; we assume homogeneous values for these
  // components, i.e. a natural boundary condition that requires no specific
  // action (it appears as a zero term in the right hand side of the weak
  // form).
  //
  // For the temperature <i>T</i>, we assume no thermal energy flux,
  // i.e. $\mathbf{n} \cdot \kappa \nabla T=0$. This, again, is a boundary
  // condition that does not require us to do anything in particular.
  //
  // Secondly, we have to set initial conditions for the temperature (no
  // initial conditions are required for the velocity and pressure, since the
  // Stokes equations for the quasi-stationary case we consider here have no
  // time derivatives of the velocity or pressure). Here, we choose a very
  // simple test case, where the initial temperature is zero, and all dynamics
  // are driven by the temperature right hand side.
  //
  // Thirdly, we need to define the right hand side of the temperature
  // equation. We choose it to be constant within three circles (or spheres in
  // 3d) somewhere at the bottom of the domain, as explained in the
  // introduction, and zero outside.
  //
  // Finally, or maybe firstly, at the top of this namespace, we define the
  // various material constants we need ($\eta,\kappa$, density $\rho$ and the
  // thermal expansion coefficient $\beta$):
  namespace EquationData
  {
    const double eta = 1;
    const double kappa = 1e-6;
    const double beta = 10;
    const double density = 1;


    template <int dim>
    class TemperatureInitialValues : public Function<dim>
    {
    public:
      TemperatureInitialValues () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
                                 Vector<double>   &value) const;
    };


    template <int dim>
    double
    TemperatureInitialValues<dim>::value (const Point<dim> &,
                                          const unsigned int) const
    {
      return 0;
    }


    template <int dim>
    void
    TemperatureInitialValues<dim>::vector_value (const Point<dim> &p,
                                                 Vector<double>   &values) const
    {
      for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = TemperatureInitialValues<dim>::value (p, c);
    }


    template <int dim>
    class TemperatureRightHandSide : public Function<dim>
    {
    public:
      TemperatureRightHandSide () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
                                 Vector<double>   &value) const;
    };


    template <int dim>
    double
    TemperatureRightHandSide<dim>::value (const Point<dim>  &p,
                                          const unsigned int component) const
    {
      Assert (component == 0,
              ExcMessage ("Invalid operation for a scalar function."));

      Assert ((dim==2) || (dim==3), ExcNotImplemented());

      static const Point<dim> source_centers[3]
        = { (dim == 2 ? Point<dim>(.3,.1) : Point<dim>(.3,.5,.1)),
            (dim == 2 ? Point<dim>(.45,.1) : Point<dim>(.45,.5,.1)),
            (dim == 2 ? Point<dim>(.75,.1) : Point<dim>(.75,.5,.1))
          };
      static const double source_radius
        = (dim == 2 ? 1./32 : 1./8);

      return ((source_centers[0].distance (p) < source_radius)
              ||
              (source_centers[1].distance (p) < source_radius)
              ||
              (source_centers[2].distance (p) < source_radius)
              ?
              1
              :
              0);
    }


    template <int dim>
    void
    TemperatureRightHandSide<dim>::vector_value (const Point<dim> &p,
                                                 Vector<double>   &values) const
    {
      for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = TemperatureRightHandSide<dim>::value (p, c);
    }
  }



  // @sect3{Linear solvers and preconditioners}

  // This section introduces some objects that are used for the solution of
  // the linear equations of the Stokes system that we need to solve in each
  // time step. Many of the ideas used here are the same as in step-20, where
  // Schur complement based preconditioners and solvers have been introduced,
  // with the actual interface taken from step-22 (in particular the
  // discussion in the "Results" section of step-22, in which we introduce
  // alternatives to the direct Schur complement approach). Note, however,
  // that here we don't use the Schur complement to solve the Stokes
  // equations, though an approximate Schur complement (the mass matrix on the
  // pressure space) appears in the preconditioner.
  namespace LinearSolvers
  {

    // @sect4{The <code>InverseMatrix</code> class template}

    // This class is an interface to calculate the action of an "inverted"
    // matrix on a vector (using the <code>vmult</code> operation) in the same
    // way as the corresponding class in step-22: when the product of an
    // object of this class is requested, we solve a linear equation system
    // with that matrix using the CG method, accelerated by a preconditioner
    // of (templated) class <code>Preconditioner</code>.
    //
    // In a minor deviation from the implementation of the same class in
    // step-22 (and step-20), we make the <code>vmult</code> function take any
    // kind of vector type (it will yield compiler errors, however, if the
    // matrix does not allow a matrix-vector product with this kind of
    // vector).
    //
    // Secondly, we catch any exceptions that the solver may have thrown. The
    // reason is as follows: When debugging a program like this one
    // occasionally makes a mistake of passing an indefinite or nonsymmetric
    // matrix or preconditioner to the current class. The solver will, in that
    // case, not converge and throw a run-time exception. If not caught here
    // it will propagate up the call stack and may end up in
    // <code>main()</code> where we output an error message that will say that
    // the CG solver failed. The question then becomes: Which CG solver? The
    // one that inverted the mass matrix? The one that inverted the top left
    // block with the Laplace operator? Or a CG solver in one of the several
    // other nested places where we use linear solvers in the current code? No
    // indication about this is present in a run-time exception because it
    // doesn't store the stack of calls through which we got to the place
    // where the exception was generated.
    //
    // So rather than letting the exception propagate freely up to
    // <code>main()</code> we realize that there is little that an outer
    // function can do if the inner solver fails and rather convert the
    // run-time exception into an assertion that fails and triggers a call to
    // <code>abort()</code>, allowing us to trace back in a debugger how we
    // got to the current place.
    template <class Matrix, class Preconditioner>
    class InverseMatrix : public Subscriptor
    {
    public:
      InverseMatrix (const Matrix         &m,
                     const Preconditioner &preconditioner);


      template <typename VectorType>
      void vmult (VectorType       &dst,
                  const VectorType &src) const;

    private:
      const SmartPointer<const Matrix> matrix;
      const Preconditioner &preconditioner;
    };


    template <class Matrix, class Preconditioner>
    InverseMatrix<Matrix,Preconditioner>::
    InverseMatrix (const Matrix &m,
                   const Preconditioner &preconditioner)
      :
      matrix (&m),
      preconditioner (preconditioner)
    {}



    template <class Matrix, class Preconditioner>
    template <typename VectorType>
    void
    InverseMatrix<Matrix,Preconditioner>::
    vmult (VectorType       &dst,
           const VectorType &src) const
    {
      SolverControl solver_control (src.size(), 1e-7*src.l2_norm());
      SolverCG<VectorType> cg (solver_control);

      dst = 0;

      try
        {
          cg.solve (*matrix, dst, src, preconditioner);
        }
      catch (std::exception &e)
        {
          Assert (false, ExcMessage(e.what()));
        }
    }

    // @sect4{Schur complement preconditioner}

    // This is the implementation of the Schur complement preconditioner as
    // described in detail in the introduction. As opposed to step-20 and
    // step-22, we solve the block system all-at-once using GMRES, and use the
    // Schur complement of the block structured matrix to build a good
    // preconditioner instead.
    //
    // Let's have a look at the ideal preconditioner matrix
    // $P=\left(\begin{array}{cc} A & 0 \\ B & -S \end{array}\right)$
    // described in the introduction. If we apply this matrix in the solution
    // of a linear system, convergence of an iterative GMRES solver will be
    // governed by the matrix @f{eqnarray*} P^{-1}\left(\begin{array}{cc} A &
    // B^T \\ B & 0 \end{array}\right) = \left(\begin{array}{cc} I & A^{-1}
    // B^T \\ 0 & I \end{array}\right), @f} which indeed is very simple. A
    // GMRES solver based on exact matrices would converge in one iteration,
    // since all eigenvalues are equal (any Krylov method takes at most as
    // many iterations as there are distinct eigenvalues). Such a
    // preconditioner for the blocked Stokes system has been proposed by
    // Silvester and Wathen ("Fast iterative solution of stabilised Stokes
    // systems part II.  Using general block preconditioners", SIAM
    // J. Numer. Anal., 31 (1994), pp. 1352-1367).
    //
    // Replacing <i>P</i> by $\tilde{P}$ keeps that spirit alive: the product
    // $P^{-1} A$ will still be close to a matrix with eigenvalues 1 with a
    // distribution that does not depend on the problem size. This lets us
    // hope to be able to get a number of GMRES iterations that is
    // problem-size independent.
    //
    // The deal.II users who have already gone through the step-20 and step-22
    // tutorials can certainly imagine how we're going to implement this.  We
    // replace the exact inverse matrices in $P^{-1}$ by some approximate
    // inverses built from the InverseMatrix class, and the inverse Schur
    // complement will be approximated by the pressure mass matrix $M_p$
    // (weighted by $\eta^{-1}$ as mentioned in the introduction). As pointed
    // out in the results section of step-22, we can replace the exact inverse
    // of <i>A</i> by just the application of a preconditioner, in this case
    // on a vector Laplace matrix as was explained in the introduction. This
    // does increase the number of (outer) GMRES iterations, but is still
    // significantly cheaper than an exact inverse, which would require
    // between 20 and 35 CG iterations for <em>each</em> outer solver step
    // (using the AMG preconditioner).
    //
    // Having the above explanations in mind, we define a preconditioner class
    // with a <code>vmult</code> functionality, which is all we need for the
    // interaction with the usual solver functions further below in the
    // program code.
    //
    // First the declarations. These are similar to the definition of the
    // Schur complement in step-20, with the difference that we need some more
    // preconditioners in the constructor and that the matrices we use here
    // are built upon Trilinos:
    template <class PreconditionerA, class PreconditionerMp>
    class BlockSchurPreconditioner : public Subscriptor
    {
    public:
      BlockSchurPreconditioner (
        const TrilinosWrappers::BlockSparseMatrix     &S,
        const InverseMatrix<TrilinosWrappers::SparseMatrix,
        PreconditionerMp>         &Mpinv,
        const PreconditionerA                         &Apreconditioner);

      void vmult (TrilinosWrappers::BlockVector       &dst,
                  const TrilinosWrappers::BlockVector &src) const;

    private:
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> stokes_matrix;
      const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix,
            PreconditionerMp > > m_inverse;
      const PreconditionerA &a_preconditioner;

      mutable TrilinosWrappers::Vector tmp;
    };



    template <class PreconditionerA, class PreconditionerMp>
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix  &S,
                             const InverseMatrix<TrilinosWrappers::SparseMatrix,
                             PreconditionerMp>      &Mpinv,
                             const PreconditionerA                      &Apreconditioner)
      :
      stokes_matrix           (&S),
      m_inverse               (&Mpinv),
      a_preconditioner        (Apreconditioner),
      tmp                     (stokes_matrix->block(1,1).m())
    {}


    // Next is the <code>vmult</code> function. We implement the action of
    // $P^{-1}$ as described above in three successive steps.  In formulas, we
    // want to compute $Y=P^{-1}X$ where $X,Y$ are both vectors with two block
    // components.
    //
    // The first step multiplies the velocity part of the vector by a
    // preconditioner of the matrix <i>A</i>, i.e. we compute $Y_0={\tilde
    // A}^{-1}X_0$.  The resulting velocity vector is then multiplied by $B$
    // and subtracted from the pressure, i.e. we want to compute $X_1-BY_0$.
    // This second step only acts on the pressure vector and is accomplished
    // by the residual function of our matrix classes, except that the sign is
    // wrong. Consequently, we change the sign in the temporary pressure
    // vector and finally multiply by the inverse pressure mass matrix to get
    // the final pressure vector, completing our work on the Stokes
    // preconditioner:
    template <class PreconditionerA, class PreconditionerMp>
    void
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    vmult (TrilinosWrappers::BlockVector       &dst,
           const TrilinosWrappers::BlockVector &src) const
    {
      a_preconditioner.vmult (dst.block(0), src.block(0));
      stokes_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
      tmp *= -1;
      m_inverse->vmult (dst.block(1), tmp);
    }
  }



  // @sect3{The <code>BoussinesqFlowProblem</code> class template}

  // The definition of the class that defines the top-level logic of solving
  // the time-dependent Boussinesq problem is mainly based on the step-22
  // tutorial program. The main differences are that now we also have to solve
  // for the temperature equation, which forces us to have a second DoFHandler
  // object for the temperature variable as well as matrices, right hand
  // sides, and solution vectors for the current and previous time steps. As
  // mentioned in the introduction, all linear algebra objects are going to
  // use wrappers of the corresponding Trilinos functionality.
  //
  // The member functions of this class are reminiscent of step-21, where we
  // also used a staggered scheme that first solve the flow equations (here
  // the Stokes equations, in step-21 Darcy flow) and then update the advected
  // quantity (here the temperature, there the saturation). The functions that
  // are new are mainly concerned with determining the time step, as well as
  // the proper size of the artificial viscosity stabilization.
  //
  // The last three variables indicate whether the various matrices or
  // preconditioners need to be rebuilt the next time the corresponding build
  // functions are called. This allows us to move the corresponding
  // <code>if</code> into the respective function and thereby keeping our main
  // <code>run()</code> function clean and easy to read.
  template <int dim>
  class BoussinesqFlowProblem
  {
  public:
    BoussinesqFlowProblem ();
    void run ();

  private:
    void setup_dofs ();
    void assemble_stokes_preconditioner ();
    void build_stokes_preconditioner ();
    void assemble_stokes_system ();
    void assemble_temperature_system (const double maximal_velocity);
    void assemble_temperature_matrix ();
    double get_maximal_velocity () const;
    std::pair<double,double> get_extrapolated_temperature_range () const;
    void solve ();
    void output_results () const;
    void refine_mesh (const unsigned int max_grid_level);

    double
    compute_viscosity(const std::vector<double>          &old_temperature,
                      const std::vector<double>          &old_old_temperature,
                      const std::vector<Tensor<1,dim> >  &old_temperature_grads,
                      const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
                      const std::vector<double>          &old_temperature_laplacians,
                      const std::vector<double>          &old_old_temperature_laplacians,
                      const std::vector<Tensor<1,dim> >  &old_velocity_values,
                      const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                      const std::vector<double>          &gamma_values,
                      const double                        global_u_infty,
                      const double                        global_T_variation,
                      const double                        cell_diameter) const;


    Triangulation<dim>                  triangulation;
    double                              global_Omega_diameter;

    const unsigned int                  stokes_degree;
    FESystem<dim>                       stokes_fe;
    DoFHandler<dim>                     stokes_dof_handler;
    ConstraintMatrix                    stokes_constraints;

    std::vector<types::global_dof_index> stokes_block_sizes;
    TrilinosWrappers::BlockSparseMatrix stokes_matrix;
    TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix;

    TrilinosWrappers::BlockVector       stokes_solution;
    TrilinosWrappers::BlockVector       old_stokes_solution;
    TrilinosWrappers::BlockVector       stokes_rhs;


    const unsigned int                  temperature_degree;
    FE_Q<dim>                           temperature_fe;
    DoFHandler<dim>                     temperature_dof_handler;
    ConstraintMatrix                    temperature_constraints;

    TrilinosWrappers::SparseMatrix      temperature_mass_matrix;
    TrilinosWrappers::SparseMatrix      temperature_stiffness_matrix;
    TrilinosWrappers::SparseMatrix      temperature_matrix;

    TrilinosWrappers::Vector            temperature_solution;
    TrilinosWrappers::Vector            old_temperature_solution;
    TrilinosWrappers::Vector            old_old_temperature_solution;
    TrilinosWrappers::Vector            temperature_rhs;


    double                              time_step;
    double                              old_time_step;
    unsigned int                        timestep_number;

    std_cxx11::shared_ptr<TrilinosWrappers::PreconditionAMG> Amg_preconditioner;
    std_cxx11::shared_ptr<TrilinosWrappers::PreconditionIC>  Mp_preconditioner;

    bool                                rebuild_stokes_matrix;
    bool                                rebuild_temperature_matrices;
    bool                                rebuild_stokes_preconditioner;
  };


  // @sect3{BoussinesqFlowProblem class implementation}

  // @sect4{BoussinesqFlowProblem::BoussinesqFlowProblem}
  //
  // The constructor of this class is an extension of the constructor in
  // step-22. We need to add the various variables that concern the
  // temperature. As discussed in the introduction, we are going to use
  // $Q_2\times Q_1$ (Taylor-Hood) elements again for the Stokes part, and
  // $Q_2$ elements for the temperature. However, by using variables that
  // store the polynomial degree of the Stokes and temperature finite
  // elements, it is easy to consistently modify the degree of the elements as
  // well as all quadrature formulas used on them downstream. Moreover, we
  // initialize the time stepping as well as the options for matrix assembly
  // and preconditioning:
  template <int dim>
  BoussinesqFlowProblem<dim>::BoussinesqFlowProblem ()
    :
    triangulation (Triangulation<dim>::maximum_smoothing),

    stokes_degree (1),
    stokes_fe (FE_Q<dim>(stokes_degree+1), dim,
               FE_Q<dim>(stokes_degree), 1),
    stokes_dof_handler (triangulation),

    temperature_degree (2),
    temperature_fe (temperature_degree),
    temperature_dof_handler (triangulation),

    time_step (0),
    old_time_step (0),
    timestep_number (0),
    rebuild_stokes_matrix (true),
    rebuild_temperature_matrices (true),
    rebuild_stokes_preconditioner (true)
  {}



  // @sect4{BoussinesqFlowProblem::get_maximal_velocity}

  // Starting the real functionality of this class is a helper function that
  // determines the maximum ($L_\infty$) velocity in the domain (at the
  // quadrature points, in fact). How it works should be relatively obvious to
  // all who have gotten to this point of the tutorial. Note that since we are
  // only interested in the velocity, rather than using
  // <code>stokes_fe_values.get_function_values</code> to get the values of
  // the entire Stokes solution (velocities and pressures) we use
  // <code>stokes_fe_values[velocities].get_function_values</code> to extract
  // only the velocities part. This has the additional benefit that we get it
  // as a Tensor<1,dim>, rather than some components in a Vector<double>,
  // allowing us to process it right away using the <code>norm()</code>
  // function to get the magnitude of the velocity.
  //
  // The only point worth thinking about a bit is how to choose the quadrature
  // points we use here. Since the goal of this function is to find the
  // maximal velocity over a domain by looking at quadrature points on each
  // cell. So we should ask how we should best choose these quadrature points
  // on each cell. To this end, recall that if we had a single $Q_1$ field
  // (rather than the vector-valued field of higher order) then the maximum
  // would be attained at a vertex of the mesh. In other words, we should use
  // the QTrapez class that has quadrature points only at the vertices of
  // cells.
  //
  // For higher order shape functions, the situation is more complicated: the
  // maxima and minima may be attained at points between the support points of
  // shape functions (for the usual $Q_p$ elements the support points are the
  // equidistant Lagrange interpolation points); furthermore, since we are
  // looking for the maximum magnitude of a vector-valued quantity, we can
  // even less say with certainty where the set of potential maximal points
  // are. Nevertheless, intuitively if not provably, the Lagrange
  // interpolation points appear to be a better choice than the Gauss points.
  //
  // There are now different methods to produce a quadrature formula with
  // quadrature points equal to the interpolation points of the finite
  // element. One option would be to use the
  // FiniteElement::get_unit_support_points() function, reduce the output to a
  // unique set of points to avoid duplicate function evaluations, and create
  // a Quadrature object using these points. Another option, chosen here, is
  // to use the QTrapez class and combine it with the QIterated class that
  // repeats the QTrapez formula on a number of sub-cells in each coordinate
  // direction. To cover all support points, we need to iterate it
  // <code>stokes_degree+1</code> times since this is the polynomial degree of
  // the Stokes element in use:
  template <int dim>
  double BoussinesqFlowProblem<dim>::get_maximal_velocity () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             stokes_degree+1);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (stokes_fe, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);
    double max_velocity = 0;

    const FEValuesExtractors::Vector velocities (0);

    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        fe_values[velocities].get_function_values (stokes_solution,
                                                   velocity_values);

        for (unsigned int q=0; q<n_q_points; ++q)
          max_velocity = std::max (max_velocity, velocity_values[q].norm());
      }

    return max_velocity;
  }




  // @sect4{BoussinesqFlowProblem::get_extrapolated_temperature_range}

  // Next a function that determines the minimum and maximum temperature at
  // quadrature points inside $\Omega$ when extrapolated from the two previous
  // time steps to the current one. We need this information in the
  // computation of the artificial viscosity parameter $\nu$ as discussed in
  // the introduction.
  //
  // The formula for the extrapolated temperature is
  // $\left(1+\frac{k_n}{k_{n-1}} \right)T^{n-1} + \frac{k_n}{k_{n-1}}
  // T^{n-2}$. The way to compute it is to loop over all quadrature points and
  // update the maximum and minimum value if the current value is
  // bigger/smaller than the previous one. We initialize the variables that
  // store the max and min before the loop over all quadrature points by the
  // smallest and the largest number representable as a double. Then we know
  // for a fact that it is larger/smaller than the minimum/maximum and that
  // the loop over all quadrature points is ultimately going to update the
  // initial value with the correct one.
  //
  // The only other complication worth mentioning here is that in the first
  // time step, $T^{k-2}$ is not yet available of course. In that case, we can
  // only use $T^{k-1}$ which we have from the initial temperature. As
  // quadrature points, we use the same choice as in the previous function
  // though with the difference that now the number of repetitions is
  // determined by the polynomial degree of the temperature field.
  template <int dim>
  std::pair<double,double>
  BoussinesqFlowProblem<dim>::get_extrapolated_temperature_range () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             temperature_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (temperature_fe, quadrature_formula,
                             update_values);
    std::vector<double> old_temperature_values(n_q_points);
    std::vector<double> old_old_temperature_values(n_q_points);

    if (timestep_number != 0)
      {
        double min_temperature = std::numeric_limits<double>::max(),
               max_temperature = -std::numeric_limits<double>::max();

        typename DoFHandler<dim>::active_cell_iterator
        cell = temperature_dof_handler.begin_active(),
        endc = temperature_dof_handler.end();
        for (; cell!=endc; ++cell)
          {
            fe_values.reinit (cell);
            fe_values.get_function_values (old_temperature_solution,
                                           old_temperature_values);
            fe_values.get_function_values (old_old_temperature_solution,
                                           old_old_temperature_values);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double temperature =
                  (1. + time_step/old_time_step) * old_temperature_values[q]-
                  time_step/old_time_step * old_old_temperature_values[q];

                min_temperature = std::min (min_temperature, temperature);
                max_temperature = std::max (max_temperature, temperature);
              }
          }

        return std::make_pair(min_temperature, max_temperature);
      }
    else
      {
        double min_temperature = std::numeric_limits<double>::max(),
               max_temperature = -std::numeric_limits<double>::max();

        typename DoFHandler<dim>::active_cell_iterator
        cell = temperature_dof_handler.begin_active(),
        endc = temperature_dof_handler.end();
        for (; cell!=endc; ++cell)
          {
            fe_values.reinit (cell);
            fe_values.get_function_values (old_temperature_solution,
                                           old_temperature_values);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double temperature = old_temperature_values[q];

                min_temperature = std::min (min_temperature, temperature);
                max_temperature = std::max (max_temperature, temperature);
              }
          }

        return std::make_pair(min_temperature, max_temperature);
      }
  }



  // @sect4{BoussinesqFlowProblem::compute_viscosity}

  // The last of the tool functions computes the artificial viscosity
  // parameter $\nu|_K$ on a cell $K$ as a function of the extrapolated
  // temperature, its gradient and Hessian (second derivatives), the velocity,
  // the right hand side $\gamma$ all on the quadrature points of the current
  // cell, and various other parameters as described in detail in the
  // introduction.
  //
  // There are some universal constants worth mentioning here. First, we need
  // to fix $\beta$; we choose $\beta=0.017\cdot dim$, a choice discussed in
  // detail in the results section of this tutorial program. The second is the
  // exponent $\alpha$; $\alpha=1$ appears to work fine for the current
  // program, even though some additional benefit might be expected from
  // choosing $\alpha = 2$. Finally, there is one thing that requires special
  // casing: In the first time step, the velocity equals zero, and the formula
  // for $\nu|_K$ is not defined. In that case, we return $\nu|_K=5\cdot 10^3
  // \cdot h_K$, a choice admittedly more motivated by heuristics than
  // anything else (it is in the same order of magnitude, however, as the
  // value returned for most cells on the second time step).
  //
  // The rest of the function should be mostly obvious based on the material
  // discussed in the introduction:
  template <int dim>
  double
  BoussinesqFlowProblem<dim>::
  compute_viscosity (const std::vector<double>          &old_temperature,
                     const std::vector<double>          &old_old_temperature,
                     const std::vector<Tensor<1,dim> >  &old_temperature_grads,
                     const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
                     const std::vector<double>          &old_temperature_laplacians,
                     const std::vector<double>          &old_old_temperature_laplacians,
                     const std::vector<Tensor<1,dim> >  &old_velocity_values,
                     const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                     const std::vector<double>          &gamma_values,
                     const double                        global_u_infty,
                     const double                        global_T_variation,
                     const double                        cell_diameter) const
  {
    const double beta = 0.017 * dim;
    const double alpha = 1;

    if (global_u_infty == 0)
      return 5e-3 * cell_diameter;

    const unsigned int n_q_points = old_temperature.size();

    double max_residual = 0;
    double max_velocity = 0;

    for (unsigned int q=0; q < n_q_points; ++q)
      {
        const Tensor<1,dim> u = (old_velocity_values[q] +
                                 old_old_velocity_values[q]) / 2;

        const double dT_dt = (old_temperature[q] - old_old_temperature[q])
                             / old_time_step;
        const double u_grad_T = u * (old_temperature_grads[q] +
                                     old_old_temperature_grads[q]) / 2;

        const double kappa_Delta_T = EquationData::kappa
                                     * (old_temperature_laplacians[q] +
                                        old_old_temperature_laplacians[q]) / 2;

        const double residual
          = std::abs((dT_dt + u_grad_T - kappa_Delta_T - gamma_values[q]) *
                     std::pow((old_temperature[q]+old_old_temperature[q]) / 2,
                              alpha-1.));

        max_residual = std::max (residual,        max_residual);
        max_velocity = std::max (std::sqrt (u*u), max_velocity);
      }

    const double c_R = std::pow (2., (4.-2*alpha)/dim);
    const double global_scaling = c_R * global_u_infty * global_T_variation *
                                  std::pow(global_Omega_diameter, alpha - 2.);

    return (beta *
            max_velocity *
            std::min (cell_diameter,
                      std::pow(cell_diameter,alpha) *
                      max_residual / global_scaling));
  }



  // @sect4{BoussinesqFlowProblem::setup_dofs}
  //
  // This is the function that sets up the DoFHandler objects we have here
  // (one for the Stokes part and one for the temperature part) as well as set
  // to the right sizes the various objects required for the linear algebra in
  // this program. Its basic operations are similar to what we do in step-22.
  //
  // The body of the function first enumerates all degrees of freedom for the
  // Stokes and temperature systems. For the Stokes part, degrees of freedom
  // are then sorted to ensure that velocities precede pressure DoFs so that
  // we can partition the Stokes matrix into a $2\times 2$ matrix. As a
  // difference to step-22, we do not perform any additional DoF
  // renumbering. In that program, it paid off since our solver was heavily
  // dependent on ILU's, whereas we use AMG here which is not sensitive to the
  // DoF numbering. The IC preconditioner for the inversion of the pressure
  // mass matrix would of course take advantage of a Cuthill-McKee like
  // renumbering, but its costs are low compared to the velocity portion, so
  // the additional work does not pay off.
  //
  // We then proceed with the generation of the hanging node constraints that
  // arise from adaptive grid refinement for both DoFHandler objects. For the
  // velocity, we impose no-flux boundary conditions $\mathbf{u}\cdot
  // \mathbf{n}=0$ by adding constraints to the object that already stores the
  // hanging node constraints matrix. The second parameter in the function
  // describes the first of the velocity components in the total dof vector,
  // which is zero here. The variable <code>no_normal_flux_boundaries</code>
  // denotes the boundary indicators for which to set the no flux boundary
  // conditions; here, this is boundary indicator zero.
  //
  // After having done so, we count the number of degrees of freedom in the
  // various blocks:
  template <int dim>
  void BoussinesqFlowProblem<dim>::setup_dofs ()
  {
    std::vector<unsigned int> stokes_sub_blocks (dim+1,0);
    stokes_sub_blocks[dim] = 1;

    {
      stokes_dof_handler.distribute_dofs (stokes_fe);
      DoFRenumbering::component_wise (stokes_dof_handler, stokes_sub_blocks);

      stokes_constraints.clear ();
      DoFTools::make_hanging_node_constraints (stokes_dof_handler,
                                               stokes_constraints);
      std::set<types::boundary_id> no_normal_flux_boundaries;
      no_normal_flux_boundaries.insert (0);
      VectorTools::compute_no_normal_flux_constraints (stokes_dof_handler, 0,
                                                       no_normal_flux_boundaries,
                                                       stokes_constraints);
      stokes_constraints.close ();
    }
    {
      temperature_dof_handler.distribute_dofs (temperature_fe);

      temperature_constraints.clear ();
      DoFTools::make_hanging_node_constraints (temperature_dof_handler,
                                               temperature_constraints);
      temperature_constraints.close ();
    }

    std::vector<types::global_dof_index> stokes_dofs_per_block (2);
    DoFTools::count_dofs_per_block (stokes_dof_handler, stokes_dofs_per_block,
                                    stokes_sub_blocks);

    const unsigned int n_u = stokes_dofs_per_block[0],
                       n_p = stokes_dofs_per_block[1],
                       n_T = temperature_dof_handler.n_dofs();

    std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << " (on "
              << triangulation.n_levels()
              << " levels)"
              << std::endl
              << "Number of degrees of freedom: "
              << n_u + n_p + n_T
              << " (" << n_u << '+' << n_p << '+'<< n_T <<')'
              << std::endl
              << std::endl;

    // The next step is to create the sparsity pattern for the Stokes and
    // temperature system matrices as well as the preconditioner matrix from
    // which we build the Stokes preconditioner. As in step-22, we choose to
    // create the pattern not as in the first few tutorial programs, but by
    // using the blocked version of CompressedSimpleSparsityPattern.  The
    // reason for doing this is mainly memory, that is, the SparsityPattern
    // class would consume too much memory when used in three spatial
    // dimensions as we intend to do for this program.
    //
    // So, we first release the memory stored in the matrices, then set up an
    // object of type BlockCompressedSimpleSparsityPattern consisting of
    // $2\times 2$ blocks (for the Stokes system matrix and preconditioner) or
    // CompressedSimpleSparsityPattern (for the temperature part). We then
    // fill these objects with the nonzero pattern, taking into account that
    // for the Stokes system matrix, there are no entries in the
    // pressure-pressure block (but all velocity vector components couple with
    // each other and with the pressure). Similarly, in the Stokes
    // preconditioner matrix, only the diagonal blocks are nonzero, since we
    // use the vector Laplacian as discussed in the introduction. This
    // operator only couples each vector component of the Laplacian with
    // itself, but not with the other vector components. (Application of the
    // constraints resulting from the no-flux boundary conditions will couple
    // vector components at the boundary again, however.)
    //
    // When generating the sparsity pattern, we directly apply the constraints
    // from hanging nodes and no-flux boundary conditions. This approach was
    // already used in step-27, but is different from the one in early
    // tutorial programs where we first built the original sparsity pattern
    // and only then added the entries resulting from constraints. The reason
    // for doing so is that later during assembly we are going to distribute
    // the constraints immediately when transferring local to global
    // dofs. Consequently, there will be no data written at positions of
    // constrained degrees of freedom, so we can let the
    // DoFTools::make_sparsity_pattern function omit these entries by setting
    // the last Boolean flag to <code>false</code>. Once the sparsity pattern
    // is ready, we can use it to initialize the Trilinos matrices. Since the
    // Trilinos matrices store the sparsity pattern internally, there is no
    // need to keep the sparsity pattern around after the initialization of
    // the matrix.
    stokes_block_sizes.resize (2);
    stokes_block_sizes[0] = n_u;
    stokes_block_sizes[1] = n_p;
    {
      stokes_matrix.clear ();

      BlockCompressedSimpleSparsityPattern csp (2,2);

      csp.block(0,0).reinit (n_u, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(1,1).reinit (n_p, n_p);

      csp.collect_sizes ();

      Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);

      for (unsigned int c=0; c<dim+1; ++c)
        for (unsigned int d=0; d<dim+1; ++d)
          if (! ((c==dim) && (d==dim)))
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      DoFTools::make_sparsity_pattern (stokes_dof_handler, coupling, csp,
                                       stokes_constraints, false);

      stokes_matrix.reinit (csp);
    }

    {
      Amg_preconditioner.reset ();
      Mp_preconditioner.reset ();
      stokes_preconditioner_matrix.clear ();

      BlockCompressedSimpleSparsityPattern csp (2,2);

      csp.block(0,0).reinit (n_u, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(1,1).reinit (n_p, n_p);

      csp.collect_sizes ();

      Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
      for (unsigned int c=0; c<dim+1; ++c)
        for (unsigned int d=0; d<dim+1; ++d)
          if (c == d)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      DoFTools::make_sparsity_pattern (stokes_dof_handler, coupling, csp,
                                       stokes_constraints, false);

      stokes_preconditioner_matrix.reinit (csp);
    }

    // The creation of the temperature matrix (or, rather, matrices, since we
    // provide a temperature mass matrix and a temperature stiffness matrix,
    // that will be added together for time discretization) follows the
    // generation of the Stokes matrix &ndash; except that it is much easier
    // here since we do not need to take care of any blocks or coupling
    // between components. Note how we initialize the three temperature
    // matrices: We only use the sparsity pattern for reinitialization of the
    // first matrix, whereas we use the previously generated matrix for the
    // two remaining reinits. The reason for doing so is that reinitialization
    // from an already generated matrix allows Trilinos to reuse the sparsity
    // pattern instead of generating a new one for each copy. This saves both
    // some time and memory.
    {
      temperature_mass_matrix.clear ();
      temperature_stiffness_matrix.clear ();
      temperature_matrix.clear ();

      CompressedSimpleSparsityPattern csp (n_T, n_T);
      DoFTools::make_sparsity_pattern (temperature_dof_handler, csp,
                                       temperature_constraints, false);

      temperature_matrix.reinit (csp);
      temperature_mass_matrix.reinit (temperature_matrix);
      temperature_stiffness_matrix.reinit (temperature_matrix);
    }

    // Lastly, we set the vectors for the Stokes solutions $\mathbf u^{n-1}$
    // and $\mathbf u^{n-2}$, as well as for the temperatures $T^{n}$,
    // $T^{n-1}$ and $T^{n-2}$ (required for time stepping) and all the system
    // right hand sides to their correct sizes and block structure:
    stokes_solution.reinit (stokes_block_sizes);
    old_stokes_solution.reinit (stokes_block_sizes);
    stokes_rhs.reinit (stokes_block_sizes);

    temperature_solution.reinit (temperature_dof_handler.n_dofs());
    old_temperature_solution.reinit (temperature_dof_handler.n_dofs());
    old_old_temperature_solution.reinit (temperature_dof_handler.n_dofs());

    temperature_rhs.reinit (temperature_dof_handler.n_dofs());
  }



  // @sect4{BoussinesqFlowProblem::assemble_stokes_preconditioner}
  //
  // This function assembles the matrix we use for preconditioning the Stokes
  // system. What we need are a vector Laplace matrix on the velocity
  // components and a mass matrix weighted by $\eta^{-1}$ on the pressure
  // component. We start by generating a quadrature object of appropriate
  // order, the FEValues object that can give values and gradients at the
  // quadrature points (together with quadrature weights). Next we create data
  // structures for the cell matrix and the relation between local and global
  // DoFs. The vectors <code>grad_phi_u</code> and <code>phi_p</code> are
  // going to hold the values of the basis functions in order to faster build
  // up the local matrices, as was already done in step-22. Before we start
  // the loop over all active cells, we have to specify which components are
  // pressure and which are velocity.
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::assemble_stokes_preconditioner ()
  {
    stokes_preconditioner_matrix = 0;

    const QGauss<dim> quadrature_formula(stokes_degree+2);
    FEValues<dim>     stokes_fe_values (stokes_fe, quadrature_formula,
                                        update_JxW_values |
                                        update_values |
                                        update_gradients);

    const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<Tensor<2,dim> > grad_phi_u (dofs_per_cell);
    std::vector<double>         phi_p      (dofs_per_cell);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        stokes_fe_values.reinit (cell);
        local_matrix = 0;

        // The creation of the local matrix is rather simple. There are only a
        // Laplace term (on the velocity) and a mass matrix weighted by
        // $\eta^{-1}$ to be generated, so the creation of the local matrix is
        // done in two lines. Once the local matrix is ready (loop over rows
        // and columns in the local matrix on each quadrature point), we get
        // the local DoF indices and write the local information into the
        // global matrix. We do this as in step-27, i.e. we directly apply the
        // constraints from hanging nodes locally. By doing so, we don't have
        // to do that afterwards, and we don't also write into entries of the
        // matrix that will actually be set to zero again later when
        // eliminating constraints.
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                grad_phi_u[k] = stokes_fe_values[velocities].gradient(k,q);
                phi_p[k]      = stokes_fe_values[pressure].value (k, q);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                local_matrix(i,j) += (EquationData::eta *
                                      scalar_product (grad_phi_u[i], grad_phi_u[j])
                                      +
                                      (1./EquationData::eta) *
                                      phi_p[i] * phi_p[j])
                                     * stokes_fe_values.JxW(q);
          }

        cell->get_dof_indices (local_dof_indices);
        stokes_constraints.distribute_local_to_global (local_matrix,
                                                       local_dof_indices,
                                                       stokes_preconditioner_matrix);
      }
  }



  // @sect4{BoussinesqFlowProblem::build_stokes_preconditioner}
  //
  // This function generates the inner preconditioners that are going to be
  // used for the Schur complement block preconditioner. Since the
  // preconditioners need only to be regenerated when the matrices change,
  // this function does not have to do anything in case the matrices have not
  // changed (i.e., the flag <code>rebuild_stokes_preconditioner</code> has
  // the value <code>false</code>). Otherwise its first task is to call
  // <code>assemble_stokes_preconditioner</code> to generate the
  // preconditioner matrices.
  //
  // Next, we set up the preconditioner for the velocity-velocity matrix
  // <i>A</i>. As explained in the introduction, we are going to use an AMG
  // preconditioner based on a vector Laplace matrix $\hat{A}$ (which is
  // spectrally close to the Stokes matrix <i>A</i>). Usually, the
  // TrilinosWrappers::PreconditionAMG class can be seen as a good black-box
  // preconditioner which does not need any special knowledge. In this case,
  // however, we have to be careful: since we build an AMG for a vector
  // problem, we have to tell the preconditioner setup which dofs belong to
  // which vector component. We do this using the function
  // DoFTools::extract_constant_modes, a function that generates a set of
  // <code>dim</code> vectors, where each one has ones in the respective
  // component of the vector problem and zeros elsewhere. Hence, these are the
  // constant modes on each component, which explains the name of the
  // variable.
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::build_stokes_preconditioner ()
  {
    if (rebuild_stokes_preconditioner == false)
      return;

    std::cout << "   Rebuilding Stokes preconditioner..." << std::flush;

    assemble_stokes_preconditioner ();

    Amg_preconditioner = std_cxx11::shared_ptr<TrilinosWrappers::PreconditionAMG>
                         (new TrilinosWrappers::PreconditionAMG());

    std::vector<std::vector<bool> > constant_modes;
    FEValuesExtractors::Vector velocity_components(0);
    DoFTools::extract_constant_modes (stokes_dof_handler,
                                      stokes_fe.component_mask(velocity_components),
                                      constant_modes);
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
    amg_data.constant_modes = constant_modes;

    // Next, we set some more options of the AMG preconditioner. In
    // particular, we need to tell the AMG setup that we use quadratic basis
    // functions for the velocity matrix (this implies more nonzero elements
    // in the matrix, so that a more robust algorithm needs to be chosen
    // internally). Moreover, we want to be able to control how the coarsening
    // structure is build up. The way the Trilinos smoothed aggregation AMG
    // does this is to look which matrix entries are of similar size as the
    // diagonal entry in order to algebraically build a coarse-grid
    // structure. By setting the parameter <code>aggregation_threshold</code>
    // to 0.02, we specify that all entries that are more than two percent of
    // size of some diagonal pivots in that row should form one coarse grid
    // point. This parameter is rather ad hoc, and some fine-tuning of it can
    // influence the performance of the preconditioner. As a rule of thumb,
    // larger values of <code>aggregation_threshold</code> will decrease the
    // number of iterations, but increase the costs per iteration. A look at
    // the Trilinos documentation will provide more information on these
    // parameters. With this data set, we then initialize the preconditioner
    // with the matrix we want it to apply to.
    //
    // Finally, we also initialize the preconditioner for the inversion of the
    // pressure mass matrix. This matrix is symmetric and well-behaved, so we
    // can chose a simple preconditioner. We stick with an incomplete Cholesky
    // (IC) factorization preconditioner, which is designed for symmetric
    // matrices. We could have also chosen an SSOR preconditioner with
    // relaxation factor around 1.2, but IC is cheaper for our example. We
    // wrap the preconditioners into a <code>std_cxx11::shared_ptr</code>
    // pointer, which makes it easier to recreate the preconditioner next time
    // around since we do not have to care about destroying the previously
    // used object.
    amg_data.elliptic = true;
    amg_data.higher_order_elements = true;
    amg_data.smoother_sweeps = 2;
    amg_data.aggregation_threshold = 0.02;
    Amg_preconditioner->initialize(stokes_preconditioner_matrix.block(0,0),
                                   amg_data);

    Mp_preconditioner = std_cxx11::shared_ptr<TrilinosWrappers::PreconditionIC>
                        (new TrilinosWrappers::PreconditionIC());
    Mp_preconditioner->initialize(stokes_preconditioner_matrix.block(1,1));

    std::cout << std::endl;

    rebuild_stokes_preconditioner = false;
  }



  // @sect4{BoussinesqFlowProblem::assemble_stokes_system}
  //
  // The time lag scheme we use for advancing the coupled Stokes-temperature
  // system forces us to split up the assembly (and the solution of linear
  // systems) into two step. The first one is to create the Stokes system
  // matrix and right hand side, and the second is to create matrix and right
  // hand sides for the temperature dofs, which depends on the result of the
  // linear system for the velocity.
  //
  // This function is called at the beginning of each time step. In the first
  // time step or if the mesh has changed, indicated by the
  // <code>rebuild_stokes_matrix</code>, we need to assemble the Stokes
  // matrix; on the other hand, if the mesh hasn't changed and the matrix is
  // already available, this is not necessary and all we need to do is
  // assemble the right hand side vector which changes in each time step.
  //
  // Regarding the technical details of implementation, not much has changed
  // from step-22. We reset matrix and vector, create a quadrature formula on
  // the cells, and then create the respective FEValues object. For the update
  // flags, we require basis function derivatives only in case of a full
  // assembly, since they are not needed for the right hand side; as always,
  // choosing the minimal set of flags depending on what is currently needed
  // makes the call to FEValues::reinit further down in the program more
  // efficient.
  //
  // There is one thing that needs to be commented &ndash; since we have a
  // separate finite element and DoFHandler for the temperature, we need to
  // generate a second FEValues object for the proper evaluation of the
  // temperature solution. This isn't too complicated to realize here: just
  // use the temperature structures and set an update flag for the basis
  // function values which we need for evaluation of the temperature
  // solution. The only important part to remember here is that the same
  // quadrature formula is used for both FEValues objects to ensure that we
  // get matching information when we loop over the quadrature points of the
  // two objects.
  //
  // The declarations proceed with some shortcuts for array sizes, the
  // creation of the local matrix and right hand side as well as the vector
  // for the indices of the local dofs compared to the global system.
  template <int dim>
  void BoussinesqFlowProblem<dim>::assemble_stokes_system ()
  {
    std::cout << "   Assembling..." << std::flush;

    if (rebuild_stokes_matrix == true)
      stokes_matrix=0;

    stokes_rhs=0;

    const QGauss<dim> quadrature_formula (stokes_degree+2);
    FEValues<dim>     stokes_fe_values (stokes_fe, quadrature_formula,
                                        update_values    |
                                        update_quadrature_points  |
                                        update_JxW_values |
                                        (rebuild_stokes_matrix == true
                                         ?
                                         update_gradients
                                         :
                                         UpdateFlags(0)));

    FEValues<dim>     temperature_fe_values (temperature_fe, quadrature_formula,
                                             update_values);

    const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs    (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    // Next we need a vector that will contain the values of the temperature
    // solution at the previous time level at the quadrature points to
    // assemble the source term in the right hand side of the momentum
    // equation. Let's call this vector <code>old_solution_values</code>.
    //
    // The set of vectors we create next hold the evaluations of the basis
    // functions as well as their gradients and symmetrized gradients that
    // will be used for creating the matrices. Putting these into their own
    // arrays rather than asking the FEValues object for this information each
    // time it is needed is an optimization to accelerate the assembly
    // process, see step-22 for details.
    //
    // The last two declarations are used to extract the individual blocks
    // (velocity, pressure, temperature) from the total FE system.
    std::vector<double>               old_temperature_values(n_q_points);

    std::vector<Tensor<1,dim> >          phi_u       (dofs_per_cell);
    std::vector<SymmetricTensor<2,dim> > grads_phi_u (dofs_per_cell);
    std::vector<double>                  div_phi_u   (dofs_per_cell);
    std::vector<double>                  phi_p       (dofs_per_cell);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    // Now start the loop over all cells in the problem. We are working on two
    // different DoFHandlers for this assembly routine, so we must have two
    // different cell iterators for the two objects in use. This might seem a
    // bit peculiar, since both the Stokes system and the temperature system
    // use the same grid, but that's the only way to keep degrees of freedom
    // in sync. The first statements within the loop are again all very
    // familiar, doing the update of the finite element data as specified by
    // the update flags, zeroing out the local arrays and getting the values
    // of the old solution at the quadrature points. Then we are ready to loop
    // over the quadrature points on the cell.
    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
    temperature_cell = temperature_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++temperature_cell)
      {
        stokes_fe_values.reinit (cell);
        temperature_fe_values.reinit (temperature_cell);

        local_matrix = 0;
        local_rhs = 0;

        temperature_fe_values.get_function_values (old_temperature_solution,
                                                   old_temperature_values);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            const double old_temperature = old_temperature_values[q];

            // Next we extract the values and gradients of basis functions
            // relevant to the terms in the inner products. As shown in
            // step-22 this helps accelerate assembly.
            //
            // Once this is done, we start the loop over the rows and columns
            // of the local matrix and feed the matrix with the relevant
            // products. The right hand side is filled with the forcing term
            // driven by temperature in direction of gravity (which is
            // vertical in our example).  Note that the right hand side term
            // is always generated, whereas the matrix contributions are only
            // updated when it is requested by the
            // <code>rebuild_matrices</code> flag.
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                phi_u[k] = stokes_fe_values[velocities].value (k,q);
                if (rebuild_stokes_matrix)
                  {
                    grads_phi_u[k] = stokes_fe_values[velocities].symmetric_gradient(k,q);
                    div_phi_u[k]   = stokes_fe_values[velocities].divergence (k, q);
                    phi_p[k]       = stokes_fe_values[pressure].value (k, q);
                  }
              }

            if (rebuild_stokes_matrix)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  local_matrix(i,j) += (EquationData::eta * 2 *
                                        (grads_phi_u[i] * grads_phi_u[j])
                                        - div_phi_u[i] * phi_p[j]
                                        - phi_p[i] * div_phi_u[j])
                                       * stokes_fe_values.JxW(q);

            const Point<dim> gravity = -( (dim == 2) ? (Point<dim> (0,1)) :
                                          (Point<dim> (0,0,1)) );
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              local_rhs(i) += (-EquationData::density *
                               EquationData::beta *
                               gravity * phi_u[i] * old_temperature)*
                              stokes_fe_values.JxW(q);
          }

        // The last step in the loop over all cells is to enter the local
        // contributions into the global matrix and vector structures to the
        // positions specified in <code>local_dof_indices</code>.  Again, we
        // let the ConstraintMatrix class do the insertion of the cell matrix
        // elements to the global matrix, which already condenses the hanging
        // node constraints.
        cell->get_dof_indices (local_dof_indices);

        if (rebuild_stokes_matrix == true)
          stokes_constraints.distribute_local_to_global (local_matrix,
                                                         local_rhs,
                                                         local_dof_indices,
                                                         stokes_matrix,
                                                         stokes_rhs);
        else
          stokes_constraints.distribute_local_to_global (local_rhs,
                                                         local_dof_indices,
                                                         stokes_rhs);
      }

    rebuild_stokes_matrix = false;

    std::cout << std::endl;
  }




  // @sect4{BoussinesqFlowProblem::assemble_temperature_matrix}
  //
  // This function assembles the matrix in the temperature equation. The
  // temperature matrix consists of two parts, a mass matrix and the time step
  // size times a stiffness matrix given by a Laplace term times the amount of
  // diffusion. Since the matrix depends on the time step size (which varies
  // from one step to another), the temperature matrix needs to be updated
  // every time step. We could simply regenerate the matrices in every time
  // step, but this is not really efficient since mass and Laplace matrix do
  // only change when we change the mesh. Hence, we do this more efficiently
  // by generating two separate matrices in this function, one for the mass
  // matrix and one for the stiffness (diffusion) matrix. We will then sum up
  // the matrix plus the stiffness matrix times the time step size once we
  // know the actual time step.
  //
  // So the details for this first step are very simple. In case we need to
  // rebuild the matrix (i.e., the mesh has changed), we zero the data
  // structures, get a quadrature formula and a FEValues object, and create
  // local matrices, local dof indices and evaluation structures for the basis
  // functions.
  template <int dim>
  void BoussinesqFlowProblem<dim>::assemble_temperature_matrix ()
  {
    if (rebuild_temperature_matrices == false)
      return;

    temperature_mass_matrix = 0;
    temperature_stiffness_matrix = 0;

    QGauss<dim>   quadrature_formula (temperature_degree+2);
    FEValues<dim> temperature_fe_values (temperature_fe, quadrature_formula,
                                         update_values    | update_gradients |
                                         update_JxW_values);

    const unsigned int   dofs_per_cell   = temperature_fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   local_mass_matrix (dofs_per_cell, dofs_per_cell);
    FullMatrix<double>   local_stiffness_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    std::vector<double>         phi_T       (dofs_per_cell);
    std::vector<Tensor<1,dim> > grad_phi_T  (dofs_per_cell);

    // Now, let's start the loop over all cells in the triangulation. We need
    // to zero out the local matrices, update the finite element evaluations,
    // and then loop over the rows and columns of the matrices on each
    // quadrature point, where we then create the mass matrix and the
    // stiffness matrix (Laplace terms times the diffusion
    // <code>EquationData::kappa</code>. Finally, we let the constraints
    // object insert these values into the global matrix, and directly
    // condense the constraints into the matrix.
    typename DoFHandler<dim>::active_cell_iterator
    cell = temperature_dof_handler.begin_active(),
    endc = temperature_dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        local_mass_matrix = 0;
        local_stiffness_matrix = 0;

        temperature_fe_values.reinit (cell);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                grad_phi_T[k] = temperature_fe_values.shape_grad (k,q);
                phi_T[k]      = temperature_fe_values.shape_value (k, q);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  local_mass_matrix(i,j)
                  += (phi_T[i] * phi_T[j]
                      *
                      temperature_fe_values.JxW(q));
                  local_stiffness_matrix(i,j)
                  += (EquationData::kappa * grad_phi_T[i] * grad_phi_T[j]
                      *
                      temperature_fe_values.JxW(q));
                }
          }

        cell->get_dof_indices (local_dof_indices);

        temperature_constraints.distribute_local_to_global (local_mass_matrix,
                                                            local_dof_indices,
                                                            temperature_mass_matrix);
        temperature_constraints.distribute_local_to_global (local_stiffness_matrix,
                                                            local_dof_indices,
                                                            temperature_stiffness_matrix);
      }

    rebuild_temperature_matrices = false;
  }



  // @sect4{BoussinesqFlowProblem::assemble_temperature_system}
  //
  // This function does the second part of the assembly work on the
  // temperature matrix, the actual addition of pressure mass and stiffness
  // matrix (where the time step size comes into play), as well as the
  // creation of the velocity-dependent right hand side. The declarations for
  // the right hand side assembly in this function are pretty much the same as
  // the ones used in the other assembly routines, except that we restrict
  // ourselves to vectors this time. We are going to calculate residuals on
  // the temperature system, which means that we have to evaluate second
  // derivatives, specified by the update flag <code>update_hessians</code>.
  //
  // The temperature equation is coupled to the Stokes system by means of the
  // fluid velocity. These two parts of the solution are associated with
  // different DoFHandlers, so we again need to create a second FEValues
  // object for the evaluation of the velocity at the quadrature points.
  template <int dim>
  void BoussinesqFlowProblem<dim>::
  assemble_temperature_system (const double maximal_velocity)
  {
    const bool use_bdf2_scheme = (timestep_number != 0);

    if (use_bdf2_scheme == true)
      {
        temperature_matrix.copy_from (temperature_mass_matrix);
        temperature_matrix *= (2*time_step + old_time_step) /
                              (time_step + old_time_step);
        temperature_matrix.add (time_step, temperature_stiffness_matrix);
      }
    else
      {
        temperature_matrix.copy_from (temperature_mass_matrix);
        temperature_matrix.add (time_step, temperature_stiffness_matrix);
      }

    temperature_rhs = 0;

    const QGauss<dim> quadrature_formula(temperature_degree+2);
    FEValues<dim>     temperature_fe_values (temperature_fe, quadrature_formula,
                                             update_values    |
                                             update_gradients |
                                             update_hessians  |
                                             update_quadrature_points  |
                                             update_JxW_values);
    FEValues<dim>     stokes_fe_values (stokes_fe, quadrature_formula,
                                        update_values);

    const unsigned int   dofs_per_cell   = temperature_fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    // Next comes the declaration of vectors to hold the old and older
    // solution values (as a notation for time levels <i>n-1</i> and
    // <i>n-2</i>, respectively) and gradients at quadrature points of the
    // current cell. We also declare an object to hold the temperature right
    // hand side values (<code>gamma_values</code>), and we again use
    // shortcuts for the temperature basis functions. Eventually, we need to
    // find the temperature extrema and the diameter of the computational
    // domain which will be used for the definition of the stabilization
    // parameter (we got the maximal velocity as an input to this function).
    std::vector<Tensor<1,dim> > old_velocity_values (n_q_points);
    std::vector<Tensor<1,dim> > old_old_velocity_values (n_q_points);
    std::vector<double>         old_temperature_values (n_q_points);
    std::vector<double>         old_old_temperature_values(n_q_points);
    std::vector<Tensor<1,dim> > old_temperature_grads(n_q_points);
    std::vector<Tensor<1,dim> > old_old_temperature_grads(n_q_points);
    std::vector<double>         old_temperature_laplacians(n_q_points);
    std::vector<double>         old_old_temperature_laplacians(n_q_points);

    EquationData::TemperatureRightHandSide<dim>  temperature_right_hand_side;
    std::vector<double> gamma_values (n_q_points);

    std::vector<double>         phi_T      (dofs_per_cell);
    std::vector<Tensor<1,dim> > grad_phi_T (dofs_per_cell);

    const std::pair<double,double>
    global_T_range = get_extrapolated_temperature_range();

    const FEValuesExtractors::Vector velocities (0);

    // Now, let's start the loop over all cells in the triangulation. Again,
    // we need two cell iterators that walk in parallel through the cells of
    // the two involved DoFHandler objects for the Stokes and temperature
    // part. Within the loop, we first set the local rhs to zero, and then get
    // the values and derivatives of the old solution functions at the
    // quadrature points, since they are going to be needed for the definition
    // of the stabilization parameters and as coefficients in the equation,
    // respectively. Note that since the temperature has its own DoFHandler
    // and FEValues object we get the entire solution at the quadrature point
    // (which is the scalar temperature field only anyway) whereas for the
    // Stokes part we restrict ourselves to extracting the velocity part (and
    // ignoring the pressure part) by using
    // <code>stokes_fe_values[velocities].get_function_values</code>.
    typename DoFHandler<dim>::active_cell_iterator
    cell = temperature_dof_handler.begin_active(),
    endc = temperature_dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
    stokes_cell = stokes_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++stokes_cell)
      {
        local_rhs = 0;

        temperature_fe_values.reinit (cell);
        stokes_fe_values.reinit (stokes_cell);

        temperature_fe_values.get_function_values (old_temperature_solution,
                                                   old_temperature_values);
        temperature_fe_values.get_function_values (old_old_temperature_solution,
                                                   old_old_temperature_values);

        temperature_fe_values.get_function_gradients (old_temperature_solution,
                                                      old_temperature_grads);
        temperature_fe_values.get_function_gradients (old_old_temperature_solution,
                                                      old_old_temperature_grads);

        temperature_fe_values.get_function_laplacians (old_temperature_solution,
                                                       old_temperature_laplacians);
        temperature_fe_values.get_function_laplacians (old_old_temperature_solution,
                                                       old_old_temperature_laplacians);

        temperature_right_hand_side.value_list (temperature_fe_values.get_quadrature_points(),
                                                gamma_values);

        stokes_fe_values[velocities].get_function_values (stokes_solution,
                                                          old_velocity_values);
        stokes_fe_values[velocities].get_function_values (old_stokes_solution,
                                                          old_old_velocity_values);

        // Next, we calculate the artificial viscosity for stabilization
        // according to the discussion in the introduction using the dedicated
        // function. With that at hand, we can get into the loop over
        // quadrature points and local rhs vector components. The terms here
        // are quite lengthy, but their definition follows the time-discrete
        // system developed in the introduction of this program. The BDF-2
        // scheme needs one more term from the old time step (and involves
        // more complicated factors) than the backward Euler scheme that is
        // used for the first time step. When all this is done, we distribute
        // the local vector into the global one (including hanging node
        // constraints).
        const double nu
          = compute_viscosity (old_temperature_values,
                               old_old_temperature_values,
                               old_temperature_grads,
                               old_old_temperature_grads,
                               old_temperature_laplacians,
                               old_old_temperature_laplacians,
                               old_velocity_values,
                               old_old_velocity_values,
                               gamma_values,
                               maximal_velocity,
                               global_T_range.second - global_T_range.first,
                               cell->diameter());

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                grad_phi_T[k] = temperature_fe_values.shape_grad (k,q);
                phi_T[k]      = temperature_fe_values.shape_value (k, q);
              }

            const double T_term_for_rhs
              = (use_bdf2_scheme ?
                 (old_temperature_values[q] *
                  (1 + time_step/old_time_step)
                  -
                  old_old_temperature_values[q] *
                  (time_step * time_step) /
                  (old_time_step * (time_step + old_time_step)))
                 :
                 old_temperature_values[q]);

            const Tensor<1,dim> ext_grad_T
              = (use_bdf2_scheme ?
                 (old_temperature_grads[q] *
                  (1 + time_step/old_time_step)
                  -
                  old_old_temperature_grads[q] *
                  time_step/old_time_step)
                 :
                 old_temperature_grads[q]);

            const Tensor<1,dim> extrapolated_u
              = (use_bdf2_scheme ?
                 (old_velocity_values[q] *
                  (1 + time_step/old_time_step)
                  -
                  old_old_velocity_values[q] *
                  time_step/old_time_step)
                 :
                 old_velocity_values[q]);

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              local_rhs(i) += (T_term_for_rhs * phi_T[i]
                               -
                               time_step *
                               extrapolated_u * ext_grad_T * phi_T[i]
                               -
                               time_step *
                               nu * ext_grad_T * grad_phi_T[i]
                               +
                               time_step *
                               gamma_values[q] * phi_T[i])
                              *
                              temperature_fe_values.JxW(q);
          }

        cell->get_dof_indices (local_dof_indices);
        temperature_constraints.distribute_local_to_global (local_rhs,
                                                            local_dof_indices,
                                                            temperature_rhs);
      }
  }




  // @sect4{BoussinesqFlowProblem::solve}
  //
  // This function solves the linear systems of equations. Following the
  // introduction, we start with the Stokes system, where we need to generate
  // our block Schur preconditioner. Since all the relevant actions are
  // implemented in the class <code>BlockSchurPreconditioner</code>, all we
  // have to do is to initialize the class appropriately. What we need to pass
  // down is an <code>InverseMatrix</code> object for the pressure mass
  // matrix, which we set up using the respective class together with the IC
  // preconditioner we already generated, and the AMG preconditioner for the
  // velocity-velocity matrix. Note that both <code>Mp_preconditioner</code>
  // and <code>Amg_preconditioner</code> are only pointers, so we use
  // <code>*</code> to pass down the actual preconditioner objects.
  //
  // Once the preconditioner is ready, we create a GMRES solver for the block
  // system. Since we are working with Trilinos data structures, we have to
  // set the respective template argument in the solver. GMRES needs to
  // internally store temporary vectors for each iteration (see the discussion
  // in the results section of step-22) &ndash; the more vectors it can use,
  // the better it will generally perform. To keep memory demands in check, we
  // set the number of vectors to 100. This means that up to 100 solver
  // iterations, every temporary vector can be stored. If the solver needs to
  // iterate more often to get the specified tolerance, it will work on a
  // reduced set of vectors by restarting at every 100 iterations.
  //
  // With this all set up, we solve the system and distribute the constraints
  // in the Stokes system, i.e. hanging nodes and no-flux boundary condition,
  // in order to have the appropriate solution values even at constrained
  // dofs. Finally, we write the number of iterations to the screen.
  template <int dim>
  void BoussinesqFlowProblem<dim>::solve ()
  {
    std::cout << "   Solving..." << std::endl;

    {
      const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix,
            TrilinosWrappers::PreconditionIC>
            mp_inverse (stokes_preconditioner_matrix.block(1,1), *Mp_preconditioner);

      const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
            TrilinosWrappers::PreconditionIC>
            preconditioner (stokes_matrix, mp_inverse, *Amg_preconditioner);

      SolverControl solver_control (stokes_matrix.m(),
                                    1e-6*stokes_rhs.l2_norm());

      SolverGMRES<TrilinosWrappers::BlockVector>
      gmres (solver_control,
             SolverGMRES<TrilinosWrappers::BlockVector >::AdditionalData(100));

      for (unsigned int i=0; i<stokes_solution.size(); ++i)
        if (stokes_constraints.is_constrained(i))
          stokes_solution(i) = 0;

      gmres.solve(stokes_matrix, stokes_solution, stokes_rhs, preconditioner);

      stokes_constraints.distribute (stokes_solution);

      std::cout << "   "
                << solver_control.last_step()
                << " GMRES iterations for Stokes subsystem."
                << std::endl;
    }

    // Once we know the Stokes solution, we can determine the new time step
    // from the maximal velocity. We have to do this to satisfy the CFL
    // condition since convection terms are treated explicitly in the
    // temperature equation, as discussed in the introduction. The exact form
    // of the formula used here for the time step is discussed in the results
    // section of this program.
    //
    // There is a snatch here. The formula contains a division by the maximum
    // value of the velocity. However, at the start of the computation, we
    // have a constant temperature field (we start with a constant
    // temperature, and it will be nonconstant only after the first time step
    // during which the source acts). Constant temperature means that no
    // buoyancy acts, and so the velocity is zero. Dividing by it will not
    // likely lead to anything good.
    //
    // To avoid the resulting infinite time step, we ask whether the maximal
    // velocity is very small (in particular smaller than the values we
    // encounter during any of the following time steps) and if so rather than
    // dividing by zero we just divide by a small value, resulting in a large
    // but finite time step.
    old_time_step = time_step;
    const double maximal_velocity = get_maximal_velocity();

    if (maximal_velocity >= 0.01)
      time_step = 1./(1.7*dim*std::sqrt(1.*dim)) /
                  temperature_degree *
                  GridTools::minimal_cell_diameter(triangulation) /
                  maximal_velocity;
    else
      time_step = 1./(1.7*dim*std::sqrt(1.*dim)) /
                  temperature_degree *
                  GridTools::minimal_cell_diameter(triangulation) /
                  .01;

    std::cout << "   " << "Time step: " << time_step
              << std::endl;

    temperature_solution = old_temperature_solution;

    // Next we set up the temperature system and the right hand side using the
    // function <code>assemble_temperature_system()</code>.  Knowing the
    // matrix and right hand side of the temperature equation, we set up a
    // preconditioner and a solver. The temperature matrix is a mass matrix
    // (with eigenvalues around one) plus a Laplace matrix (with eigenvalues
    // between zero and $ch^{-2}$) times a small number proportional to the
    // time step $k_n$. Hence, the resulting symmetric and positive definite
    // matrix has eigenvalues in the range $[1,1+k_nh^{-2}]$ (up to
    // constants). This matrix is only moderately ill conditioned even for
    // small mesh sizes and we get a reasonably good preconditioner by simple
    // means, for example with an incomplete Cholesky decomposition
    // preconditioner (IC) as we also use for preconditioning the pressure
    // mass matrix solver. As a solver, we choose the conjugate gradient
    // method CG. As before, we tell the solver to use Trilinos vectors via
    // the template argument <code>TrilinosWrappers::Vector</code>.  Finally,
    // we solve, distribute the hanging node constraints and write out the
    // number of iterations.
    assemble_temperature_system (maximal_velocity);
    {

      SolverControl solver_control (temperature_matrix.m(),
                                    1e-8*temperature_rhs.l2_norm());
      SolverCG<TrilinosWrappers::Vector> cg (solver_control);

      TrilinosWrappers::PreconditionIC preconditioner;
      preconditioner.initialize (temperature_matrix);

      cg.solve (temperature_matrix, temperature_solution,
                temperature_rhs, preconditioner);

      temperature_constraints.distribute (temperature_solution);

      std::cout << "   "
                << solver_control.last_step()
                << " CG iterations for temperature."
                << std::endl;

      // At the end of this function, we step through the vector and read out
      // the maximum and minimum temperature value, which we also want to
      // output. This will come in handy when determining the correct constant
      // in the choice of time step as discuss in the results section of this
      // program.
      double min_temperature = temperature_solution(0),
             max_temperature = temperature_solution(0);
      for (unsigned int i=0; i<temperature_solution.size(); ++i)
        {
          min_temperature = std::min<double> (min_temperature,
                                              temperature_solution(i));
          max_temperature = std::max<double> (max_temperature,
                                              temperature_solution(i));
        }

      std::cout << "   Temperature range: "
                << min_temperature << ' ' << max_temperature
                << std::endl;
    }
  }



  // @sect4{BoussinesqFlowProblem::output_results}
  //
  // This function writes the solution to a VTK output file for visualization,
  // which is done every tenth time step. This is usually quite a simple task,
  // since the deal.II library provides functions that do almost all the job
  // for us. There is one new function compared to previous examples: We want
  // to visualize both the Stokes solution and the temperature as one data
  // set, but we have done all the calculations based on two different
  // DoFHandler objects. Luckily, the DataOut class is prepared to deal with
  // it. All we have to do is to not attach one single DoFHandler at the
  // beginning and then use that for all added vector, but specify the
  // DoFHandler to each vector separately. The rest is done as in step-22. We
  // create solution names (that are going to appear in the visualization
  // program for the individual components). The first <code>dim</code>
  // components are the vector velocity, and then we have pressure for the
  // Stokes part, whereas temperature is scalar. This information is read out
  // using the DataComponentInterpretation helper class. Next, we actually
  // attach the data vectors with their DoFHandler objects, build patches
  // according to the degree of freedom, which are (sub-) elements that
  // describe the data for visualization programs. Finally, we set a file name
  // (that includes the time step number) and write the vtk file.
  template <int dim>
  void BoussinesqFlowProblem<dim>::output_results ()  const
  {
    if (timestep_number % 10 != 0)
      return;

    std::vector<std::string> stokes_names (dim, "velocity");
    stokes_names.push_back ("p");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    stokes_component_interpretation
    (dim+1, DataComponentInterpretation::component_is_scalar);
    for (unsigned int i=0; i<dim; ++i)
      stokes_component_interpretation[i]
        = DataComponentInterpretation::component_is_part_of_vector;

    DataOut<dim> data_out;
    data_out.add_data_vector (stokes_dof_handler, stokes_solution,
                              stokes_names, stokes_component_interpretation);
    data_out.add_data_vector (temperature_dof_handler, temperature_solution,
                              "T");
    data_out.build_patches (std::min(stokes_degree, temperature_degree));

    std::ostringstream filename;
    filename << "solution-" << Utilities::int_to_string(timestep_number, 4) << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);
  }



  // @sect4{BoussinesqFlowProblem::refine_mesh}
  //
  // This function takes care of the adaptive mesh refinement. The three tasks
  // this function performs is to first find out which cells to
  // refine/coarsen, then to actually do the refinement and eventually
  // transfer the solution vectors between the two different grids. The first
  // task is simply achieved by using the well-established Kelly error
  // estimator on the temperature (it is the temperature we're mainly
  // interested in for this program, and we need to be accurate in regions of
  // high temperature gradients, also to not have too much numerical
  // diffusion). The second task is to actually do the remeshing. That
  // involves only basic functions as well, such as the
  // <code>refine_and_coarsen_fixed_fraction</code> that refines those cells
  // with the largest estimated error that together make up 80 per cent of the
  // error, and coarsens those cells with the smallest error that make up for
  // a combined 10 per cent of the error.
  //
  // If implemented like this, we would get a program that will not make much
  // progress: Remember that we expect temperature fields that are nearly
  // discontinuous (the diffusivity $\kappa$ is very small after all) and
  // consequently we can expect that a freely adapted mesh will refine further
  // and further into the areas of large gradients. This decrease in mesh size
  // will then be accompanied by a decrease in time step, requiring an
  // exceedingly large number of time steps to solve to a given final time. It
  // will also lead to meshes that are much better at resolving
  // discontinuities after several mesh refinement cycles than in the
  // beginning.
  //
  // In particular to prevent the decrease in time step size and the
  // correspondingly large number of time steps, we limit the maximal
  // refinement depth of the mesh. To this end, after the refinement indicator
  // has been applied to the cells, we simply loop over all cells on the
  // finest level and unselect them from refinement if they would result in
  // too high a mesh level.
  template <int dim>
  void BoussinesqFlowProblem<dim>::refine_mesh (const unsigned int max_grid_level)
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate (temperature_dof_handler,
                                        QGauss<dim-1>(temperature_degree+1),
                                        typename FunctionMap<dim>::type(),
                                        temperature_solution,
                                        estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
                                                       estimated_error_per_cell,
                                                       0.8, 0.1);
    if (triangulation.n_levels() > max_grid_level)
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active(max_grid_level);
           cell != triangulation.end(); ++cell)
        cell->clear_refine_flag ();

    // As part of mesh refinement we need to transfer the solution vectors
    // from the old mesh to the new one. To this end we use the
    // SolutionTransfer class and we have to prepare the solution vectors that
    // should be transferred to the new grid (we will lose the old grid once
    // we have done the refinement so the transfer has to happen concurrently
    // with refinement). What we definitely need are the current and the old
    // temperature (BDF-2 time stepping requires two old solutions). Since the
    // SolutionTransfer objects only support to transfer one object per dof
    // handler, we need to collect the two temperature solutions in one data
    // structure. Moreover, we choose to transfer the Stokes solution, too,
    // since we need the velocity at two previous time steps, of which only
    // one is calculated on the fly.
    //
    // Consequently, we initialize two SolutionTransfer objects for the Stokes
    // and temperature DoFHandler objects, by attaching them to the old dof
    // handlers. With this at place, we can prepare the triangulation and the
    // data vectors for refinement (in this order).
    std::vector<TrilinosWrappers::Vector> x_temperature (2);
    x_temperature[0] = temperature_solution;
    x_temperature[1] = old_temperature_solution;
    TrilinosWrappers::BlockVector x_stokes = stokes_solution;

    SolutionTransfer<dim,TrilinosWrappers::Vector>
    temperature_trans(temperature_dof_handler);
    SolutionTransfer<dim,TrilinosWrappers::BlockVector>
    stokes_trans(stokes_dof_handler);

    triangulation.prepare_coarsening_and_refinement();
    temperature_trans.prepare_for_coarsening_and_refinement(x_temperature);
    stokes_trans.prepare_for_coarsening_and_refinement(x_stokes);

    // Now everything is ready, so do the refinement and recreate the dof
    // structure on the new grid, and initialize the matrix structures and the
    // new vectors in the <code>setup_dofs</code> function. Next, we actually
    // perform the interpolation of the solutions between the grids. We create
    // another copy of temporary vectors for temperature (now corresponding to
    // the new grid), and let the interpolate function do the job. Then, the
    // resulting array of vectors is written into the respective vector member
    // variables. For the Stokes vector, everything is just the same &ndash;
    // except that we do not need another temporary vector since we just
    // interpolate a single vector. In the end, we have to tell the program
    // that the matrices and preconditioners need to be regenerated, since the
    // mesh has changed.
    triangulation.execute_coarsening_and_refinement ();
    setup_dofs ();

    std::vector<TrilinosWrappers::Vector> tmp (2);
    tmp[0].reinit (temperature_solution);
    tmp[1].reinit (temperature_solution);
    temperature_trans.interpolate(x_temperature, tmp);

    temperature_solution = tmp[0];
    old_temperature_solution = tmp[1];

    stokes_trans.interpolate (x_stokes, stokes_solution);

    rebuild_stokes_matrix         = true;
    rebuild_temperature_matrices  = true;
    rebuild_stokes_preconditioner = true;
  }



  // @sect4{BoussinesqFlowProblem::run}
  //
  // This function performs all the essential steps in the Boussinesq
  // program. It starts by setting up a grid (depending on the spatial
  // dimension, we choose some different level of initial refinement and
  // additional adaptive refinement steps, and then create a cube in
  // <code>dim</code> dimensions and set up the dofs for the first time. Since
  // we want to start the time stepping already with an adaptively refined
  // grid, we perform some pre-refinement steps, consisting of all assembly,
  // solution and refinement, but without actually advancing in time. Rather,
  // we use the vilified <code>goto</code> statement to jump out of the time
  // loop right after mesh refinement to start all over again on the new mesh
  // beginning at the <code>start_time_iteration</code> label.
  //
  // Before we start, we project the initial values to the grid and obtain the
  // first data for the <code>old_temperature_solution</code> vector. Then, we
  // initialize time step number and time step and start the time loop.
  template <int dim>
  void BoussinesqFlowProblem<dim>::run ()
  {
    const unsigned int initial_refinement = (dim == 2 ? 4 : 2);
    const unsigned int n_pre_refinement_steps = (dim == 2 ? 4 : 3);


    GridGenerator::hyper_cube (triangulation);
    global_Omega_diameter = GridTools::diameter (triangulation);

    triangulation.refine_global (initial_refinement);

    setup_dofs();

    unsigned int pre_refinement_step = 0;

start_time_iteration:

    VectorTools::project (temperature_dof_handler,
                          temperature_constraints,
                          QGauss<dim>(temperature_degree+2),
                          EquationData::TemperatureInitialValues<dim>(),
                          old_temperature_solution);

    timestep_number           = 0;
    time_step = old_time_step = 0;

    double time = 0;

    do
      {
        std::cout << "Timestep " << timestep_number
                  << ":  t=" << time
                  << std::endl;

        // The first steps in the time loop are all obvious &ndash; we
        // assemble the Stokes system, the preconditioner, the temperature
        // matrix (matrices and preconditioner do actually only change in case
        // we've remeshed before), and then do the solve. Before going on with
        // the next time step, we have to check whether we should first finish
        // the pre-refinement steps or if we should remesh (every fifth time
        // step), refining up to a level that is consistent with initial
        // refinement and pre-refinement steps. Last in the loop is to advance
        // the solutions, i.e. to copy the solutions to the next "older" time
        // level.
        assemble_stokes_system ();
        build_stokes_preconditioner ();
        assemble_temperature_matrix ();

        solve ();

        output_results ();

        std::cout << std::endl;

        if ((timestep_number == 0) &&
            (pre_refinement_step < n_pre_refinement_steps))
          {
            refine_mesh (initial_refinement + n_pre_refinement_steps);
            ++pre_refinement_step;
            goto start_time_iteration;
          }
        else if ((timestep_number > 0) && (timestep_number % 5 == 0))
          refine_mesh (initial_refinement + n_pre_refinement_steps);

        time += time_step;
        ++timestep_number;

        old_stokes_solution          = stokes_solution;
        old_old_temperature_solution = old_temperature_solution;
        old_temperature_solution     = temperature_solution;
      }
    // Do all the above until we arrive at time 100.
    while (time <= 100);
  }
}



// @sect3{The <code>main</code> function}
//
// The main function looks almost the same as in all other programs.
//
// There is one difference we have to be careful about. This program uses
// Trilinos and, typically, Trilinos is configured so that it can run in
// %parallel using MPI. This doesn't mean that it <i>has</i> to run in
// %parallel, and in fact this program (unlike step-32) makes no attempt at
// all to do anything in %parallel using MPI. Nevertheless, Trilinos wants the
// MPI system to be initialized. We do that be creating an object of type
// Utilities::MPI::MPI_InitFinalize that initializes MPI (if available) using
// the arguments given to main() (i.e., <code>argc</code> and
// <code>argv</code>) and de-initializes it again when the object goes out of
// scope.
int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step31;

      deallog.depth_console (0);

      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

      BoussinesqFlowProblem<2> flow_problem;
      flow_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
