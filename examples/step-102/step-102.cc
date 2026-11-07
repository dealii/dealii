/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2010 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Authors: Chih-Che Chueh, University of Victoria, 2010
 *          Wolfgang Bangerth, Texas A&M University, 2010
 */


// @sect3{Include files}

// The first step, as always, is to include the functionality of a number of
// deal.II and C++ header files.
//
// The list includes vector, matrix, and preconditioner wrappers for Trilinos,
// as well as the SUNDIALS IDA interface used for time integration.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/sundials/ida.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <numeric>


// At the end of this top-matter, we open a namespace for the current project
// into which all the following material will go, and then import all deal.II
// names into this namespace:
namespace Step102
{
  using namespace dealii;


  // @sect3{Boundary and initial value classes}

  // The following part is taken directly from step-21 so there is no need to
  // repeat the descriptions found there.
  template <int dim>
  class PressureBoundaryValues : public Function<dim>
  {
  public:
    PressureBoundaryValues()
      : Function<dim>(1)
    {}

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };


  template <int dim>
  double
  PressureBoundaryValues<dim>::value(const Point<dim> &p,
                                     const unsigned int /*component*/) const
  {
    return 1 - p[0];
  }


  template <int dim>
  class SaturationBoundaryValues : public Function<dim>
  {
  public:
    SaturationBoundaryValues()
      : Function<dim>(1)
    {}

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };



  template <int dim>
  double
  SaturationBoundaryValues<dim>::value(const Point<dim> &p,
                                       const unsigned int /*component*/) const
  {
    if (p[0] == 0)
      return 1;
    else
      return 0;
  }


  template <int dim>
  class SaturationInitialValues : public Function<dim>
  {
  public:
    SaturationInitialValues()
      : Function<dim>(1)
    {}

    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

    virtual void vector_value(const Point<dim> &p,
                              Vector<double>   &value) const override;
  };


  template <int dim>
  double
  SaturationInitialValues<dim>::value(const Point<dim> & /*p*/,
                                      const unsigned int /*component*/) const
  {
    return 0.2;
  }


  template <int dim>
  void SaturationInitialValues<dim>::vector_value(const Point<dim> &p,
                                                  Vector<double> &values) const
  {
    for (unsigned int c = 0; c < this->n_components; ++c)
      values(c) = SaturationInitialValues<dim>::value(p, c);
  }


  template <int dim>
  class InitialValues : public Function<dim>
  {
  public:
    InitialValues()
      : Function<dim>(dim + 2)
    {}

    virtual void vector_value(const Point<dim> &p,
                              Vector<double>   &values) const override;
  };


  template <int dim>
  void InitialValues<dim>::vector_value(const Point<dim> &p,
                                        Vector<double>   &values) const
  {
    AssertDimension(values.size(), dim + 2);

    values          = 0;
    values(dim + 1) = SaturationInitialValues<dim>().value(p);
  }


  // @sect3{Permeability models}

  // In this tutorial, we still use the two permeability models previously
  // used in step-21 so we again refrain from commenting in detail about them.
  namespace SingleCurvingCrack
  {
    template <int dim>
    class KInverse : public TensorFunction<2, dim>
    {
    public:
      KInverse()
        : TensorFunction<2, dim>()
      {}

      virtual void
      value_list(const std::vector<Point<dim>> &points,
                 std::vector<Tensor<2, dim>>   &values) const override;
    };


    template <int dim>
    void KInverse<dim>::value_list(const std::vector<Point<dim>> &points,
                                   std::vector<Tensor<2, dim>>   &values) const
    {
      AssertDimension(points.size(), values.size());

      for (unsigned int p = 0; p < points.size(); ++p)
        {
          values[p].clear();

          const double distance_to_flowline =
            std::fabs(points[p][1] - 0.5 - 0.1 * std::sin(10 * points[p][0]));

          const double permeability =
            std::max(std::exp(-(distance_to_flowline * distance_to_flowline) /
                              (0.1 * 0.1)),
                     0.01);

          for (unsigned int d = 0; d < dim; ++d)
            values[p][d][d] = 1. / permeability;
        }
    }
  } // namespace SingleCurvingCrack


  namespace RandomMedium
  {
    template <int dim>
    class KInverse : public TensorFunction<2, dim>
    {
    public:
      KInverse()
        : TensorFunction<2, dim>()
      {}

      virtual void
      value_list(const std::vector<Point<dim>> &points,
                 std::vector<Tensor<2, dim>>   &values) const override;

    private:
      static std::vector<Point<dim>> centers;
    };



    template <int dim>
    std::vector<Point<dim>> KInverse<dim>::centers = []() {
      const unsigned int N =
        (dim == 2 ? 40 : (dim == 3 ? 100 : throw ExcNotImplemented()));

      std::vector<Point<dim>> centers_list(N);
      for (unsigned int i = 0; i < N; ++i)
        for (unsigned int d = 0; d < dim; ++d)
          centers_list[i][d] = static_cast<double>(rand()) / RAND_MAX;

      return centers_list;
    }();



    template <int dim>
    void KInverse<dim>::value_list(const std::vector<Point<dim>> &points,
                                   std::vector<Tensor<2, dim>>   &values) const
    {
      AssertDimension(points.size(), values.size());

      for (unsigned int p = 0; p < points.size(); ++p)
        {
          values[p].clear();

          double permeability = 0;
          for (unsigned int i = 0; i < centers.size(); ++i)
            permeability +=
              std::exp(-(points[p] - centers[i]).norm_square() / (0.05 * 0.05));

          const double normalized_permeability =
            std::clamp(permeability, 0.01, 4.);

          for (unsigned int d = 0; d < dim; ++d)
            values[p][d][d] = 1. / normalized_permeability;
        }
    }
  } // namespace RandomMedium


  // @sect3{Physical quantities}

  // The implementations of all the physical quantities such as total mobility
  // $\lambda_t$ and fractional flow of water $F$ are taken from step-21 so
  // again we don't have do any comment about them. Compared to step-21 we
  // have added checks that the saturation passed to these functions is in
  // fact within the physically valid range. Furthermore, given that the
  // wetting phase moves at speed $\mathbf u F'(S)$ it is clear that $F'(S)$
  // must be greater or equal to zero, so we assert that as well to make sure
  // that our calculations to get at the formula for the derivative made
  // sense.
  double bounded_saturation(const double S)
  {
    Assert(std::isfinite(S), ExcMessage("Saturation is not finite."));
    return std::clamp(S, 0.0, 1.0);
  }


  double mobility_inverse(const double S, const double viscosity)
  {
    const double bounded_S = bounded_saturation(S);
    return 1.0 / (1.0 / viscosity * bounded_S * bounded_S +
                  (1 - bounded_S) * (1 - bounded_S));
  }


  double fractional_flow(const double S, const double viscosity)
  {
    const double bounded_S = bounded_saturation(S);

    return bounded_S * bounded_S /
           (bounded_S * bounded_S +
            viscosity * (1 - bounded_S) * (1 - bounded_S));
  }


  double fractional_flow_derivative(const double S, const double viscosity)
  {
    const double bounded_S = bounded_saturation(S);

    const double temp =
      (bounded_S * bounded_S + viscosity * (1 - bounded_S) * (1 - bounded_S));

    const double numerator =
      2.0 * bounded_S * temp -
      bounded_S * bounded_S *
        (2.0 * bounded_S - 2.0 * viscosity * (1 - bounded_S));
    const double denominator = Utilities::fixed_power<2>(temp);

    const double F_prime = numerator / denominator;

    Assert(F_prime >= 0, ExcInternalError());

    return F_prime;
  }


  // @sect3{Helper classes for solvers and preconditioners}

  // In this first part we define a number of classes that we need in the
  // construction of linear solvers and preconditioners. This part is
  // essentially the same as that used in step-31. The only difference is that
  // the original variable name stokes_matrix is replaced by another name
  // darcy_matrix to match our problem.
  namespace LinearSolvers
  {
    template <class MatrixType, class PreconditionerType>
    class InverseMatrix : public EnableObserverPointer
    {
    public:
      InverseMatrix(const MatrixType         &m,
                    const PreconditionerType &preconditioner);


      template <typename VectorType>
      void vmult(VectorType &dst, const VectorType &src) const;

    private:
      const ObserverPointer<const MatrixType> matrix;
      const PreconditionerType               &preconditioner;
    };


    template <class MatrixType, class PreconditionerType>
    InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix(
      const MatrixType         &m,
      const PreconditionerType &preconditioner)
      : matrix(&m)
      , preconditioner(preconditioner)
    {}



    template <class MatrixType, class PreconditionerType>
    template <typename VectorType>
    void InverseMatrix<MatrixType, PreconditionerType>::vmult(
      VectorType       &dst,
      const VectorType &src) const
    {
      SolverControl        solver_control(src.size(), 1e-7 * src.l2_norm());
      SolverCG<VectorType> cg(solver_control);

      dst = 0;

      try
        {
          cg.solve(*matrix, dst, src, preconditioner);
        }
      catch (std::exception &e)
        {
          Assert(false, ExcMessage(e.what()));
        }
    }

    template <class PreconditionerTypeA, class PreconditionerTypeMp>
    class BlockSchurPreconditioner : public EnableObserverPointer
    {
    public:
      BlockSchurPreconditioner(
        const TrilinosWrappers::BlockSparseMatrix &S,
        const InverseMatrix<TrilinosWrappers::SparseMatrix,
                            PreconditionerTypeMp> &Mpinv,
        const PreconditionerTypeA                 &Apreconditioner);

      void vmult(TrilinosWrappers::MPI::BlockVector       &dst,
                 const TrilinosWrappers::MPI::BlockVector &src) const;

    private:
      const ObserverPointer<const TrilinosWrappers::BlockSparseMatrix>
        darcy_matrix;
      const ObserverPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix,
                                                PreconditionerTypeMp>>
                                 m_inverse;
      const PreconditionerTypeA &a_preconditioner;

      mutable TrilinosWrappers::MPI::Vector tmp;
    };

    template <class PreconditionerTypeA,
              class PreconditionerTypeMp,
              class PreconditionerTypeS>
    class BlockDiagonalPreconditioner : public EnableObserverPointer
    {
    public:
      BlockDiagonalPreconditioner(
        const TrilinosWrappers::BlockSparseMatrix &S,
        const std::vector<IndexSet>               &darcy_partitioning,
        const InverseMatrix<TrilinosWrappers::SparseMatrix,
                            PreconditionerTypeMp> &Mpinv,
        const PreconditionerTypeA                 &Apreconditioner,
        const PreconditionerTypeS                 &Spreconditioner);

      void vmult(TrilinosWrappers::MPI::BlockVector       &dst,
                 const TrilinosWrappers::MPI::BlockVector &src) const;

    private:
      const ObserverPointer<const TrilinosWrappers::BlockSparseMatrix>
        system_matrix;
      const ObserverPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix,
                                                PreconditionerTypeMp>>
                                 m_inverse;
      const PreconditionerTypeA &a_preconditioner;
      const PreconditionerTypeS &s_preconditioner;

      mutable TrilinosWrappers::MPI::Vector      tmp;
      mutable TrilinosWrappers::MPI::BlockVector darcy_src;
      mutable TrilinosWrappers::MPI::BlockVector darcy_dst;
    };



    template <class PreconditionerTypeA, class PreconditionerTypeMp>
    BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::
      BlockSchurPreconditioner(
        const TrilinosWrappers::BlockSparseMatrix &S,
        const InverseMatrix<TrilinosWrappers::SparseMatrix,
                            PreconditionerTypeMp> &Mpinv,
        const PreconditionerTypeA                 &Apreconditioner)
      : darcy_matrix(&S)
      , m_inverse(&Mpinv)
      , a_preconditioner(Apreconditioner)
      , tmp(complete_index_set(darcy_matrix->block(1, 1).m()))
    {}


    template <class PreconditionerTypeA, class PreconditionerTypeMp>
    void
    BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::vmult(
      TrilinosWrappers::MPI::BlockVector       &dst,
      const TrilinosWrappers::MPI::BlockVector &src) const
    {
      a_preconditioner.vmult(dst.block(0), src.block(0));
      darcy_matrix->block(1, 0).residual(tmp, dst.block(0), src.block(1));
      tmp *= -1;
      m_inverse->vmult(dst.block(1), tmp);
    }


    template <class PreconditionerTypeA,
              class PreconditionerTypeMp,
              class PreconditionerTypeS>
    BlockDiagonalPreconditioner<PreconditionerTypeA,
                                PreconditionerTypeMp,
                                PreconditionerTypeS>::
      BlockDiagonalPreconditioner(
        const TrilinosWrappers::BlockSparseMatrix &S,
        const std::vector<IndexSet>               &darcy_partitioning,
        const InverseMatrix<TrilinosWrappers::SparseMatrix,
                            PreconditionerTypeMp> &Mpinv,
        const PreconditionerTypeA                 &Apreconditioner,
        const PreconditionerTypeS                 &Spreconditioner)
      : system_matrix(&S)
      , m_inverse(&Mpinv)
      , a_preconditioner(Apreconditioner)
      , s_preconditioner(Spreconditioner)
      , tmp(complete_index_set(system_matrix->block(1, 1).m()))
      , darcy_src(darcy_partitioning, MPI_COMM_WORLD)
      , darcy_dst(darcy_partitioning, MPI_COMM_WORLD)
    {}


    template <class PreconditionerTypeA,
              class PreconditionerTypeMp,
              class PreconditionerTypeS>
    void BlockDiagonalPreconditioner<PreconditionerTypeA,
                                     PreconditionerTypeMp,
                                     PreconditionerTypeS>::
      vmult(TrilinosWrappers::MPI::BlockVector       &dst,
            const TrilinosWrappers::MPI::BlockVector &src) const
    {
      darcy_src.block(0) = src.block(0);
      darcy_src.block(1) = src.block(1);
      darcy_dst          = 0;

      a_preconditioner.vmult(darcy_dst.block(0), darcy_src.block(0));
      system_matrix->block(1, 0).residual(tmp,
                                          darcy_dst.block(0),
                                          darcy_src.block(1));
      tmp *= -1;
      m_inverse->vmult(darcy_dst.block(1), tmp);

      dst          = 0;
      dst.block(0) = darcy_dst.block(0);
      dst.block(1) = darcy_dst.block(1);
      s_preconditioner.vmult(dst.block(2), src.block(2));
    }
  } // namespace LinearSolvers


  // @sect3{The TwoPhaseFlowProblem class}

  // This class implements a mixed finite element discretization of the
  // Buckley-Leverett problem with unknowns velocity, pressure, and saturation
  // on a single DoFHandler. Time integration is handled by SUNDIALS IDA:
  // velocity and pressure are treated as algebraic variables, while
  // saturation is treated as the differential component. The active code path
  // therefore centers around assembling the semidiscrete residual and an
  // approximate Jacobian for the coupled DAE system, together with the mesh
  // refinement and output callbacks needed by IDA.
  //
  // We keep a second constraint object for the Darcy preconditioner. It is
  // used to assemble the pressure block with homogeneous Dirichlet boundary
  // conditions so that the Schur-complement approximation remains positive
  // definite.
  template <int dim>
  class TwoPhaseFlowProblem
  {
  public:
    TwoPhaseFlowProblem(const unsigned int degree);
    void run();

  private:
    using TimeStepper = SUNDIALS::IDA<TrilinosWrappers::MPI::BlockVector>;

    static constexpr unsigned int velocity_block   = 0;
    static constexpr unsigned int pressure_block   = 1;
    static constexpr unsigned int saturation_block = 2;

    void setup_dofs();
    void setup_time_stepper();
    void assemble_darcy_preconditioner(
      const TrilinosWrappers::MPI::BlockVector &state);
    void
    build_darcy_preconditioner(const TrilinosWrappers::MPI::BlockVector &state);
    void assemble_darcy_system(const TrilinosWrappers::MPI::BlockVector &state);
    void assemble_saturation_matrix();
    void refine_mesh(const unsigned int                  min_grid_level,
                     const unsigned int                  max_grid_level,
                     TrilinosWrappers::MPI::BlockVector &y,
                     TrilinosWrappers::MPI::BlockVector &y_dot);
    void output_results(const TrilinosWrappers::MPI::BlockVector &state) const;
    void ida_residual(const double                              t,
                      const TrilinosWrappers::MPI::BlockVector &y,
                      const TrilinosWrappers::MPI::BlockVector &y_dot,
                      TrilinosWrappers::MPI::BlockVector       &residual);
    void ida_setup_jacobian(const double                              t,
                            const TrilinosWrappers::MPI::BlockVector &y,
                            const TrilinosWrappers::MPI::BlockVector &y_dot,
                            const double                              alpha);
    void ida_solve_with_jacobian(const TrilinosWrappers::MPI::BlockVector &rhs,
                                 TrilinosWrappers::MPI::BlockVector       &dst,
                                 const double tolerance);
    IndexSet ida_differential_components() const;
    bool     ida_solver_should_restart(const double                        t,
                                       TrilinosWrappers::MPI::BlockVector &sol,
                                       TrilinosWrappers::MPI::BlockVector &sol_dot);
    void     ida_output_step(const double                              t,
                             const TrilinosWrappers::MPI::BlockVector &sol,
                             const TrilinosWrappers::MPI::BlockVector &sol_dot,
                             const unsigned int                        step_number);

    // We follow with a number of helper functions that are used in a variety
    // of places throughout the program:
    double get_max_u_F_prime(const TrilinosWrappers::MPI::BlockVector &y) const;
    void   project_back_saturation(TrilinosWrappers::MPI::BlockVector &y);
    unsigned int component_to_block(const unsigned int component) const;
    unsigned int n_local_dofs_for_blocks(
      const std::vector<unsigned int> &selected_blocks) const;
    std::vector<unsigned int> local_dof_indices_for_blocks(
      const std::vector<unsigned int> &selected_blocks) const;
    types::global_dof_index get_block_offset(const unsigned int block) const;
    void                    extract_block_local_dof_indices(
                         const std::vector<types::global_dof_index> &local_dof_indices,
                         const std::vector<unsigned int>            &selected_blocks,
                         std::vector<types::global_dof_index>       &extracted_indices,
                         const bool                                  rebase_indices) const;
    double compute_viscosity(
      const std::vector<double>         &old_saturation,
      const std::vector<double>         &old_old_saturation,
      const std::vector<Tensor<1, dim>> &old_saturation_grads,
      const std::vector<Tensor<1, dim>> &old_old_saturation_grads,
      const std::vector<Tensor<1, dim>> &present_darcy_values,
      const double                       global_max_u_F_prime,
      const double                       global_S_variation,
      const double                       cell_diameter) const;
    // The remaining members describe the mesh, finite element spaces, linear
    // algebra objects, and IDA integration state.
    Triangulation<dim> triangulation;
    double             global_Omega_diameter;

    const unsigned int degree;

    const unsigned int        darcy_degree;
    const unsigned int        saturation_degree;
    const FESystem<dim>       fe;
    DoFHandler<dim>           dof_handler;
    AffineConstraints<double> constraints;

    AffineConstraints<double>            darcy_preconditioner_constraints;
    std::vector<types::global_dof_index> dofs_per_block;
    std::vector<IndexSet>                state_partitioning;
    std::vector<IndexSet>                darcy_partitioning;

    TrilinosWrappers::BlockSparseMatrix system_matrix;
    TrilinosWrappers::BlockSparseMatrix darcy_preconditioner_matrix;

    TrilinosWrappers::MPI::BlockVector darcy_rhs;

    TrilinosWrappers::MPI::BlockVector saturation_rhs;

    const double saturation_refinement_threshold;

    double       time;
    const double end_time;

    unsigned int timestep_number;

    const double       viscosity;
    const double       porosity;
    const double       refinement_interval;
    const unsigned int initial_refinement;
    const unsigned int max_refinement_level;

    std::shared_ptr<TrilinosWrappers::PreconditionIC> top_left_preconditioner;
    std::shared_ptr<TrilinosWrappers::PreconditionIC>
      bottom_right_preconditioner;
    std::shared_ptr<TrilinosWrappers::PreconditionIC> saturation_preconditioner;

    std::unique_ptr<TimeStepper> time_stepper;

    double next_refinement_time;

    // At the very end we declare a variable that denotes the material
    // model. Compared to step-21, we do this here as a member variable since
    // we will want to use it in a variety of places and so having a central
    // place where such a variable is declared will make it simpler to replace
    // one class by another (e.g. replace RandomMedium::KInverse by
    // SingleCurvingCrack::KInverse).
    const RandomMedium::KInverse<dim> k_inverse;
  };


  // @sect3{TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem}

  // We use a mixed element with Taylor-Hood velocity-pressure components and
  // one scalar saturation component. The constructor only stores the
  // configuration parameters; actual vectors, matrices, and the IDA object
  // are sized later in setup_dofs() and setup_time_stepper().
  template <int dim>
  TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem(const unsigned int degree)
    : triangulation(Triangulation<dim>::maximum_smoothing)
    , global_Omega_diameter(std::numeric_limits<double>::quiet_NaN())
    , degree(degree)
    , darcy_degree(degree)
    , saturation_degree(degree + 1)
    , fe(FESystem<dim>(FE_Q<dim>(darcy_degree + 1) ^ dim,
                       FE_Q<dim>(darcy_degree)),
         1,
         FE_Q<dim>(saturation_degree),
         1)
    , dof_handler(triangulation)

    , saturation_refinement_threshold(0.5)
    , time(0)
    , end_time(10)
    , timestep_number(0)
    , viscosity(0.2)
    , porosity(1.0)
    , refinement_interval(0.02)
    , initial_refinement(dim == 2 ? 5 : 2)
    , max_refinement_level(initial_refinement + (dim == 2 ? 3 : 2))
    , next_refinement_time(refinement_interval)
  {}


  // @sect3{TwoPhaseFlowProblem<dim>::setup_dofs}

  // Set up the single mixed DoFHandler, compute the block sizes for
  // `(u,p,S)`, build the full mixed system matrix used by the IDA Jacobian,
  // and create the auxiliary 2x2 Darcy preconditioner matrix. The full
  // problem uses one mixed constraint object; the auxiliary Darcy
  // preconditioner gets its own pressure-constrained view.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::setup_dofs()
  {
    std::vector<unsigned int> block_component(dim + 2, velocity_block);
    block_component[dim]     = pressure_block;
    block_component[dim + 1] = saturation_block;

    dof_handler.distribute_dofs(fe);
    DoFRenumbering::Cuthill_McKee(dof_handler);
    DoFRenumbering::component_wise(dof_handler, block_component);

    constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    constraints.close();


    {
      darcy_preconditioner_constraints.clear();

      const FEValuesExtractors::Scalar pressure(dim);

      DoFTools::make_hanging_node_constraints(dof_handler,
                                              darcy_preconditioner_constraints);
      DoFTools::make_zero_boundary_constraints(dof_handler,
                                               darcy_preconditioner_constraints,
                                               fe.component_mask(pressure));

      darcy_preconditioner_constraints.close();
    }


    const std::vector<types::global_dof_index> dofs_per_component =
      DoFTools::count_dofs_per_fe_component(dof_handler);

    dofs_per_block.resize(3);
    dofs_per_block[velocity_block] =
      std::accumulate(dofs_per_component.begin(),
                      dofs_per_component.begin() + dim,
                      types::global_dof_index(0));
    dofs_per_block[pressure_block]   = dofs_per_component[dim];
    dofs_per_block[saturation_block] = dofs_per_component[dim + 1];

    const types::global_dof_index n_u = dofs_per_block[velocity_block],
                                  n_p = dofs_per_block[pressure_block],
                                  n_s = dofs_per_block[saturation_block];

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << " (on " << triangulation.n_levels() << " levels)" << std::endl
              << "Number of degrees of freedom: " << n_u + n_p + n_s << " ("
              << n_u << '+' << n_p << '+' << n_s << ')' << std::endl
              << std::endl;

    {
      system_matrix.clear();
      BlockDynamicSparsityPattern dsp({n_u, n_p, n_s}, {n_u, n_p, n_s});
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
      system_matrix.reinit(dsp);
    }

    {
      top_left_preconditioner.reset();
      bottom_right_preconditioner.reset();
      saturation_preconditioner.reset();
      darcy_preconditioner_matrix.clear();

      DynamicSparsityPattern full_dsp(dof_handler.n_dofs(),
                                      dof_handler.n_dofs());

      Table<2, DoFTools::Coupling> coupling(dim + 2, dim + 2);
      for (unsigned int c = 0; c < dim + 2; ++c)
        for (unsigned int d = 0; d < dim + 2; ++d)
          if ((c != dim + 1) && (d != dim + 1) && (c == d))
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      DoFTools::make_sparsity_pattern(
        dof_handler, coupling, full_dsp, constraints, false);

      BlockDynamicSparsityPattern dsp({n_u, n_p}, {n_u, n_p});
      for (types::global_dof_index row = 0; row < n_u + n_p; ++row)
        for (auto entry = full_dsp.begin(row); entry != full_dsp.end(row);
             ++entry)
          if (entry->column() < n_u + n_p)
            dsp.add(row, entry->column());

      darcy_preconditioner_matrix.reinit(dsp);
    }
    state_partitioning = {complete_index_set(n_u),
                          complete_index_set(n_p),
                          complete_index_set(n_s)};

    darcy_partitioning = {complete_index_set(n_u), complete_index_set(n_p)};
    darcy_rhs.reinit(state_partitioning, MPI_COMM_WORLD);

    saturation_rhs.reinit(state_partitioning, MPI_COMM_WORLD);
  }


  template <int dim>
  void TwoPhaseFlowProblem<dim>::setup_time_stepper()
  {
    typename TimeStepper::AdditionalData data;
    data.initial_time      = 0.0;
    data.final_time        = end_time;
    data.initial_step_size = 1e-4;
    data.output_period     = .01;
    data.ic_type           = TimeStepper::AdditionalData::use_y_diff;
    data.reset_type        = TimeStepper::AdditionalData::none;

    std::cout << "IDA configuration: ic_type=" << data.ic_type
              << ", reset_type=" << data.reset_type << std::endl;

    time_stepper = std::make_unique<TimeStepper>(data, MPI_COMM_WORLD);

    time_stepper->reinit_vector = [&](TrilinosWrappers::MPI::BlockVector &v) {
      v.reinit(state_partitioning, MPI_COMM_WORLD);
    };

    time_stepper->residual =
      [&](const double                              t,
          const TrilinosWrappers::MPI::BlockVector &y,
          const TrilinosWrappers::MPI::BlockVector &y_dot,
          TrilinosWrappers::MPI::BlockVector       &residual) {
        ida_residual(t, y, y_dot, residual);
      };

    time_stepper->setup_jacobian =
      [&](const double                              t,
          const TrilinosWrappers::MPI::BlockVector &y,
          const TrilinosWrappers::MPI::BlockVector &y_dot,
          const double alpha) { ida_setup_jacobian(t, y, y_dot, alpha); };

    time_stepper->solve_with_jacobian =
      [&](const TrilinosWrappers::MPI::BlockVector &rhs,
          TrilinosWrappers::MPI::BlockVector       &dst,
          const double                              tolerance) {
        ida_solve_with_jacobian(rhs, dst, tolerance);
      };

    time_stepper->differential_components = [&]() {
      return ida_differential_components();
    };

    time_stepper->solver_should_restart =
      [&](const double                        t,
          TrilinosWrappers::MPI::BlockVector &sol,
          TrilinosWrappers::MPI::BlockVector &sol_dot) {
        return ida_solver_should_restart(t, sol, sol_dot);
      };

    time_stepper->output_step =
      [&](const double                              t,
          const TrilinosWrappers::MPI::BlockVector &sol,
          const TrilinosWrappers::MPI::BlockVector &sol_dot,
          const unsigned int                        step_number) {
        ida_output_step(t, sol, sol_dot, step_number);
      };
  }


  template <int dim>
  unsigned int TwoPhaseFlowProblem<dim>::component_to_block(
    const unsigned int component) const
  {
    if (component < dim)
      return velocity_block;
    if (component == dim)
      return pressure_block;

    AssertDimension(component, dim + 1);
    return saturation_block;
  }


  template <int dim>
  unsigned int TwoPhaseFlowProblem<dim>::n_local_dofs_for_blocks(
    const std::vector<unsigned int> &selected_blocks) const
  {
    return local_dof_indices_for_blocks(selected_blocks).size();
  }


  template <int dim>
  std::vector<unsigned int>
  TwoPhaseFlowProblem<dim>::local_dof_indices_for_blocks(
    const std::vector<unsigned int> &selected_blocks) const
  {
    std::vector<unsigned int> indices;

    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      {
        const unsigned int component = fe.system_to_component_index(i).first;
        const unsigned int block     = component_to_block(component);

        if (std::find(selected_blocks.begin(), selected_blocks.end(), block) !=
            selected_blocks.end())
          indices.push_back(i);
      }

    return indices;
  }


  template <int dim>
  types::global_dof_index
  TwoPhaseFlowProblem<dim>::get_block_offset(const unsigned int block) const
  {
    return std::accumulate(dofs_per_block.begin(),
                           dofs_per_block.begin() + block,
                           types::global_dof_index(0));
  }


  template <int dim>
  void TwoPhaseFlowProblem<dim>::extract_block_local_dof_indices(
    const std::vector<types::global_dof_index> &local_dof_indices,
    const std::vector<unsigned int>            &selected_blocks,
    std::vector<types::global_dof_index>       &extracted_indices,
    const bool                                  rebase_indices) const
  {
    extracted_indices.clear();

    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        const unsigned int component = fe.system_to_component_index(i).first;
        const unsigned int block     = component_to_block(component);

        if (std::find(selected_blocks.begin(), selected_blocks.end(), block) !=
            selected_blocks.end())
          extracted_indices.push_back(
            rebase_indices ? (local_dof_indices[i] - get_block_offset(block)) :
                             local_dof_indices[i]);
      }
  }


  // @sect3{Assembling matrices and preconditioners}

  // The next few functions are devoted to setting up the various system and
  // preconditioner matrices and right hand sides that we have to deal with in
  // this program.

  // @sect4{TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner}

  // This function assembles the matrix we use for preconditioning the Darcy
  // system. What we need are a vector @ref GlossMassMatrix "mass matrix" weighted by
  // $\left(\mathbf{K} \lambda_t\right)^{-1}$ on the velocity components and a
  // mass matrix weighted by $\left(\mathbf{K} \lambda_t\right)$ on the
  // pressure component. We start by generating a quadrature object of
  // appropriate order, the FEValues object that can give values and gradients
  // at the quadrature points (together with quadrature weights). Next we
  // create data structures for the cell matrix and the relation between local
  // and global DoFs. The vectors phi_u and grad_phi_p are going to hold the
  // values of the basis functions in order to faster build up the local
  // matrices, as was already done in step-22. Before we start the loop over
  // all active cells, we have to specify which components are pressure and
  // which are velocity.
  //
  // The creation of the local matrix is rather simple. There are only a term
  // weighted by $\left(\mathbf{K} \lambda_t\right)^{-1}$ (on the velocity)
  // and a Laplace matrix weighted by $\left(\mathbf{K} \lambda_t\right)$ to
  // be generated, so the creation of the local matrix is done in essentially
  // two lines. Since the material model functions at the top of this file
  // only provide the inverses of the permeability and mobility, we have to
  // compute $\mathbf K$ and $\lambda_t$ by hand from the given values, once
  // per quadrature point.
  //
  // Once the local matrix is ready (loop over rows and columns in the local
  // matrix on each quadrature point), we get the local DoF indices and write
  // the local information into the global matrix. We do this by directly
  // applying the constraints (i.e. darcy_preconditioner_constraints) that
  // takes care of hanging node and zero Dirichlet boundary condition
  // constraints. By doing so, we don't have to do that afterwards, and we
  // later don't have to use AffineConstraints::condense and
  // MatrixTools::apply_boundary_values, both functions that would need to
  // modify matrix and vector entries and so are difficult to write for the
  // Trilinos classes where we don't immediately have access to individual
  // memory locations.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner(
    const TrilinosWrappers::MPI::BlockVector &state)
  {
    std::cout << "   Rebuilding darcy preconditioner..." << std::endl;

    darcy_preconditioner_matrix = 0;

    const QGauss<dim> quadrature_formula(darcy_degree + 2);
    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_JxW_values | update_values |
                              update_gradients | update_quadrature_points);

    const std::vector<unsigned int> selected_local_dof_indices =
      local_dof_indices_for_blocks({velocity_block, pressure_block});
    const unsigned int dofs_per_cell = selected_local_dof_indices.size();
    const unsigned int n_q_points    = quadrature_formula.size();

    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);

    std::vector<double> old_saturation_values(n_q_points);

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_joint_dof_indices(
      fe.n_dofs_per_cell());
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
    std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);
    const FEValuesExtractors::Scalar saturation(dim + 1);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        local_matrix = 0;

        fe_values[saturation].get_function_values(state, old_saturation_values);

        k_inverse.value_list(fe_values.get_quadrature_points(),
                             k_inverse_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double old_s = old_saturation_values[q];

            const double inverse_mobility = mobility_inverse(old_s, viscosity);
            const double mobility         = 1.0 / inverse_mobility;
            const Tensor<2, dim> permeability = invert(k_inverse_values[q]);

            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                const unsigned int system_index = selected_local_dof_indices[k];
                phi_u[k]      = fe_values[velocities].value(system_index, q);
                grad_phi_p[k] = fe_values[pressure].gradient(system_index, q);
              }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  local_matrix(i, j) +=
                    (k_inverse_values[q] * inverse_mobility * phi_u[i] *
                       phi_u[j] +
                     permeability * mobility * grad_phi_p[i] * grad_phi_p[j]) *
                    fe_values.JxW(q);
                }
          }

        cell->get_dof_indices(local_joint_dof_indices);
        extract_block_local_dof_indices(local_joint_dof_indices,
                                        {velocity_block, pressure_block},
                                        local_dof_indices,
                                        false);
        darcy_preconditioner_constraints.distribute_local_to_global(
          local_matrix, local_dof_indices, darcy_preconditioner_matrix);
      }
  }


  // @sect4{TwoPhaseFlowProblem<dim>::build_darcy_preconditioner}

  // After calling the above functions to assemble the preconditioner matrix,
  // this function generates the inner preconditioners that are going to be
  // used for the Schur complement block preconditioner. The preconditioners
  // need to be regenerated at every saturation time step since they depend on
  // the saturation $S$ that varies with time.
  //
  // In here, we set up the preconditioner for the velocity-velocity matrix
  // $\mathbf{M}^{\mathbf{u}}$ and the Schur complement $\mathbf{S}$. As
  // explained in the introduction, we are going to use an IC preconditioner
  // based on the vector matrix $\mathbf{M}^{\mathbf{u}}$ and another based on
  // the scalar Laplace matrix $\tilde{\mathbf{S}}^p$ (which is spectrally
  // close to the Schur complement of the Darcy matrix). Usually, the
  // TrilinosWrappers::PreconditionIC class can be seen as a good black-box
  // preconditioner which does not need any special knowledge of the matrix
  // structure and/or the operator that's behind it.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::build_darcy_preconditioner(
    const TrilinosWrappers::MPI::BlockVector &state)
  {
    assemble_darcy_preconditioner(state);

    top_left_preconditioner =
      std::make_shared<TrilinosWrappers::PreconditionIC>();
    top_left_preconditioner->initialize(
      darcy_preconditioner_matrix.block(0, 0));

    bottom_right_preconditioner =
      std::make_shared<TrilinosWrappers::PreconditionIC>();
    bottom_right_preconditioner->initialize(
      darcy_preconditioner_matrix.block(1, 1));
  }


  // @sect4{TwoPhaseFlowProblem<dim>::assemble_darcy_system}

  // Assemble the Darcy part of the mixed residual and Jacobian. Even though
  // the full problem lives in one `(u,p,S)` system matrix, this function only
  // fills the `(u,p) x (u,p)` blocks and the corresponding right hand side.
  // Saturation enters here only through the mobility coefficient evaluated
  // from the current mixed state.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::assemble_darcy_system(
    const TrilinosWrappers::MPI::BlockVector &state)
  {
    system_matrix.block(velocity_block, velocity_block) = 0;
    system_matrix.block(velocity_block, pressure_block) = 0;
    system_matrix.block(pressure_block, velocity_block) = 0;
    system_matrix.block(pressure_block, pressure_block) = 0;
    darcy_rhs                                           = 0;

    const QGauss<dim>     quadrature_formula(darcy_degree + 2);
    const QGauss<dim - 1> face_quadrature_formula(darcy_degree + 2);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const Functions::ZeroFunction<dim> pressure_right_hand_side;
    const PressureBoundaryValues<dim>  pressure_boundary_values;

    std::vector<double>         pressure_rhs_values(n_q_points);
    std::vector<double>         boundary_values(n_face_q_points);
    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points);

    // Evaluate the saturation block at quadrature points to obtain the
    // mobility coefficient for the Darcy operator. We keep local arrays for
    // basis values and divergences to avoid repeatedly querying FEValues.
    std::vector<double> old_saturation_values(n_q_points);

    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
    std::vector<double>         div_phi_u(dofs_per_cell);
    std::vector<double>         phi_p(dofs_per_cell);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);
    const FEValuesExtractors::Scalar saturation(dim + 1);

    // The cell loop evaluates the mixed state on each cell, forms only the
    // velocity-pressure part of the local matrix, and inserts those
    // contributions into the full mixed matrix with the global mixed
    // numbering.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        local_matrix = 0;
        local_rhs    = 0;

        fe_values[saturation].get_function_values(state, old_saturation_values);

        pressure_right_hand_side.value_list(fe_values.get_quadrature_points(),
                                            pressure_rhs_values);
        k_inverse.value_list(fe_values.get_quadrature_points(),
                             k_inverse_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                phi_u[k]     = fe_values[velocities].value(k, q);
                div_phi_u[k] = fe_values[velocities].divergence(k, q);
                phi_p[k]     = fe_values[pressure].value(k, q);
              }
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                if (component_to_block(fe.system_to_component_index(i).first) ==
                    saturation_block)
                  continue;

                const double old_s = old_saturation_values[q];
                for (unsigned int j = 0; j <= i; ++j)
                  {
                    if (component_to_block(
                          fe.system_to_component_index(j).first) ==
                        saturation_block)
                      continue;

                    local_matrix(i, j) +=
                      (phi_u[i] * k_inverse_values[q] *
                         mobility_inverse(old_s, viscosity) * phi_u[j] -
                       div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                      fe_values.JxW(q);
                  }

                local_rhs(i) +=
                  (-phi_p[i] * pressure_rhs_values[q]) * fe_values.JxW(q);
              }
          }

        for (const auto &face : cell->face_iterators())
          if (face->at_boundary())
            {
              fe_face_values.reinit(cell, face);

              pressure_boundary_values.value_list(
                fe_face_values.get_quadrature_points(), boundary_values);

              for (unsigned int q = 0; q < n_face_q_points; ++q)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    if (component_to_block(
                          fe.system_to_component_index(i).first) ==
                        saturation_block)
                      continue;

                    const Tensor<1, dim> phi_i_u =
                      fe_face_values[velocities].value(i, q);

                    local_rhs(i) +=
                      -(phi_i_u * fe_face_values.normal_vector(q) *
                        boundary_values[q] * fe_face_values.JxW(q));
                  }
            }

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
            local_matrix(i, j) = local_matrix(j, i);

        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          local_matrix, local_rhs, local_dof_indices, system_matrix, darcy_rhs);
      }
  }


  // @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_matrix}

  // Assemble the saturation mass matrix into the `(S,S)` block of the full
  // mixed system matrix. This block is reused as part of the IDA Jacobian.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::assemble_saturation_matrix()
  {
    const QGauss<dim> quadrature_formula(saturation_degree + 2);

    FEValues<dim>                    fe_values(fe,
                            quadrature_formula,
                            update_values | update_JxW_values);
    const FEValuesExtractors::Scalar saturation(dim + 1);
    const unsigned int               dofs_per_cell = fe.n_dofs_per_cell();

    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        local_matrix = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            if (component_to_block(fe.system_to_component_index(i).first) ==
                saturation_block)
              {
                const double phi_i_s = fe_values[saturation].value(i, q);
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  if (component_to_block(
                        fe.system_to_component_index(j).first) ==
                      saturation_block)
                    {
                      const double phi_j_s = fe_values[saturation].value(j, q);
                      local_matrix(i, j) +=
                        porosity * phi_i_s * phi_j_s * fe_values.JxW(q);
                    }
              }
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(local_matrix,
                                               local_dof_indices,
                                               system_matrix);
      }
  }



  // @sect3{TwoPhaseFlowProblem<dim>::refine_mesh}

  // The next function does the refinement and coarsening of the mesh. It does
  // its work in three blocks: (i) Compute refinement indicators by looking at
  // the gradient of the current saturation solution. (ii) Flagging those
  // cells for refinement and coarsening where the gradient is larger or
  // smaller than a certain threshold, preserving minimal and maximal levels
  // of mesh refinement. (iii) Transferring the current IDA state and its time
  // derivative from the old to the new mesh. None of this is particularly
  // difficult.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::refine_mesh(
    const unsigned int                  min_grid_level,
    const unsigned int                  max_grid_level,
    TrilinosWrappers::MPI::BlockVector &y,
    TrilinosWrappers::MPI::BlockVector &y_dot)
  {
    Vector<double> refinement_indicators(triangulation.n_active_cells());
    const FEValuesExtractors::Scalar saturation(dim + 1);
    {
      const QMidpoint<dim> quadrature_formula;
      FEValues<dim>        fe_values(fe, quadrature_formula, update_gradients);
      std::vector<Tensor<1, dim>> grad_saturation(1);

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          const unsigned int cell_no = cell->active_cell_index();
          fe_values.reinit(cell);
          fe_values[saturation].get_function_gradients(y, grad_saturation);

          refinement_indicators(cell_no) = grad_saturation[0].norm();
        }
    }

    {
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          const unsigned int cell_no = cell->active_cell_index();
          cell->clear_coarsen_flag();
          cell->clear_refine_flag();

          if ((static_cast<unsigned int>(cell->level()) < max_grid_level) &&
              (std::fabs(refinement_indicators(cell_no)) >
               saturation_refinement_threshold))
            cell->set_refine_flag();
          else if ((static_cast<unsigned int>(cell->level()) >
                    min_grid_level) &&
                   (std::fabs(refinement_indicators(cell_no)) <
                    0.5 * saturation_refinement_threshold))
            cell->set_coarsen_flag();
        }
    }

    triangulation.prepare_coarsening_and_refinement();

    {
      std::vector<TrilinosWrappers::MPI::BlockVector> x_solution(2);
      x_solution[0] = y;
      x_solution[1] = y_dot;

      SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector>
        solution_soltrans(dof_handler);


      triangulation.prepare_coarsening_and_refinement();
      solution_soltrans.prepare_for_coarsening_and_refinement(x_solution);

      triangulation.execute_coarsening_and_refinement();
      setup_dofs();

      std::vector<TrilinosWrappers::MPI::BlockVector> tmp_solution(2);
      tmp_solution[0].reinit(state_partitioning, MPI_COMM_WORLD);
      tmp_solution[1].reinit(state_partitioning, MPI_COMM_WORLD);
      solution_soltrans.interpolate(tmp_solution);

      y     = tmp_solution[0];
      y_dot = tmp_solution[1];

      constraints.distribute(y);
      constraints.distribute(y_dot);
    }
  }



  // @sect3{TwoPhaseFlowProblem<dim>::output_results}

  // Write the current mixed state to VTU output.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::output_results(
    const TrilinosWrappers::MPI::BlockVector &state) const
  {
    Vector<double> joint_solution(dof_handler.n_dofs());
    for (unsigned int block = 0; block < state.n_blocks(); ++block)
      for (unsigned int i = 0; i < state.block(block).size(); ++i)
        joint_solution(get_block_offset(block) + i) = state.block(block)(i);

    std::vector<std::string> joint_solution_names(dim, "velocity");
    joint_solution_names.emplace_back("pressure");
    joint_solution_names.emplace_back("saturation");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(joint_solution,
                             joint_solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);

    data_out.build_patches();

    std::string filename =
      "solution-" + Utilities::int_to_string(timestep_number, 5) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
  }


  template <int dim>
  void TwoPhaseFlowProblem<dim>::ida_residual(
    const double                              t,
    const TrilinosWrappers::MPI::BlockVector &y,
    const TrilinosWrappers::MPI::BlockVector &y_dot,
    TrilinosWrappers::MPI::BlockVector       &residual)
  {
    static unsigned int residual_call = 0;
    ++residual_call;
    std::cout << "IDA residual call " << residual_call << ": t=" << t
              << std::endl;

    residual = 0;

    TrilinosWrappers::MPI::BlockVector constrained_y(y);
    TrilinosWrappers::MPI::BlockVector constrained_y_dot(y_dot);
    constraints.distribute(constrained_y);
    constraints.distribute(constrained_y_dot);

    assemble_darcy_system(constrained_y);

    TrilinosWrappers::MPI::BlockVector darcy_state;
    darcy_state.reinit({complete_index_set(dofs_per_block[velocity_block]),
                        complete_index_set(dofs_per_block[pressure_block])},
                       MPI_COMM_WORLD);
    darcy_state.block(0) = constrained_y.block(velocity_block);
    darcy_state.block(1) = constrained_y.block(pressure_block);

    TrilinosWrappers::MPI::BlockVector darcy_residual;
    darcy_residual.reinit(darcy_rhs);
    darcy_residual = 0;
    system_matrix.block(velocity_block, velocity_block)
      .vmult(darcy_residual.block(velocity_block), darcy_state.block(0));
    system_matrix.block(velocity_block, pressure_block)
      .vmult_add(darcy_residual.block(velocity_block), darcy_state.block(1));
    system_matrix.block(pressure_block, velocity_block)
      .vmult(darcy_residual.block(pressure_block), darcy_state.block(0));
    system_matrix.block(pressure_block, pressure_block)
      .vmult_add(darcy_residual.block(pressure_block), darcy_state.block(1));
    darcy_residual -= darcy_rhs;

    residual.block(velocity_block) = darcy_residual.block(velocity_block);
    residual.block(pressure_block) = darcy_residual.block(pressure_block);

    const QGauss<dim>     quadrature_formula(saturation_degree + 2);
    const QGauss<dim - 1> face_quadrature_formula(saturation_degree + 2);

    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar saturation(dim + 1);

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    Vector<double>     local_residual(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> velocity_values(n_q_points);
    std::vector<double>         saturation_values(n_q_points);
    std::vector<double>         saturation_dot_values(n_q_points);
    std::vector<Tensor<1, dim>> saturation_gradients(n_q_points);

    std::vector<Tensor<1, dim>> velocity_values_face(n_face_q_points);
    std::vector<double>         saturation_values_face(n_face_q_points);
    std::vector<double>         boundary_saturation(n_face_q_points);

    SaturationBoundaryValues<dim> saturation_boundary_values;

    double global_max_u_F_prime = 0.0;
    double min_saturation       = std::numeric_limits<double>::max();
    double max_saturation       = -std::numeric_limits<double>::max();

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        fe_values[velocities].get_function_values(constrained_y,
                                                  velocity_values);
        fe_values[saturation].get_function_values(constrained_y,
                                                  saturation_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double bounded_s = bounded_saturation(saturation_values[q]);
            min_saturation         = std::min(min_saturation, bounded_s);
            max_saturation         = std::max(max_saturation, bounded_s);
            global_max_u_F_prime =
              std::max(global_max_u_F_prime,
                       velocity_values[q].norm() *
                         fractional_flow_derivative(bounded_s, viscosity));
          }
      }

    const double global_S_variation =
      std::max(max_saturation - min_saturation, 1e-12);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        local_residual = 0;

        fe_values[velocities].get_function_values(constrained_y,
                                                  velocity_values);
        fe_values[saturation].get_function_values(constrained_y,
                                                  saturation_values);
        fe_values[saturation].get_function_values(constrained_y_dot,
                                                  saturation_dot_values);
        fe_values[saturation].get_function_gradients(constrained_y,
                                                     saturation_gradients);

        const double nu = compute_viscosity(saturation_values,
                                            saturation_values,
                                            saturation_gradients,
                                            saturation_gradients,
                                            velocity_values,
                                            global_max_u_F_prime,
                                            global_S_variation,
                                            cell->diameter());

        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            if (component_to_block(fe.system_to_component_index(i).first) ==
                saturation_block)
              {
                const double phi_i_s = fe_values[saturation].value(i, q);
                const Tensor<1, dim> grad_phi_i_s =
                  fe_values[saturation].gradient(i, q);

                local_residual(i) +=
                  (porosity * saturation_dot_values[q] * phi_i_s -
                   fractional_flow(saturation_values[q], viscosity) *
                     velocity_values[q] * grad_phi_i_s +
                   nu * saturation_gradients[q] * grad_phi_i_s) *
                  fe_values.JxW(q);
              }

        for (const auto &face : cell->face_iterators())
          if (face->at_boundary())
            {
              fe_face_values.reinit(cell, face);

              fe_face_values[velocities].get_function_values(
                constrained_y, velocity_values_face);
              fe_face_values[saturation].get_function_values(
                constrained_y, saturation_values_face);
              saturation_boundary_values.value_list(
                fe_face_values.get_quadrature_points(), boundary_saturation);

              for (unsigned int q = 0; q < n_face_q_points; ++q)
                {
                  const double normal_flux =
                    velocity_values_face[q] * fe_face_values.normal_vector(q);
                  const bool   is_outflow = (normal_flux >= 0.0);
                  const double upwind_saturation =
                    (is_outflow ? saturation_values_face[q] :
                                  boundary_saturation[q]);

                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    if (component_to_block(
                          fe.system_to_component_index(i).first) ==
                        saturation_block)
                      local_residual(i) +=
                        normal_flux *
                        fractional_flow(upwind_saturation, viscosity) *
                        fe_face_values[saturation].value(i, q) *
                        fe_face_values.JxW(q);
                }
            }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(local_residual,
                                               local_dof_indices,
                                               residual);
      }
  }


  template <int dim>
  void TwoPhaseFlowProblem<dim>::ida_setup_jacobian(
    const double                              t,
    const TrilinosWrappers::MPI::BlockVector &y,
    const TrilinosWrappers::MPI::BlockVector &y_dot,
    const double                              alpha)
  {
    static unsigned int jacobian_call = 0;
    ++jacobian_call;
    std::cout << "IDA jacobian setup " << jacobian_call << ": t=" << t
              << ", alpha=" << alpha << std::endl;

    TrilinosWrappers::MPI::BlockVector constrained_y(y);
    constraints.distribute(constrained_y);

    assemble_darcy_system(constrained_y);
    build_darcy_preconditioner(constrained_y);

    system_matrix.block(saturation_block, saturation_block) = 0;

    const QGauss<dim>     quadrature_formula(saturation_degree + 2);
    const QGauss<dim - 1> face_quadrature_formula(saturation_degree + 2);

    FEValues<dim>     fe_values(fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values(fe,
                                     face_quadrature_formula,
                                     update_values | update_normal_vectors |
                                       update_quadrature_points |
                                       update_JxW_values);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar saturation(dim + 1);

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell();
    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> velocity_values(n_q_points);
    std::vector<double>         saturation_values(n_q_points);
    std::vector<Tensor<1, dim>> saturation_gradients(n_q_points);
    std::vector<double>         saturation_dot_values(n_q_points);
    std::vector<Tensor<1, dim>> velocity_values_face(n_face_q_points);
    std::vector<double>         saturation_values_face(n_face_q_points);

    double global_max_u_F_prime = 0.0;
    double min_saturation       = std::numeric_limits<double>::max();
    double max_saturation       = -std::numeric_limits<double>::max();

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        fe_values[velocities].get_function_values(constrained_y,
                                                  velocity_values);
        fe_values[saturation].get_function_values(constrained_y,
                                                  saturation_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double bounded_s = bounded_saturation(saturation_values[q]);
            min_saturation         = std::min(min_saturation, bounded_s);
            max_saturation         = std::max(max_saturation, bounded_s);
            global_max_u_F_prime =
              std::max(global_max_u_F_prime,
                       velocity_values[q].norm() *
                         fractional_flow_derivative(bounded_s, viscosity));
          }
      }

    const double global_S_variation =
      std::max(max_saturation - min_saturation, 1e-12);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);
        local_matrix = 0;

        fe_values[velocities].get_function_values(constrained_y,
                                                  velocity_values);
        fe_values[saturation].get_function_values(constrained_y,
                                                  saturation_values);
        fe_values[saturation].get_function_gradients(constrained_y,
                                                     saturation_gradients);
        fe_values[saturation].get_function_values(y_dot, saturation_dot_values);

        const double nu = compute_viscosity(saturation_values,
                                            saturation_values,
                                            saturation_gradients,
                                            saturation_gradients,
                                            velocity_values,
                                            global_max_u_F_prime,
                                            global_S_variation,
                                            cell->diameter());

        for (unsigned int q = 0; q < n_q_points; ++q)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            if (component_to_block(fe.system_to_component_index(i).first) ==
                saturation_block)
              {
                const double phi_i_s = fe_values[saturation].value(i, q);
                const Tensor<1, dim> grad_phi_i_s =
                  fe_values[saturation].gradient(i, q);

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  if (component_to_block(
                        fe.system_to_component_index(j).first) ==
                      saturation_block)
                    {
                      const double phi_j_s = fe_values[saturation].value(j, q);
                      local_matrix(i, j) +=
                        (alpha * porosity * phi_i_s * phi_j_s -
                         fractional_flow_derivative(saturation_values[q],
                                                    viscosity) *
                           phi_j_s * velocity_values[q] * grad_phi_i_s +
                         nu * fe_values[saturation].gradient(j, q) *
                           grad_phi_i_s) *
                        fe_values.JxW(q);
                    }
              }

        for (const auto &face : cell->face_iterators())
          if (face->at_boundary())
            {
              fe_face_values.reinit(cell, face);
              fe_face_values[velocities].get_function_values(
                constrained_y, velocity_values_face);
              fe_face_values[saturation].get_function_values(
                constrained_y, saturation_values_face);

              for (unsigned int q = 0; q < n_face_q_points; ++q)
                {
                  const double normal_flux =
                    velocity_values_face[q] * fe_face_values.normal_vector(q);

                  if (normal_flux >= 0.0)
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      if (component_to_block(
                            fe.system_to_component_index(i).first) ==
                          saturation_block)
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                          if (component_to_block(
                                fe.system_to_component_index(j).first) ==
                              saturation_block)
                            local_matrix(i, j) +=
                              normal_flux *
                              fractional_flow_derivative(
                                saturation_values_face[q], viscosity) *
                              fe_face_values[saturation].value(j, q) *
                              fe_face_values[saturation].value(i, q) *
                              fe_face_values.JxW(q);
                }
            }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(local_matrix,
                                               local_dof_indices,
                                               system_matrix);
      }

    saturation_preconditioner =
      std::make_shared<TrilinosWrappers::PreconditionIC>();
    saturation_preconditioner->initialize(
      system_matrix.block(saturation_block, saturation_block));
  }


  template <int dim>
  void TwoPhaseFlowProblem<dim>::ida_solve_with_jacobian(
    const TrilinosWrappers::MPI::BlockVector &rhs,
    TrilinosWrappers::MPI::BlockVector       &dst,
    const double                              tolerance)
  {
    static unsigned int linear_solve_call = 0;
    ++linear_solve_call;
    std::cout << "IDA linear solve " << linear_solve_call
              << ": tolerance=" << tolerance << ", |rhs|=" << rhs.l2_norm()
              << std::endl;

    dst = 0;

    const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix,
                                       TrilinosWrappers::PreconditionIC>
      mp_inverse(darcy_preconditioner_matrix.block(1, 1),
                 *bottom_right_preconditioner);

    const LinearSolvers::BlockDiagonalPreconditioner<
      TrilinosWrappers::PreconditionIC,
      TrilinosWrappers::PreconditionIC,
      TrilinosWrappers::PreconditionIC>
      preconditioner(system_matrix,
                     darcy_partitioning,
                     mp_inverse,
                     *top_left_preconditioner,
                     *saturation_preconditioner);

    SolverControl solver_control(system_matrix.m(),
                                 std::max(tolerance, 1e-12 * rhs.l2_norm()));
    SolverGMRES<TrilinosWrappers::MPI::BlockVector> gmres(
      solver_control,
      SolverGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(100));
    gmres.solve(system_matrix, dst, rhs, preconditioner);
  }


  template <int dim>
  IndexSet TwoPhaseFlowProblem<dim>::ida_differential_components() const
  {
    IndexSet differential_components(dof_handler.n_dofs());
    differential_components.add_range(get_block_offset(saturation_block),
                                      get_block_offset(saturation_block) +
                                        dofs_per_block[saturation_block]);
    return differential_components;
  }


  template <int dim>
  bool TwoPhaseFlowProblem<dim>::ida_solver_should_restart(
    const double                        t,
    TrilinosWrappers::MPI::BlockVector &sol,
    TrilinosWrappers::MPI::BlockVector &sol_dot)
  {
    if (t + 1e-12 < next_refinement_time)
      return false;

    std::cout << "IDA restart for mesh refinement at t=" << t << std::endl;
    std::cout << "  restart input blocks: y=" << sol.n_blocks()
              << ", y_dot=" << sol_dot.n_blocks() << std::endl;

    time = t;
    constraints.distribute(sol);
    constraints.distribute(sol_dot);

    refine_mesh(initial_refinement, max_refinement_level, sol, sol_dot);

    std::cout << "  restart output blocks: y=" << sol.n_blocks()
              << ", y_dot=" << sol_dot.n_blocks() << std::endl;

    while (next_refinement_time <= t + 1e-12)
      next_refinement_time += refinement_interval;

    return true;
  }


  template <int dim>
  void TwoPhaseFlowProblem<dim>::ida_output_step(
    const double                              t,
    const TrilinosWrappers::MPI::BlockVector &sol,
    const TrilinosWrappers::MPI::BlockVector &sol_dot,
    const unsigned int                        step_number)
  {
    std::cout << "IDA output step " << step_number << ": t=" << t
              << ", |y|=" << sol.l2_norm() << ", |y_dot|=" << sol_dot.l2_norm()
              << std::endl;

    time            = t;
    timestep_number = step_number;

    TrilinosWrappers::MPI::BlockVector output_state(sol);
    constraints.distribute(output_state);
    output_results(output_state);
  }



  // @sect3{Tool functions}

  // @sect4{TwoPhaseFlowProblem<dim>::project_back_saturation}

  // Clamp the saturation block to the physically meaningful interval `[0,1]`.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::project_back_saturation(
    TrilinosWrappers::MPI::BlockVector &y)
  {
    for (unsigned int i = 0; i < y.block(saturation_block).size(); ++i)
      if (y.block(saturation_block)(i) < 0.0)
        y.block(saturation_block)(i) = 0.0;
      else if (y.block(saturation_block)(i) > 1)
        y.block(saturation_block)(i) = 1;
  }



  // @sect4{TwoPhaseFlowProblem<dim>::get_max_u_F_prime}
  //
  // Compute an estimate of `||u F'(S)||_{L_infty}` used to scale the
  // artificial viscosity.
  template <int dim>
  double TwoPhaseFlowProblem<dim>::get_max_u_F_prime(
    const TrilinosWrappers::MPI::BlockVector &y) const
  {
    const QGauss<dim>  quadrature_formula(darcy_degree + 2);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values(fe, quadrature_formula, update_values);
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar saturation(dim + 1);

    std::vector<Tensor<1, dim>> darcy_solution_values(n_q_points);
    std::vector<double>         saturation_values(n_q_points);

    double max_velocity_times_dF_dS = 0;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        fe_values[velocities].get_function_values(y, darcy_solution_values);
        fe_values[saturation].get_function_values(y, saturation_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double dF_dS =
              fractional_flow_derivative(saturation_values[q], viscosity);

            max_velocity_times_dF_dS =
              std::max(max_velocity_times_dF_dS,
                       darcy_solution_values[q].norm() * dF_dS);
          }
      }

    return max_velocity_times_dF_dS;
  }


  // @sect4{TwoPhaseFlowProblem<dim>::compute_viscosity}
  //
  // Compute the cellwise artificial viscosity used in the saturation
  // residual. In the current IDA formulation this is a lagged quantity built
  // from the current state and used as stabilization for the transport part.
  template <int dim>
  double TwoPhaseFlowProblem<dim>::compute_viscosity(
    const std::vector<double>         &old_saturation,
    const std::vector<double>         &old_old_saturation,
    const std::vector<Tensor<1, dim>> &old_saturation_grads,
    const std::vector<Tensor<1, dim>> &old_old_saturation_grads,
    const std::vector<Tensor<1, dim>> &present_darcy_values,
    const double                       global_max_u_F_prime,
    const double                       global_S_variation,
    const double                       cell_diameter) const
  {
    const double beta  = .4 * dim;
    const double alpha = 1;

    if (global_max_u_F_prime == 0)
      return 5e-3 * cell_diameter;

    const unsigned int n_q_points = old_saturation.size();

    double max_residual             = 0;
    double max_velocity_times_dF_dS = 0;

    const bool use_dF_dS = true;

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        const Tensor<1, dim> u     = present_darcy_values[q];
        const double         dS_dt = 0.0;

        const double dF_dS = fractional_flow_derivative(
          (old_saturation[q] + old_old_saturation[q]) / 2.0, viscosity);

        const double u_grad_S =
          u * dF_dS * (old_saturation_grads[q] + old_old_saturation_grads[q]) /
          2.0;

        const double residual =
          std::abs((dS_dt + u_grad_S) *
                   std::pow((old_saturation[q] + old_old_saturation[q]) / 2,
                            alpha - 1.));

        max_residual = std::max(residual, max_residual);
        max_velocity_times_dF_dS =
          std::max(std::sqrt(u * u) * (use_dF_dS ? std::max(dF_dS, 1.) : 1),
                   max_velocity_times_dF_dS);
      }

    const double c_R            = 1.0;
    const double global_scaling = c_R * porosity *
                                  (global_max_u_F_prime)*global_S_variation /
                                  std::pow(global_Omega_diameter, alpha - 2.);

    return (beta *
            (max_velocity_times_dF_dS)*std::min(cell_diameter,
                                                std::pow(cell_diameter, alpha) *
                                                  max_residual /
                                                  global_scaling));
  }


  // @sect3{TwoPhaseFlowProblem<dim>::run}

  // Set up the mesh and mixed finite element spaces, initialize the
  // `(u,p,S)` state, and hand the problem over to IDA. Output and adaptive
  // refinement are driven by the IDA callbacks configured in
  // `setup_time_stepper()`.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::run()
  {
    TrilinosWrappers::MPI::BlockVector y;
    TrilinosWrappers::MPI::BlockVector y_dot;

    GridGenerator::hyper_cube(triangulation, 0, 1);
    triangulation.refine_global(initial_refinement);
    global_Omega_diameter = GridTools::diameter(triangulation);

    setup_dofs();
    y.reinit(state_partitioning, MPI_COMM_WORLD);
    y_dot.reinit(state_partitioning, MPI_COMM_WORLD);
    setup_time_stepper();

    VectorTools::project(dof_handler,
                         constraints,
                         QGauss<dim>(saturation_degree + 2),
                         InitialValues<dim>(),
                         y);
    constraints.distribute(y);

    y_dot = 0;

    time                 = 0;
    timestep_number      = 0;
    next_refinement_time = refinement_interval;

    output_results(y);
    time_stepper->solve_dae(y, y_dot);
  }
} // namespace Step102



// @sect3{The <code>main()</code> function}
//
// The main function only initializes MPI and runs the serial example.
int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step102;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);

      // This program can only be run in serial. Otherwise, throw an exception.
      AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1,
                  ExcMessage(
                    "This program can only be run in serial, use ./step-102"));

      TwoPhaseFlowProblem<2> two_phase_flow_problem(1);
      two_phase_flow_problem.run();
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
