/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2019 - 2025 by the deal.II authors
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
 * Authors: Bruno Turcksin, Daniel Arndt, Oak Ridge National Laboratory, 2019
 */

// First include the necessary files from the deal.II library known from the
// previous tutorials.
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// The following ones include the data structures for the
// implementation of matrix-free methods on GPU:
#include <deal.II/matrix_free/portable_fe_evaluation.h>
#include <deal.II/matrix_free/portable_matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <fstream>


// As usual, we enclose everything into a namespace of its own:
namespace Step64
{
  using namespace dealii;


  // @sect3{Class <code>VaryingCoefficientFunctor</code>}

  // Next, we define a class that implements the varying coefficients
  // we want to use in the Helmholtz operator. Later, we want to pass
  // an object of this type to a Portable::MatrixFree
  // object that expects the class to have an `operator()` that fills the
  // values provided in the constructor for a given cell. This operator
  // needs to run on the device, so it needs to be marked as
  // `DEAL_II_HOST_DEVICE` for the compiler.
  template <int dim, int fe_degree>
  class VaryingCoefficientFunctor
  {
  public:
    VaryingCoefficientFunctor(double *coefficient)
      : coef(coefficient)
    {}

    DEAL_II_HOST_DEVICE void
    operator()(const typename Portable::MatrixFree<dim, double>::Data *gpu_data,
               const unsigned int                                      cell,
               const unsigned int                                      q) const;

    // Since Portable::MatrixFree::Data doesn't know about the size of its
    // arrays, we need to store the number of quadrature points and the
    // number of degrees of freedom in this class to do necessary index
    // conversions.
    static const unsigned int n_local_dofs = Utilities::pow(fe_degree + 1, dim);
    static const unsigned int n_q_points   = Utilities::pow(fe_degree + 1, dim);

  private:
    double *coef;
  };



  // The following function implements this coefficient. Recall from
  // the introduction that we have defined it as $a(\mathbf
  // x)=\frac{10}{0.05 + 2\|\mathbf x\|^2}$
  template <int dim, int fe_degree>
  DEAL_II_HOST_DEVICE void
  VaryingCoefficientFunctor<dim, fe_degree>::operator()(
    const typename Portable::MatrixFree<dim, double>::Data *gpu_data,
    const unsigned int                                      cell,
    const unsigned int                                      q) const
  {
    const unsigned int pos = gpu_data->local_q_point_id(cell, n_q_points, q);
    const Point<dim>   q_point = gpu_data->get_quadrature_point(cell, q);

    double p_square = 0.;
    for (unsigned int i = 0; i < dim; ++i)
      {
        const double coord = q_point[i];
        p_square += coord * coord;
      }
    coef[pos] = 10. / (0.05 + 2. * p_square);
  }


  // @sect3{Class <code>HelmholtzOperatorQuad</code>}

  // The class `HelmholtzOperatorQuad` implements the evaluation of
  // the Helmholtz operator at each quadrature point. It uses a
  // similar mechanism as the MatrixFree framework introduced in
  // step-37. In contrast to there, the actual quadrature point
  // index is treated implicitly by converting the current thread
  // index. As before, the functions of this class need to run on
  // the device, so need to be marked as `DEAL_II_HOST_DEVICE` for the
  // compiler.
  template <int dim, int fe_degree>
  class HelmholtzOperatorQuad
  {
  public:
    DEAL_II_HOST_DEVICE HelmholtzOperatorQuad(double *coef)
      : coef(coef)
    {}

    DEAL_II_HOST_DEVICE void operator()(
      Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double> *fe_eval,
      const int q_point) const;



    static const unsigned int n_q_points =
      dealii::Utilities::pow(fe_degree + 1, dim);

    static const unsigned int n_local_dofs = n_q_points;

  private:
    double *coef;
  };


  // The Helmholtz problem we want to solve here reads in weak form as follows:
  // @f{eqnarray*}{
  //   (\nabla v, \nabla u)+ (v, a(\mathbf x) u) &=&(v,1) \quad \forall v.
  // @f}
  // If you have seen step-37, then it will be obvious that
  // the two terms on the left-hand side correspond to the two function calls
  // here:
  template <int dim, int fe_degree>
  DEAL_II_HOST_DEVICE void HelmholtzOperatorQuad<dim, fe_degree>::operator()(
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double> *fe_eval,
    const int q_point) const
  {
    const int cell_index = fe_eval->get_current_cell_index();
    const typename Portable::MatrixFree<dim, double>::Data *data =
      fe_eval->get_matrix_free_data();

    const unsigned int position =
      data->local_q_point_id(cell_index, n_q_points, q_point);
    auto coeff = coef[position];

    auto value = fe_eval->get_value(q_point);

    fe_eval->submit_value(coeff * value, q_point);
    fe_eval->submit_gradient(fe_eval->get_gradient(q_point), q_point);
  }


  // @sect3{Class <code>LocalHelmholtzOperator</code>}

  // Finally, we need to define a class that implements the whole operator
  // evaluation that corresponds to a matrix-vector product in matrix-based
  // approaches.
  template <int dim, int fe_degree>
  class LocalHelmholtzOperator
  {
  public:
    // Again, the Portable::MatrixFree object doesn't know about the number
    // of degrees of freedom and the number of quadrature points so we need
    // to store these for index calculations in the call operator.
    static constexpr unsigned int n_local_dofs =
      Utilities::pow(fe_degree + 1, dim);
    static constexpr unsigned int n_q_points =
      Utilities::pow(fe_degree + 1, dim);

    LocalHelmholtzOperator(double *coefficient)
      : coef(coefficient)
    {}

    DEAL_II_HOST_DEVICE void
    operator()(const typename Portable::MatrixFree<dim, double>::Data *data,
               const Portable::DeviceVector<double>                   &src,
               Portable::DeviceVector<double> &dst) const;

  private:
    double *coef;
  };


  // This is the call operator that performs the Helmholtz operator evaluation
  // on a given cell similar to the MatrixFree framework on the CPU.
  // In particular, we need access to both values and gradients of the source
  // vector and we write value and gradient information to the destination
  // vector.
  template <int dim, int fe_degree>
  DEAL_II_HOST_DEVICE void LocalHelmholtzOperator<dim, fe_degree>::operator()(
    const typename Portable::MatrixFree<dim, double>::Data *data,
    const Portable::DeviceVector<double>                   &src,
    Portable::DeviceVector<double>                         &dst) const
  {
    Portable::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double> fe_eval(
      data);
    fe_eval.read_dof_values(src);
    fe_eval.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    fe_eval.apply_for_each_quad_point(
      HelmholtzOperatorQuad<dim, fe_degree>(coef));
    fe_eval.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    fe_eval.distribute_local_to_global(dst);
  }


  // @sect3{Class <code>HelmholtzOperator</code>}

  // The `HelmholtzOperator` class acts as wrapper for
  // `LocalHelmholtzOperator` defining an interface that can be used
  // with linear solvers like SolverCG. In particular, like every
  // class that implements the interface of a linear operator, it
  // needs to have a `vmult()` function that performs the action of
  // the linear operator on a source vector.
  template <int dim, int fe_degree>
  class HelmholtzOperator : public EnableObserverPointer
  {
  public:
    HelmholtzOperator(const DoFHandler<dim>           &dof_handler,
                      const AffineConstraints<double> &constraints);

    void
    vmult(LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &dst,
          const LinearAlgebra::distributed::Vector<double, MemorySpace::Default>
            &src) const;

    void initialize_dof_vector(
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &vec)
      const;

    void compute_diagonal();

    std::shared_ptr<DiagonalMatrix<
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>>
    get_matrix_diagonal_inverse() const;

    types::global_dof_index m() const;

    types::global_dof_index n() const;

    double el(const types::global_dof_index row,
              const types::global_dof_index col) const;

  private:
    Portable::MatrixFree<dim, double>                                mf_data;
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default> coef;
    std::shared_ptr<DiagonalMatrix<
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>>
      inverse_diagonal_entries;
  };



  // The following is the implementation of the constructor of this
  // class. In the first part, we initialize the `mf_data` member
  // variable that is going to provide us with the necessary
  // information when evaluating the operator.
  //
  // In the second half, we need to store the value of the coefficient
  // for each quadrature point in every active, locally owned cell.
  // We can ask the parallel triangulation for the number of active, locally
  // owned cells but only have a DoFHandler object at hand. Since
  // DoFHandler::get_triangulation() returns a Triangulation object, not a
  // parallel::TriangulationBase object, we have to downcast the return value.
  // This is safe to do here because we know that the triangulation is a
  // parallel::distributed::Triangulation object in fact.
  template <int dim, int fe_degree>
  HelmholtzOperator<dim, fe_degree>::HelmholtzOperator(
    const DoFHandler<dim>           &dof_handler,
    const AffineConstraints<double> &constraints)
  {
    const MappingQ<dim> mapping(fe_degree);
    typename Portable::MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_values | update_gradients |
                                           update_JxW_values |
                                           update_quadrature_points;
    const QGauss<1> quad(fe_degree + 1);
    mf_data.reinit(mapping, dof_handler, constraints, quad, additional_data);


    const unsigned int n_owned_cells =
      dynamic_cast<const parallel::TriangulationBase<dim> *>(
        &dof_handler.get_triangulation())
        ->n_locally_owned_active_cells();
    coef.reinit(Utilities::pow(fe_degree + 1, dim) * n_owned_cells);

    const VaryingCoefficientFunctor<dim, fe_degree> functor(coef.get_values());
    mf_data.evaluate_coefficients(functor);
  }


  // The key step then is to use all of the previous classes to loop over
  // all cells to perform the matrix-vector product. We implement this
  // in the next function.
  //
  // When applying the Helmholtz operator, we have to be careful to handle
  // boundary conditions correctly. Since the local operator doesn't know about
  // constraints, we have to copy the correct values from the source to the
  // destination vector afterwards.
  template <int dim, int fe_degree>
  void HelmholtzOperator<dim, fe_degree>::vmult(
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>       &dst,
    const LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &src)
    const
  {
    dst = 0.;
    LocalHelmholtzOperator<dim, fe_degree> helmholtz_operator(
      coef.get_values());
    mf_data.cell_loop(helmholtz_operator, src, dst);
    mf_data.copy_constrained_values(src, dst);
  }



  template <int dim, int fe_degree>
  void HelmholtzOperator<dim, fe_degree>::initialize_dof_vector(
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &vec) const
  {
    mf_data.initialize_dof_vector(vec);
  }



  template <int dim, int fe_degree>
  void HelmholtzOperator<dim, fe_degree>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(
      new DiagonalMatrix<
        LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>());
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>
      &inverse_diagonal = inverse_diagonal_entries->get_vector();
    initialize_dof_vector(inverse_diagonal);

    HelmholtzOperatorQuad<dim, fe_degree> helmholtz_operator_quad(
      coef.get_values());

    MatrixFreeTools::compute_diagonal<dim, fe_degree, fe_degree + 1, 1, double>(
      mf_data,
      inverse_diagonal,
      helmholtz_operator_quad,
      EvaluationFlags::values | EvaluationFlags::gradients,
      EvaluationFlags::values | EvaluationFlags::gradients);

    double *raw_diagonal = inverse_diagonal.get_values();

    Kokkos::parallel_for(
      inverse_diagonal.locally_owned_size(), KOKKOS_LAMBDA(int i) {
        Assert(raw_diagonal[i] > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        raw_diagonal[i] = 1. / raw_diagonal[i];
      });
  }



  template <int dim, int fe_degree>
  std::shared_ptr<DiagonalMatrix<
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>>
  HelmholtzOperator<dim, fe_degree>::get_matrix_diagonal_inverse() const
  {
    return inverse_diagonal_entries;
  }



  template <int dim, int fe_degree>
  types::global_dof_index HelmholtzOperator<dim, fe_degree>::m() const
  {
    return mf_data.get_vector_partitioner()->size();
  }



  template <int dim, int fe_degree>
  types::global_dof_index HelmholtzOperator<dim, fe_degree>::n() const
  {
    return mf_data.get_vector_partitioner()->size();
  }



  template <int dim, int fe_degree>
  double
  HelmholtzOperator<dim, fe_degree>::el(const types::global_dof_index row,
                                        const types::global_dof_index col) const
  {
    (void)col;
    Assert(row == col, ExcNotImplemented());
    Assert(inverse_diagonal_entries.get() != nullptr &&
             inverse_diagonal_entries->m() > 0,
           ExcNotInitialized());
    return 1.0 / (*inverse_diagonal_entries)(row, row);
  }



  // @sect3{Class <code>HelmholtzProblem</code>}

  // This is the main class of this program. It defines the usual
  // framework we use for tutorial programs. The only point worth
  // commenting on is the `solve()` function and the choice of vector
  // types.
  template <int dim, int fe_degree>
  class HelmholtzProblem
  {
  public:
    HelmholtzProblem();

    void run();

  private:
    void setup_system();

    void assemble_rhs();

    void solve();

    void output_results(const unsigned int cycle) const;

    MPI_Comm mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;

    const FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double>                          constraints;
    std::unique_ptr<HelmholtzOperator<dim, fe_degree>> system_matrix_dev;

    // Since all the operations in the `solve()` function are executed on the
    // graphics card, it is necessary for the vectors used to store their values
    // on the GPU as well. LinearAlgebra::distributed::Vector can be told which
    // memory space to use. It might
    // be worth noticing that the communication between different MPI processes
    // can be improved if the MPI implementation is GPU-aware and the configure
    // flag `DEAL_II_MPI_WITH_DEVICE_SUPPORT` is enabled. (The value of this
    // flag needs to be set at the time you call `cmake` when installing
    // deal.II.)
    //
    // In addition, we also keep a solution vector with CPU storage such that we
    // can view and display the solution as usual.
    LinearAlgebra::distributed::Vector<double, MemorySpace::Host>
      ghost_solution_host;
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>
      solution_dev;
    LinearAlgebra::distributed::Vector<double, MemorySpace::Default>
      system_rhs_dev;

    ConditionalOStream pcout;
  };


  // The implementation of all the remaining functions of this class apart from
  // `Helmholtzproblem::solve()` doesn't contain anything new and we won't
  // further comment much on the overall approach.
  template <int dim, int fe_degree>
  HelmholtzProblem<dim, fe_degree>::HelmholtzProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(mpi_communicator)
    , fe(fe_degree)
    , dof_handler(triangulation)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {}



  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::setup_system()
  {
    dof_handler.distribute_dofs(fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    system_rhs_dev.reinit(locally_owned_dofs, mpi_communicator);

    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close();

    system_matrix_dev.reset(
      new HelmholtzOperator<dim, fe_degree>(dof_handler, constraints));

    ghost_solution_host.reinit(locally_owned_dofs,
                               locally_relevant_dofs,
                               mpi_communicator);
    system_matrix_dev->initialize_dof_vector(solution_dev);
    system_rhs_dev.reinit(solution_dev);
  }



  // Unlike programs such as step-4 or step-6, we will not have to
  // assemble the whole linear system but only the right hand side
  // vector. This looks in essence like we did in step-4, for example,
  // but we have to pay attention to using the right constraints
  // object when copying local contributions into the global
  // vector. In particular, we need to make sure the entries that
  // correspond to boundary nodes are properly zeroed out. This is
  // necessary for CG to converge.  (Another solution would be to
  // modify the `vmult()` function above in such a way that we pretend
  // the source vector has zero entries by just not taking them into
  // account in matrix-vector products. But the approach used here is
  // simpler.)
  //
  // At the end of the function, we can't directly copy the values
  // from the host to the device but need to use an intermediate
  // object of type LinearAlgebra::ReadWriteVector to construct the
  // correct communication pattern necessary.
  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::assemble_rhs()
  {
    LinearAlgebra::distributed::Vector<double, MemorySpace::Host>
                      system_rhs_host(locally_owned_dofs,
                      locally_relevant_dofs,
                      mpi_communicator);
    const QGauss<dim> quadrature_formula(fe_degree + 1);

    FEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_rhs = 0;

          fe_values.reinit(cell);

          for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                cell_rhs(i) += (fe_values.shape_value(i, q_index) * 1.0 *
                                fe_values.JxW(q_index));
            }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 system_rhs_host);
        }
    system_rhs_host.compress(VectorOperation::add);

    LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
    rw_vector.import_elements(system_rhs_host, VectorOperation::insert);
    system_rhs_dev.import_elements(rw_vector, VectorOperation::insert);
  }



  // This solve() function finally contains the calls to the new classes
  // previously discussed. Here we don't use any preconditioner, i.e.,
  // precondition by the identity matrix, to focus just on the peculiarities of
  // the Portable::MatrixFree framework. Of course, in a real application
  // the choice of a suitable preconditioner is crucial but we have at least the
  // same restrictions as in step-37 since matrix entries are computed on the
  // fly and not stored.
  //
  // After solving the linear system in the first part of the function, we
  // copy the solution from the device to the host to be able to view its
  // values and display it in `output_results()`. This transfer works the same
  // as at the end of the previous function.
  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::solve()
  {
    system_matrix_dev->compute_diagonal();

    using PreconditionerType = PreconditionChebyshev<
      HelmholtzOperator<dim, fe_degree>,
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>;
    typename PreconditionerType::AdditionalData additional_data;
    additional_data.smoothing_range     = 15.;
    additional_data.degree              = 5;
    additional_data.eig_cg_n_iterations = 10;
    additional_data.preconditioner =
      system_matrix_dev->get_matrix_diagonal_inverse();

    PreconditionerType preconditioner;
    preconditioner.initialize(*system_matrix_dev, additional_data);

    SolverControl solver_control(system_rhs_dev.size(),
                                 1e-12 * system_rhs_dev.l2_norm());
    SolverCG<LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>
      cg(solver_control);
    cg.solve(*system_matrix_dev, solution_dev, system_rhs_dev, preconditioner);

    pcout << "  Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
    rw_vector.import_elements(solution_dev, VectorOperation::insert);
    ghost_solution_host.import_elements(rw_vector, VectorOperation::insert);

    constraints.distribute(ghost_solution_host);

    ghost_solution_host.update_ghost_values();
  }

  // The output results function is as usual since we have already copied the
  // values back from the GPU to the CPU.
  //
  // While we're already doing something with the function, we might
  // as well compute the $L_2$ norm of the solution. We do this by
  // calling VectorTools::integrate_difference(). That function is
  // meant to compute the error by evaluating the difference between
  // the numerical solution (given by a vector of values for the
  // degrees of freedom) and an object representing the exact
  // solution. But we can easily compute the $L_2$ norm of the
  // solution by passing in a zero function instead. That is, instead
  // of evaluating the error $\|u_h-u\|_{L_2(\Omega)}$, we are just
  // evaluating $\|u_h-0\|_{L_2(\Omega)}=\|u_h\|_{L_2(\Omega)}$
  // instead.
  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::output_results(
    const unsigned int cycle) const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(ghost_solution_host, "solution");
    data_out.build_patches();

    DataOutBase::VtkFlags flags;
    flags.compression_level = DataOutBase::CompressionLevel::best_speed;
    data_out.set_flags(flags);
    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, mpi_communicator, 2);

    Vector<float> cellwise_norm(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      ghost_solution_host,
                                      Functions::ZeroFunction<dim>(),
                                      cellwise_norm,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::L2_norm);
    const double global_norm =
      VectorTools::compute_global_error(triangulation,
                                        cellwise_norm,
                                        VectorTools::L2_norm);
    pcout << "  solution norm: " << global_norm << std::endl;
  }


  // There is nothing surprising in the `run()` function either. We simply
  // compute the solution on a series of (globally) refined meshes.
  template <int dim, int fe_degree>
  void HelmholtzProblem<dim, fe_degree>::run()
  {
    for (unsigned int cycle = 0; cycle < 7 - dim; ++cycle)
      {
        pcout << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          GridGenerator::hyper_cube(triangulation, 0., 1.);
        triangulation.refine_global(1);

        setup_system();

        pcout << "   Number of active cells:       "
              << triangulation.n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

        assemble_rhs();
        solve();
        output_results(cycle);
        pcout << std::endl;
      }
  }
} // namespace Step64


// @sect3{The <code>main()</code> function}
//
// Finally for the `main()` function.
// Kokkos needs to be initialized before being used, just as MPI does.
// Utilities::MPI::MPI_InitFinalize takes care of first initializing MPI and
// then Kokkos. This implies that Kokkos can take advantage of environment
// variables set by MPI such as OMPI_COMM_WORLD_LOCAL_RANK or
// OMPI_COMM_WORLD_LOCAL_SIZE to assign GPUs to MPI processes in a round-robin
// fashion. If such environment variables are not present, Kokkos uses the first
// visible GPU on every process. This might be suboptimal if that implies that
// multiple processes use the same GPU.
int main(int argc, char *argv[])
{
  try
    {
      using namespace Step64;

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

      HelmholtzProblem<3, 3> helmholtz_problem;
      helmholtz_problem.run();
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
