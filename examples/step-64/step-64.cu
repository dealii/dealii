/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the deal.II authors
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
#include <deal.II/base/cuda.h>

#include <deal.II/matrix_free/cuda_fe_evaluation.h>
#include <deal.II/matrix_free/cuda_matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <fstream>


// As usual, we enclose everything into a namespace of its own:
namespace Step64
{
  using namespace dealii;


  // @sect3{Class <code>VaryingCoefficientFunctor</code>}

  // Next, we define a class that implements the varying coefficients
  // we want to use in the Helmholtz operator. Later, we want to pass
  // an object of this type to a CUDAWrappers::MatrixFree
  // object that expects the class to have an `operator()` that fills the
  // values provided in the constructor for a given cell. This operator
  // needs to run on the device, so it needs to be marked as `__device__`
  // for the compiler.
  template <int dim, int fe_degree>
  class VaryingCoefficientFunctor
  {
  public:
    VaryingCoefficientFunctor(double *coefficient)
      : coef(coefficient)
    {}

    __device__ void operator()(
      const unsigned int                                          cell,
      const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data);

    // Since CUDAWrappers::MatrixFree::Data doesn't know about the size of its
    // arrays, we need to store the number of quadrature points and the numbers
    // of degrees of freedom in this class to do necessary index conversions.
    static const unsigned int n_dofs_1d = fe_degree + 1;
    static const unsigned int n_local_dofs =
      dealii::Utilities::pow(n_dofs_1d, dim);
    static const unsigned int n_q_points =
      dealii::Utilities::pow(n_dofs_1d, dim);

  private:
    double *coef;
  };



  // The following function implements this coefficient. Recall from
  // the introduction that we have defined it as $a(\mathbf
  // x)=\frac{10}{0.05 + 2\|\mathbf x\|^2}$
  template <int dim, int fe_degree>
  __device__ void VaryingCoefficientFunctor<dim, fe_degree>::operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data)
  {
    const unsigned int pos = CUDAWrappers::local_q_point_id<dim, double>(
      cell, gpu_data, n_dofs_1d, n_q_points);
    const Point<dim> q_point =
      CUDAWrappers::get_quadrature_point<dim, double>(cell,
                                                      gpu_data,
                                                      n_dofs_1d);

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
  // the device, so need to be marked as `__device__` for the
  // compiler.
  template <int dim, int fe_degree>
  class HelmholtzOperatorQuad
  {
  public:
    __device__ HelmholtzOperatorQuad(double coef)
      : coef(coef)
    {}

    __device__ void
    operator()(CUDAWrappers::FEEvaluation<dim, fe_degree> *fe_eval) const;

  private:
    double coef;
  };


  // The Helmholtz problem we want to solve here reads in weak form as follows:
  // @f{eqnarray*}
  //   (\nabla v, \nabla u)+ (v, a(\mathbf x) u) &=&(v,1) \quad \forall v.
  // @f}
  // If you have seen step-37, then it will be obvious that
  // the two terms on the left-hand side correspond to the two function calls
  // here:
  template <int dim, int fe_degree>
  __device__ void HelmholtzOperatorQuad<dim, fe_degree>::
                  operator()(CUDAWrappers::FEEvaluation<dim, fe_degree> *fe_eval) const
  {
    fe_eval->submit_value(coef * fe_eval->get_value());
    fe_eval->submit_gradient(fe_eval->get_gradient());
  }


  // @sect3{Class <code>LocalHelmholtzOperator</code>}

  // Finally, we need to define a class that implements the whole operator
  // evaluation that corresponds to a matrix-vector product in matrix-based
  // approaches.
  template <int dim, int fe_degree>
  class LocalHelmholtzOperator
  {
  public:
    LocalHelmholtzOperator(double *coefficient)
      : coef(coefficient)
    {}

    __device__ void operator()(
      const unsigned int                                          cell,
      const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data,
      CUDAWrappers::SharedData<dim, double> *                     shared_data,
      const double *                                              src,
      double *                                                    dst) const;

    // Again, the CUDAWrappers::MatrixFree object doesn't know about the number
    // of degrees of freedom and the number of quadrature points so we need
    // to store these for index calculations in the call operator.
    static const unsigned int n_dofs_1d    = fe_degree + 1;
    static const unsigned int n_local_dofs = Utilities::pow(fe_degree + 1, dim);
    static const unsigned int n_q_points   = Utilities::pow(fe_degree + 1, dim);

  private:
    double *coef;
  };


  // This is the call operator that performs the Helmholtz operator evaluation
  // on a given cell similar to the MatrixFree framework on the CPU.
  // In particular, we need access to both values and gradients of the source
  // vector and we write value and gradient information to the destination
  // vector.
  template <int dim, int fe_degree>
  __device__ void LocalHelmholtzOperator<dim, fe_degree>::operator()(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, double>::Data *gpu_data,
    CUDAWrappers::SharedData<dim, double> *                     shared_data,
    const double *                                              src,
    double *                                                    dst) const
  {
    const unsigned int pos = CUDAWrappers::local_q_point_id<dim, double>(
      cell, gpu_data, n_dofs_1d, n_q_points);

    CUDAWrappers::FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double>
      fe_eval(cell, gpu_data, shared_data);
    fe_eval.read_dof_values(src);
    fe_eval.evaluate(true, true);
    fe_eval.apply_for_each_quad_point(
      HelmholtzOperatorQuad<dim, fe_degree>(coef[pos]));
    fe_eval.integrate(true, true);
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
  class HelmholtzOperator
  {
  public:
    HelmholtzOperator(const DoFHandler<dim> &          dof_handler,
                      const AffineConstraints<double> &constraints);

    void
    vmult(LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &dst,
          const LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>
            &src) const;

    void initialize_dof_vector(
      LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &vec) const;

  private:
    CUDAWrappers::MatrixFree<dim, double>       mf_data;
    LinearAlgebra::CUDAWrappers::Vector<double> coef;
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
  // parallel::Triangulation object, we have to downcast the return value. This
  // is safe to do here because we know that the triangulation is a
  // parallel:distributed::Triangulation object in fact.
  template <int dim, int fe_degree>
  HelmholtzOperator<dim, fe_degree>::HelmholtzOperator(
    const DoFHandler<dim> &          dof_handler,
    const AffineConstraints<double> &constraints)
  {
    MappingQGeneric<dim> mapping(fe_degree);
    typename CUDAWrappers::MatrixFree<dim, double>::AdditionalData
      additional_data;
    additional_data.mapping_update_flags = update_values | update_gradients |
                                           update_JxW_values |
                                           update_quadrature_points;
    const QGauss<1> quad(fe_degree + 1);
    mf_data.reinit(mapping, dof_handler, constraints, quad, additional_data);


    const unsigned int n_owned_cells =
      dynamic_cast<const parallel::Triangulation<dim> *>(
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
    LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &      dst,
    const LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &src)
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
    LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> &vec) const
  {
    mf_data.initialize_dof_vector(vec);
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

    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double>                          constraints;
    std::unique_ptr<HelmholtzOperator<dim, fe_degree>> system_matrix_dev;

    // Since all the operations in the `solve()` function are executed on the
    // graphics card, it is necessary for the vectors used to store their values
    // on the GPU as well. LinearAlgebra::distributed::Vector can be told which
    // memory space to use. There is also LinearAlgebra::CUDAWrappers::Vector
    // that always uses GPU memory storage but doesn't work with MPI. It might
    // be worth noticing that the communication between different MPI processes
    // can be improved if the MPI implementation is CUDA-aware and the configure
    // flag `DEAL_II_MPI_WITH_CUDA_SUPPORT` is enabled. (The value of this
    // flag needs to be set at the time you call `cmake` when installing
    // deal.II.)
    //
    // In addition, we also keep a solution vector with CPU storage such that we
    // can view and display the solution as usual.
    LinearAlgebra::distributed::Vector<double, MemorySpace::Host>
                                                                  ghost_solution_host;
    LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> solution_dev;
    LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>
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
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    system_rhs_dev.reinit(locally_owned_dofs, mpi_communicator);

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
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

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
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
    rw_vector.import(system_rhs_host, VectorOperation::insert);
    system_rhs_dev.import(rw_vector, VectorOperation::insert);
  }



  // This solve() function finally contains the calls to the new classes
  // previously discussed. Here we don't use any preconditioner, i.e.,
  // precondition by the identity matrix, to focus just on the peculiarities of
  // the CUDAWrappers::MatrixFree framework. Of course, in a real application
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
    PreconditionIdentity preconditioner;

    SolverControl solver_control(system_rhs_dev.size(),
                                 1e-12 * system_rhs_dev.l2_norm());
    SolverCG<LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA>> cg(
      solver_control);
    cg.solve(*system_matrix_dev, solution_dev, system_rhs_dev, preconditioner);

    pcout << "  Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);
    rw_vector.import(solution_dev, VectorOperation::insert);
    ghost_solution_host.import(rw_vector, VectorOperation::insert);

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
    flags.compression_level = DataOutBase::VtkFlags::best_speed;
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

// Finally for the `main()` function.  By default, all the MPI ranks
// will try to access the device with number 0, which we assume to be
// the GPU device associated with the CPU on which a particular MPI
// rank runs. This works, but if we are running with MPI support it
// may be that multiple MPI processes are running on the same machine
// (for example, one per CPU core) and then they would all want to
// access the same GPU on that machine. If there is only one GPU in
// the machine, there is nothing we can do about it: All MPI ranks on
// that machine need to share it. But if there are more than one GPU,
// then it is better to address different graphic cards for different
// processes. The choice below is based on the MPI process id by
// assigning GPUs round robin to GPU ranks. (To work correctly, this
// scheme assumes that the MPI ranks on one machine are
// consecutive. If that were not the case, then the rank-GPU
// association may just not be optimal.) To make this work, MPI needs
// to be initialized before using this function.
int main(int argc, char *argv[])
{
  try
    {
      using namespace Step64;

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

      int         n_devices       = 0;
      cudaError_t cuda_error_code = cudaGetDeviceCount(&n_devices);
      AssertCuda(cuda_error_code);
      const unsigned int my_mpi_id =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      const int device_id = my_mpi_id % n_devices;
      cuda_error_code     = cudaSetDevice(device_id);
      AssertCuda(cuda_error_code);

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
