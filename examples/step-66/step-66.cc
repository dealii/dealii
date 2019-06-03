/* ---------------------------------------------------------------------
 * 
 * Copyright (C) 2003 - 2018 by the deal.II authors
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
 * Authors: Fabian Castelli, Karlsruhe Institute of Technology (KIT)
 *          Timo Heister, University of Utah
 */



// First we include the typical headers of the deal.II library needed for this
// tutroial:
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_c1.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

// In particular, we need to include the headers for the matrix-free framework:
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

// And since we want to use a geometric multigrid preconditioner, we need also
// the multilevel headers:
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>


// Finally some common C++ headers for in and output:
#include <iostream>
#include <fstream>



namespace Step66
{
  using namespace dealii;
  
  
  
  // As in the previous tutorials on the matrix-free framework we define the
  // degree of the finite element space and the space dimension as constants
  // variables. For this example we solve a tow dimensional problem and use the
  // fourth order Lagrangian finite element space.
  const unsigned int degree_finite_element = 4;
  const unsigned int dimension             = 2;
  
  
  
  // sect3{Matrix-free operators}
  
  // In the beginning we define the two matrix-free operator for the Jacobian
  // and the residual. As guideline we follow the tutorials step-37 and step-48,
  // where the precise interface of the MatrixFreeOperators::Base class was
  // extensively documented.
  
  
  // @setc4{JacobianOperator}
  
  // Since we want to use the Jacobian as system matrix and pass it to the
  // linear solver as well as to the multilevel preconditioner classes we derive
  // the JacobianOperator class from the MatrixFreeOperators::Base class, such
  // that we have already the right interface. The two functions we need to
  // reimplement from the base class are the apply_add and the compute_diagonal
  // function. To allow preconditioning with float precision we define the
  // number type as template argument.
  //
  // As mentioned already in the introduction, for the Jacobian $F'$ of the
  // $n+1$-th Newton step we need to evaluate it at the last Newton step
  // $u_h^n$. To get the information of the last Newton step $u_h^n$ we do
  // pretty much the same as we in step-37, where we stored the values of a
  // coefficient function in a table once before we use the matrix-free
  // operator. Instead of a function evaluate_coefficient, we here implement a
  // function evaluate_newton_step.
  //
  // As private member functions of the JacobianOperator we implement the
  // local_apply and the local_compute_diagonal function, which we call in the
  // cell_loop for computing the matrix-vector product or the diagonal.
  template <int dim, int fe_degree, typename number>
  class JacobianOperator : public MatrixFreeOperators::Base<dim,LinearAlgebra::distributed::Vector<number> >
  {
  public:
    using value_type = number;
    
    JacobianOperator();
    
    virtual
    void
    clear() override;
    
    void
    evaluate_newton_step(const LinearAlgebra::distributed::Vector<number> &src);
    
    virtual
    void
    compute_diagonal() override;
    
  private:
    virtual
    void
    apply_add(LinearAlgebra::distributed::Vector<number> &dst, const LinearAlgebra::distributed::Vector<number> &src) const override;
    
    void
    local_apply(const MatrixFree<dim,number> &data, LinearAlgebra::distributed::Vector<number> &dst, const LinearAlgebra::distributed::Vector<number> &src, const std::pair<unsigned int,unsigned int> &cell_range) const;
    
    void
    local_compute_diagonal(const MatrixFree<dim,number> &data, LinearAlgebra::distributed::Vector<number> &dst, const unsigned int &dummy, const std::pair<unsigned int,unsigned int> &cell_range) const;
    
    Table<2, VectorizedArray<number> > nonlinear_values;
  };
  
  
  
  // The constructor of the JacobainOperator just calls the constructor of the
  // base class MatrixFreeOperators::Base, which is itself derived from the
  // Subscriptor class.
  template <int dim, int fe_degree, typename number>
  JacobianOperator<dim,fe_degree,number>::JacobianOperator()
  :
  MatrixFreeOperators::Base<dim,LinearAlgebra::distributed::Vector<number> >()
  {}
  
  
  
  // The clear function resets out the table holding the values for the
  // nonlinearity and call the clear function of the base class.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::clear()
  {
    nonlinear_values.reinit(0, 0);
    MatrixFreeOperators::Base<dim,LinearAlgebra::distributed::Vector<number> >::clear();
  }
  
  
  
  // The following evlauate_newton_step function is based on the
  // evlauate_coefficient function from step-37, however, it does not evaluate a
  // function object, it evaluates a vector representing a finite element
  // function, namely the last Newton step needed for the Jacobian. Therefore we
  // set up a FEEvaluateion object and evaluate the finite element function in
  // the quadrature points with the gather_evaluate function. Remember if we
  // only loop over cells this function is just a wrapper around functions
  // read_dof_values and evaluate. Since we store the evaluated values of the
  // finite element function in a table we do not have to call integrate in
  // combination with distribute_local_to_global or integrate scatter.
  //
  // This will work well and in the local_apply function we can use the values
  // stored in the table to apply the matrix-vector product, however, we can
  // also optimize the implementation of the Jacobian at this stage. We can
  // directly evaluate the nonlinear function std::exp(u_h^n[q]) and store these
  // values in the table. This saves in each call of the vmult function all
  // evaluations of the nonlinearity.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::evaluate_newton_step(const LinearAlgebra::distributed::Vector<number> &src)
  {
    const unsigned int n_cells = this->data->n_macro_cells();
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi(*this->data);
    
    nonlinear_values.reinit(n_cells, phi.n_q_points);
    for(unsigned int cell=0; cell<n_cells; ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(src, true, false);
      
      for(unsigned int q=0; q<phi.n_q_points; ++q)
        nonlinear_values(cell, q) = std::exp(phi.get_value(q));
    }
  }
  
  
  
  // Now in the local_apply function, which actually implements the local action
  // of the system matrix on a cell, we can use the information of the last
  // Newton step stored in the table nonlinear_values. The rest of this function
  // as basically the same in step-37. We setup the FEEvaluation object, gather
  // and evaluate the values and gradients of the input vector, submit the
  // values and gradients according to the form of the Jacobian and finally call
  // integrate_scatter to distribute the local contributions into the global
  // vector dst.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::local_apply(const MatrixFree<dim,number> &data, LinearAlgebra::distributed::Vector<number> &dst, const LinearAlgebra::distributed::Vector<number> &src, const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi(data);
    
    for(unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      AssertDimension(nonlinear_values.size(0), data.n_macro_cells());
      AssertDimension(nonlinear_values.size(1), phi.n_q_points);
      
      phi.reinit(cell);
      phi.gather_evaluate(src, true, true);
      
      for(unsigned int q=0; q<phi.n_q_points; ++q)
      {
        phi.submit_value(-nonlinear_values(cell,q)*phi.get_value(q), q);
        phi.submit_gradient(phi.get_gradient(q), q);
      }
      
      phi.integrate_scatter(true, true, dst);
    }
  }
  
  
  
  // Next we use cell_loop of the MatrixFree class to perform the actual loop
  // over all cells computing the cell contribution to the matrix-vector
  // product.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::apply_add(LinearAlgebra::distributed::Vector<number> &dst, const LinearAlgebra::distributed::Vector<number> &src) const
  {
    this->data->cell_loop(&JacobianOperator::local_apply, this, dst, src, true);
  }
  
  
  
  // The missing two functions of the JacobianOperator calculate the diagonal
  // entries of the Jacobian. The only difference compared to step-37 is the
  // calculation of the cell contribution in the local_compute_diagonal
  // function. However, therefore we only have to extend and change the
  // arguments for the submit functions in the loop over all quadrature points
  // and this can be done according to the local_apply function. So no further
  // comments to these two functions schould be necessary.
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(new DiagonalMatrix<LinearAlgebra::distributed::Vector<number> >());
    LinearAlgebra::distributed::Vector<number> &inverse_diagonal = this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);
    
    unsigned int dummy = 0;
    
    this->data->cell_loop(&JacobianOperator::local_compute_diagonal, this, inverse_diagonal, dummy);
    
    this->set_constrained_entries_to_one(inverse_diagonal);
    
    for(unsigned int i=0; i<inverse_diagonal.local_size(); ++i)
    {
      Assert(inverse_diagonal.local_element(i) > 0., ExcMessage("No diagonal entry in a positive definite operator should be zero"));
      inverse_diagonal.local_element(i) = 1./inverse_diagonal.local_element(i);
    }
  }
  
  
  
  template <int dim, int fe_degree, typename number>
  void
  JacobianOperator<dim,fe_degree,number>::local_compute_diagonal(const MatrixFree<dim,number> &data, LinearAlgebra::distributed::Vector<number> &dst, const unsigned int &, const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi(data);
    
    AlignedVector<VectorizedArray<number> > diagonal(phi.dofs_per_cell);
    
    for(unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      AssertDimension(nonlinear_values.size(0), data.n_macro_cells());
      AssertDimension(nonlinear_values.size(1), phi.n_q_points);
      
      phi.reinit(cell);
      for(unsigned int i=0; i<phi.dofs_per_cell; ++i)
      {
        for(unsigned int j=0; j<phi.dofs_per_cell; ++j)
          phi.submit_dof_value(VectorizedArray<number>(), j);
        phi.submit_dof_value(make_vectorized_array<number>(1.), i);
        
        phi.evaluate(true, true);
        for(unsigned int q=0; q<phi.n_q_points; ++q)
        {
          phi.submit_value(-nonlinear_values(cell, q)*phi.get_value(q), q);
          phi.submit_gradient(phi.get_gradient(q), q);
        }
        phi.integrate(true, true);
        diagonal[i] = phi.get_dof_value(i);
      }
      for(unsigned int i=0; i<phi.dofs_per_cell; ++i)
        phi.submit_dof_value(diagonal[i], i);
      phi.distribute_local_to_global(dst);
    }
  }
  
  
  
  // @sect4{ResidualOperator}
  
  // The next class implements the ResidualOperator, which we use to evaluate
  // the residual and to assemble the right hand side. Since we will not pass
  // objects of this class to any other objects for example linear solvers, etc.
  // we have not need to fulfill any interface requirements. The operator's main
  // functionality is implemented in the apply function. We can actually think
  // of ResidualOperator.apply(dst, src) is mathematically expressed as
  // F(src) = dst, where $F\colon\mathbb{R}^N\to\mathbb{R}^N$ is the discrete
  // residual introduced in the introduction. Further the ResidualOperator has
  // an initialization function, which stores a shared pointer to the MatrixFree
  // object handling the loop over all cells, and a local_apply function
  // implementing the calculation of the cell contribution.
  
  // ??? As functor with operator() function.
  
  template <int dim, int fe_degree>
  class ResidualOperator
  {
  public:
    ResidualOperator() = default;
    
    void
    initialize(std::shared_ptr<const MatrixFree<dim,double> > data_in);
    
    void
    apply(LinearAlgebra::distributed::Vector<double> &dst, const LinearAlgebra::distributed::Vector<double> &src) const;
          
  private:
    void
    local_apply(const MatrixFree<dim,double> &data, LinearAlgebra::distributed::Vector<double> &dst, const LinearAlgebra::distributed::Vector<double> &src, const std::pair<unsigned int,unsigned int> &cell_range) const;
    
    std::shared_ptr<const MatrixFree<dim,double> > data;
  };
  
  
  
  // The initialize function checks if given shared pointer of the MatrixFree
  // object is not empty and pass it to the class varaiable.
  template <int dim, int fe_degree>
  void
  ResidualOperator<dim,fe_degree>::initialize(std::shared_ptr<const MatrixFree<dim,double> > data_in)
  {
    Assert(data_in, ExcNotInitialized());
    
    data = data_in;
  }
  
  
  
  // The main function evaluating the residual is just a wrapper around the
  // cell_loop of the MatrixFree object.
  template <int dim, int fe_degree>
  void
  ResidualOperator<dim, fe_degree>::apply(LinearAlgebra::distributed::Vector<double> &dst, const LinearAlgebra::distributed::Vector<double> &src) const
  {
    Assert(data, ExcNotInitialized());
    
    data->cell_loop(&ResidualOperator<dim,fe_degree>::local_apply, this, dst, src, true);
  }
  
  
  
  // The heart of the ResidualOperator is the local_apply function. The
  // implementation is similar to the local_apply function of the above
  // JacobianOperator. We setup a FEEvaluation object, gather and evaluate the
  // values and gradients, submit the new values and gradients and finally
  // integrate and distribute the local contributions to the global vector.
  // Different to the Jacobian we need no additional table storing the values of
  // the old Newton step, instead we can evaluate the nonlinearity on the fly,
  // since we have to evaluate the residual in the input vector src,
  // representing the last Newton step.
  template <int dim, int fe_degree>
  void
  ResidualOperator<dim, fe_degree>::local_apply(const MatrixFree<dim> &data, LinearAlgebra::distributed::Vector<double> &dst, const LinearAlgebra::distributed::Vector<double> &src, const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree> phi(data);
    
    for(unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      phi.reinit(cell);
      phi.gather_evaluate(src, true, true);
      
      for(unsigned int q=0; q<phi.n_q_points; ++q)
      {
        phi.submit_value(-std::exp(phi.get_value(q)), q);
        phi.submit_gradient(phi.get_gradient(q), q);
      }
      
      phi.integrate_scatter(true, true, dst);
    }
  }
  
  
  
  // @sect3{GelfandProblem class}
  
  // After implementing the matrix-free operators we can now define the solver
  // class for the <i>Gelfand problem</i>. This class is based on the common
  // structure of all previous tutorial programs, in particluar it is based on
  // step-15, solving also a nonlinear problem. Since we using the matrix-free
  // framework, we need no assemble_system function any more, instead the
  // information of the matrix is rebuild in every call of the vmult function.
  // However, for the nonlinear solver we need additional functions for the
  // computation of the residual and the Newton method.
  template <int dim>
  class GelfandProblem
  {
  public:
    GelfandProblem();
    
    void
    run();
    
  private:
    void
    make_grid();
    
    void
    setup_system();
    
    void
    assemble_rhs();
    
    double
    compute_residual(const double alpha);
    
    void
    compute_update();
    
    void
    solve();
    
    double
    compute_solution_norm() const;
    
    void
    output_results(const unsigned int cycle) const;
    
    
    // For the triangulation we use directly the parallel::distributed version
    // and define an object for the C^1 mapping.
    parallel::distributed::Triangulation<dim> triangulation;
    const MappingC1<dim>                      mapping;
    
    FE_Q<dim>                                 fe;
    DoFHandler<dim>                           dof_handler;
    
    
    // For the system matrix we have defined the class JacobianOperator, which
    // is a matrix-free operator. For the assembly of the right hand side and
    // the resiudal we use the matrix-free ResidualOperator also implemented
    // above.
    ConstraintMatrix                                   constraints;
    JacobianOperator<dim,degree_finite_element,double> system_matrix;
    ResidualOperator<dim,degree_finite_element>        residual_operator;
    
    
    // The multilevel object which are also based on the matrix-free operator
    // for the Jacobian. Since we need to evaluate the Jacobian with the last
    // Newton step we also need to evaluate the level operator with last Newton
    // step for the preconditioner. For this reason we need beside the
    // mg_matrices also a MGLevelObject to store the interpolated solution
    // vector on each level. As in step-37 we use for the preconditioner float
    // precision.
    MGConstrainedDoFs                                         mg_constrained_dofs;
    typedef JacobianOperator<dim,degree_finite_element,float> LevelMatrixType;
    MGLevelObject<LevelMatrixType>                            mg_matrices;
    MGLevelObject<LinearAlgebra::distributed::Vector<float> > mg_solution;
    
    
    // Of course we need also vector holding the solution, the Newton update and
    // the system right hand side. In that way we can store the last Newton step
    // always in the solution vector and just add the update to get the next
    // Newton step.
    LinearAlgebra::distributed::Vector<double> solution;
    LinearAlgebra::distributed::Vector<double> newton_update;
    LinearAlgebra::distributed::Vector<double> system_rhs;
    
    
    // Finally we have a variable for the number of iterations of the linear
    // solver.
    unsigned int linear_iterations;
    
    
    // For the output we use in programs running in parallel with MPI the
    // ConditionalOStream class to avoid multiple output of the same data by
    // different MPI ranks.
    ConditionalOStream	pcout;
    
    
    // Finally for the time measurement we use a TimerOutput object, which
    // prints in a nice table the elapsed CPU and wall times for each function
    // after the program is finished.
    TimerOutput computing_timer;
  };
  
  
  
  // The constructor of the GelfandProblem initializes the class variables. In
  // particluar, we setup the multilevel support for the
  // parallel::distributed::Triangulation<dim>, initilaize the
  // ConditionalOStream and tell the TimerOutput taht we want to see the CPU and
  // wall time in the end of the program.
  template <int dim>
  GelfandProblem<dim>::GelfandProblem()
  :
  triangulation(MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices, parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy),
  fe(degree_finite_element),
  dof_handler(triangulation),
  pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
  computing_timer(MPI_COMM_WORLD, pcout, TimerOutput::summary, TimerOutput::cpu_and_wall_times)
  {}
  
  
  
  // @sect4{GelfandProblem::make_grid}
  
  // As computational domain we use the two-dimensional unit circle. We follow
  // the instructions for the TransfiniteInterpolationManifold class and assign
  // also a SphericalManifold for the boundaty. In the end we use a six times
  // globally refined triangulation.
  
  // ??? Issue #6663
  template <int dim>
  void
  GelfandProblem<dim>::make_grid()
  {
    TimerOutput::Scope t(computing_timer, "make grid");
    
    SphericalManifold<dim> boundary_manifold;
    TransfiniteInterpolationManifold<dim> inner_manifold;
    
    GridGenerator::hyper_ball(triangulation);
    
    triangulation.set_all_manifold_ids(1);
    triangulation.set_all_manifold_ids_on_boundary(0);
    
    triangulation.set_manifold(0, boundary_manifold);
    
    inner_manifold.initialize(triangulation);
    triangulation.set_manifold(1, inner_manifold);
    
    triangulation.refine_global(3-dim);
  }
  
  
  
  // @setc4{GelfandProblem::setup_system}
  
  // The setup_system function is quasi identical to the one in step-37. The
  // only differences are obviously the time measurement with only one
  // TimerOutput::Scope instead of measuring each part individually and more
  // important the initialization of the ResidualOperator with the MatrixFree
  // object as well as the initialization of the MGLevelObject for the
  // interpolated solution vector of the last Newton.
  //
  // Note, how we can use the same MatrixFree object twice, for the
  // JacobianOperator and the ResidualOperator.
  template <int dim>
  void
  GelfandProblem<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup system");
    
    system_matrix.clear();
    mg_matrices.clear_elements();
    
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();
    
    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
    
    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
    constraints.close();
    
    {
      typename MatrixFree<dim,double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_color;
      additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values | update_quadrature_points);
      std::shared_ptr<MatrixFree<dim,double> > system_mf_storage(new MatrixFree<dim,double>());
      system_mf_storage->reinit(dof_handler, constraints, QGauss<1>(fe.degree+1), additional_data);
      
      system_matrix.initialize(system_mf_storage);
      residual_operator.initialize(system_mf_storage);
    }
    
    system_matrix.initialize_dof_vector(solution);
    system_matrix.initialize_dof_vector(newton_update);
    system_matrix.initialize_dof_vector(system_rhs);
    
    
    const unsigned int nlevels = triangulation.n_global_levels();
    mg_matrices.resize(0, nlevels-1);
    mg_solution.resize(0, nlevels-1);
    
    std::set<types::boundary_id> dirichlet_boundary;
    dirichlet_boundary.insert(0);
    mg_constrained_dofs.initialize(dof_handler);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, dirichlet_boundary);
    
    for(unsigned int level=0; level<nlevels; ++level)
    {
      IndexSet relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(dof_handler, level, relevant_dofs);
      
      ConstraintMatrix level_constraints;
      level_constraints.reinit(relevant_dofs);
      level_constraints.add_lines(mg_constrained_dofs.get_boundary_indices(level));
      level_constraints.close();
      
      typename MatrixFree<dim,float>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme = MatrixFree<dim,float>::AdditionalData::partition_color;
      additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values | update_quadrature_points);
      additional_data.level_mg_handler = level;
      std::shared_ptr<MatrixFree<dim,float> > mg_mf_storage_level(new MatrixFree<dim,float>());
      mg_mf_storage_level->reinit(dof_handler, level_constraints, QGauss<1>(fe.degree+1), additional_data);
      
      mg_matrices[level].initialize(mg_mf_storage_level, mg_constrained_dofs, level);
      mg_matrices[level].initialize_dof_vector(mg_solution[level]);
    }
  }
  
  
  
  // @sect4{GelfandProblem::assemble_rhs}
  
  // Using the implementation of the ResidualOperator the assembly of the right
  // hand side becomes now a very easy task. We just call the main function of
  // the residual operator and multiply the result with minus one.
  //
  // Experiences show that using the cell_loop of the MatrixFree class together
  // with the FEEvaluation class is much faster than a classical implementation.
  template <int dim>
  void
  GelfandProblem<dim>::assemble_rhs()
  {
    TimerOutput::Scope t(computing_timer, "assemble right hand side");
    
    residual_operator.apply(system_rhs, solution);
    
    system_rhs *= -1.0;
  }
  
  
  
  // @sect4{GelfandProblem::compute_residual}
  
  // With the help of the ResidualOperator the next function computes the
  // residual of the nonlinear problem. The input argument alpha becomes
  // important if we would use an adaptive version of the Newton method. Then
  // for example we would compute the residual for different step lengths and
  // compare the residuals. However for our problem the full Newton step with
  // $\alpha=1$ is the best we can do. An adaptive verison of Newton's method
  // becomes interesting if we have no good initial value. Note, the theory
  // states that Newton's method converges with quadratic order, but only if we
  // have an appropriate initial value. For wrong initial values the Newton
  // method diverges even with quadratic order. A common way is then to use a
  // damped version $\alpha<1$ unitl the Newton step is good enough and the full
  // Newton step can be performed. This was also discussed in step-15.
  template <int dim>
  double
  GelfandProblem<dim>::compute_residual(const double alpha)
  {
    TimerOutput::Scope t(computing_timer, "compute residual");
    
    LinearAlgebra::distributed::Vector<double> residual;
    LinearAlgebra::distributed::Vector<double> evaluation_point;
    
    system_matrix.initialize_dof_vector(residual);
    system_matrix.initialize_dof_vector(evaluation_point);
    
    evaluation_point = solution;
    if(alpha > 1e-12)
      evaluation_point.add(alpha, newton_update);
    
    residual_operator.apply(residual, evaluation_point);
    
    return residual.l2_norm();
  }
  
  
  
  // @sect4{GelfandProblem::compute_update}
  
  // In order to compute the Newton updates in each step Newton step we solve
  // the linear system with the cg algorithm together with a geometric multigrid
  // preconditioner. For this we first setup the PreconditionMG object with a
  // Chebyshev smoother like we did in step-37.
  template <int dim>
  void
  GelfandProblem<dim>::compute_update()
  {
    TimerOutput::Scope t(computing_timer, "compute update");
    
    // We remember that the Jacobian depends on the last Newton step stored in
    // the solution vector. So we update the ghost values of the Newton step and
    // pass it to the JacobianOperator to store the information.
    solution.update_ghost_values();
    
    system_matrix.evaluate_newton_step(solution);
    
    
    // Next we have to pass the last Newton step also to the multilevel
    // operators. However, therefore we need to interpolate the Newton step to
    // all levels of the triangulation. This is done with the interpolate_to_mg
    // function of the MGTransferMatrixFree class.
    MGTransferMatrixFree<dim,float> mg_transfer(mg_constrained_dofs);
    mg_transfer.build(dof_handler);
    
    mg_transfer.interpolate_to_mg(dof_handler, mg_solution, solution);
    
    
    // Now we can setup the preconditioner. We define the smoother and pass the
    // interpolated vectors of the Newton step to the multilevel operators. 
    typedef PreconditionChebyshev<LevelMatrixType,LinearAlgebra::distributed::Vector<float> > SmootherType;
    mg::SmootherRelaxation<SmootherType, LinearAlgebra::distributed::Vector<float> > mg_smoother;
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
    smoother_data.resize(0, triangulation.n_global_levels()-1);
    for(unsigned int level = 0; level<triangulation.n_global_levels(); ++level)
    {
      if(level > 0)
      {
        smoother_data[level].smoothing_range = 15.;
        smoother_data[level].degree = 4;
        smoother_data[level].eig_cg_n_iterations = 10;
      }
      else
      {
        smoother_data[0].smoothing_range = 1e-3;
        smoother_data[0].degree = numbers::invalid_unsigned_int;
        smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
      }
      
      mg_matrices[level].evaluate_newton_step(mg_solution[level]);
      mg_matrices[level].compute_diagonal();
      
      smoother_data[level].preconditioner = mg_matrices[level].get_matrix_diagonal_inverse();
    }
    mg_smoother.initialize(mg_matrices, smoother_data);
    
    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float> > mg_coarse;
    mg_coarse.initialize(mg_smoother);
    
    mg::Matrix<LinearAlgebra::distributed::Vector<float> > mg_matrix(mg_matrices);
    
    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType> > mg_interface_matrices;
    mg_interface_matrices.resize(0, triangulation.n_global_levels()-1);
    for(unsigned int level=0; level<triangulation.n_global_levels(); ++level)
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    mg::Matrix<LinearAlgebra::distributed::Vector<float> > mg_interface(mg_interface_matrices);
    
    Multigrid<LinearAlgebra::distributed::Vector<float> > mg(mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
    mg.set_edge_matrices(mg_interface, mg_interface);
    
    PreconditionMG<dim, LinearAlgebra::distributed::Vector<float>, MGTransferMatrixFree<dim,float> > preconditioner(dof_handler, mg, mg_transfer);
    
    
    // Having now a geometric multigrid preconditioner for the Jacobian we solve
    // the linear system with the cg algorithm.
    SolverControl solver_control(100, 1.e-12);
    SolverCG<LinearAlgebra::distributed::Vector<double> > cg(solver_control);
    
    constraints.set_zero(solution);
    
    cg.solve(system_matrix, newton_update, system_rhs, preconditioner);
    
    constraints.distribute(newton_update);
    
    linear_iterations = solver_control.last_step();
  }
  
  
  
  // @sect4{GelfandProblem::solve}
  
  // Now we implement the actual Newton solver for the nonlinear problem.  
  template<int dim>
  void
  GelfandProblem<dim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "solve");
    
    
    // We define a maximal number of Newton step and tolerances for the
    // convergence criterion. Usually, for good starting values, the Newton
    // method converges in three to six steps, so maximal ten steps should be
    // totally sufficient. As tolerances we use for the residual
    // $\|F(u^n_h)\|<\text{TOL}_f = 10^{-12}$ and
    // $\|s_h^n\| < \text{TOL}_x = 10^{-10}$.
    // This seems a bit over the top, but we will see that, for our example, we
    // will achieve these tolerances after a few steps.
    const unsigned int itmax = 10;
    const double       TOLf  = 1e-12;
    const double       TOLx  = 1e-10;
    
    
    Timer solver_timer;
    solver_timer.start();
    
    
    // In a loop over the Newton steps we first assemble the right hand side for
    // the linear problem and then compute the update.
    for(unsigned int newton_step=1; newton_step<=itmax; ++newton_step)
    {
      assemble_rhs();
      
      compute_update();
      
      
      // Then we compute the errors, namely the norm of the Newton update and
      // the residual.
      const double ERRx = newton_update.l2_norm();
      const double ERRf = compute_residual(1.0);
      
      
      // Compute the next Newton step by adding the update. A short output will
      // inform us on the current Newton step.
      solution.add(1.0, newton_update);
      
      pcout << "   Nstep " << newton_step << ", errf = " << ERRf << ", errx = " << ERRx << ", it = " << linear_iterations << std::endl;
      
      
      
      // After each Newton step we check the convergence criterions. If at least
      // one of those is fulfilled we end break the loop with success. Else we
      // check if we have computed already to much Newton steps. If this is not
      // the case we start from top and compute the next newotn step.
      if(ERRf < TOLf || ERRx < TOLx)
      {
        solver_timer.stop();
        
        pcout << "Convergence step " << newton_step << " value " << ERRf << " (used wall time: " << solver_timer.wall_time() << " s)" << std::endl;
        
        break;
      }
      else if(newton_step==itmax)
      {
        solver_timer.stop();
        pcout << "WARNING: No convergence of Newton's method after " << newton_step << " steps." << std::endl;
        
        break;
      }
    }
  }
  
  
  
  // @sect4{GelfandProblem::compute_solution_norm}
  
  // The computation of the H1-seminorm of the solution can be done in the same
  // way as in step-59. We update the ghost values and use the function
  // integrate_difference from the VectorTools namespace. In the end we gather
  // all computations from all MPI ranks and return the norm.
  template <int dim>
  double
  GelfandProblem<dim>::compute_solution_norm() const 
  {
    solution.update_ghost_values();
    
    Vector<float> norm_per_cell(triangulation.n_active_cells());
    
    VectorTools::integrate_difference(mapping, dof_handler, solution, Functions::ZeroFunction<dim>(), norm_per_cell, QGauss<dim>(fe.degree+2), VectorTools::H1_seminorm);
    
    return VectorTools::compute_global_error(triangulation, norm_per_cell, VectorTools::H1_seminorm);
  }
  
  
  
  // @sect4{GelfandProblem::output_results}
  
  // The graphical output files in vtu format together with an pvtu master file
  // we generate in the same way as in step-37. Note, that we write also the
  // distribution of the triangulation into the output file as it was done in
  // step-40. So no further comments are necessary.
  template <int dim>
  void
  GelfandProblem<dim>::output_results(const unsigned int cycle) const
  {
    if(triangulation.n_global_active_cells() > 1e6)
      return;
    
    solution.update_ghost_values();
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    
    Vector<float> subdomain(triangulation.n_active_cells());
    for(unsigned int i=0; i<subdomain.size(); ++i)
    {
      subdomain(i) = triangulation.locally_owned_subdomain();
    }
    data_out.add_data_vector(subdomain, "subdomain");
    
    data_out.build_patches(fe.degree);
    std::ofstream output("solution-" + Utilities::to_string(cycle, 2) + "." + Utilities::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".vtu");
    data_out.write_vtu(output);
    
    if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
    {
      std::vector<std::string> filenames;
      for(unsigned int i=0; i<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
      {
        filenames.emplace_back("solution-" + Utilities::to_string(cycle, 2) + "." + Utilities::to_string(i, 4) + ".vtu");
      }
      std::ofstream master_output("solution-" + Utilities::to_string(cycle, 2) + ".pvtu");
      data_out.write_pvtu_record(master_output, filenames);
    }
  }
  
  
  
  // @sect4{Run function}
  
  // The last missing function of the solver class for the
  // <i>Gelfand problem</i> is the run function. In the beginning we write a
  // short information of the system and the finite element space we use. The
  // problem is solved several times on a successivley refined mesh.
  template <int dim>
  void
  GelfandProblem<dim>::run() 
  {
    {
      const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      const unsigned int n_vect_doubles = VectorizedArray<double>::n_array_elements;
      const unsigned int n_vect_bits = 8*sizeof(double)*n_vect_doubles;
      
      std::string DAT_header = "START DATE: " + Utilities::System::get_date() + ", TIME: " + Utilities::System::get_time();
      std::string MPI_header = "Running with " + std::to_string(n_ranks) + " MPI process" + (n_ranks>1 ? "es" : "");
      std::string VEC_header = "Vectorization over " + std::to_string(n_vect_doubles) + " doubles = " + std::to_string(n_vect_bits) + " bits (" + Utilities::System::get_current_vectorization_level() + "), VECTORIZATION_LEVEL=" + std::to_string(DEAL_II_COMPILER_VECTORIZATION_LEVEL);
      std::string SOL_header = "Finite element space: " + fe.get_name();
      
      pcout << std::string(80, '=') << std::endl;
      pcout << DAT_header << std::endl;
      pcout << std::string(80, '-') << std::endl;
      
      pcout << MPI_header << std::endl;
      pcout << VEC_header << std::endl;
      pcout << SOL_header << std::endl;
      
      pcout << std::string(80, '=') << std::endl;
    }
    
    
    for(unsigned int cycle=0; cycle<9-dim; ++cycle)
    {
      pcout << std::string(80,'-') << std::endl;
      pcout << "Cycle " << cycle << std::endl;
      pcout << std::string(80,'-') << std::endl;
      
      
      // The first task to actually solve the problem is to generate the
      // triangulation or to refine the triangulation existing one.
      if(cycle==0)
      {
        make_grid();
      }
      else
      {
        triangulation.refine_global(1);
        
      }
      
      
      // Now we setup the system and solve the problem. These steps are
      // accompayned by a time measurement and textual output.
      Timer timer;
      
      pcout << "Setup system..." << std::endl;
      setup_system();
      
      pcout << "   Triangulation: "
      << triangulation.n_global_active_cells() << " cells" << std::endl;
      pcout << "   DoFHandler:    "
      << dof_handler.n_dofs() << " DoFs" << std::endl;
      pcout << std::endl;
      
      
      pcout << "Solve using Newton's method..." << std::endl;
      solve();
      pcout << std::endl;
      
      
      timer.stop();
      pcout << "Time for setup+solve (CPU/Wall) "
      << timer.cpu_time() << "/" << timer.wall_time() << " s" << std::endl;
      pcout << std::endl;
      
      
      // After the problem was solved we compute the norm of the solution and
      // generate the graphical output files.
      pcout << "Output results..." << std::endl;
      const double norm = compute_solution_norm();
      output_results(cycle);
      
      pcout << "  H1 seminorm: " << norm << std::endl;
      pcout << std::endl;
      
    }
    
    // At the end of the run function we close the textual output writing the
    // end time and date.
    {
      pcout << std::string(80, '=') << std::endl;
      pcout << "END DATE: " << Utilities::System::get_date()
      << ", TIME: " << Utilities::System::get_time() << std::endl;
      pcout << std::string(80, '=') << std::endl;
    }
  }
}
