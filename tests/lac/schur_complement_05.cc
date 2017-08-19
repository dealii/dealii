/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2017 by the deal.II authors
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
 */

//
// Check the operation of the Schur Complement for Trilinos vectors
// and in parallel
// This is a lightly modified version of the test "mpi/periodicity_02.cc"
// The only difference is the linear solver routine, which now uses the
// schur complement operation first tested in "lac/schur_complement_03.cc"
//

#define PERIODIC

#include "../tests.h"
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/trilinos_linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/schur_complement.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

namespace Step22
{
  using namespace dealii;

  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem (const unsigned int degree);
    void run ();

  private:
    void setup_dofs ();
    void assemble_system ();
    void solve ();
    void get_point_value (const Point<dim> point, const int proc,
                          Vector<double> &value) const;
    void check_periodicity(const unsigned int cycle) const;
    void output_results (const unsigned int refinement_cycle) const;
    void refine_mesh ();

    const unsigned int   degree;

    MPI_Comm                                    mpi_communicator;

    HyperShellBoundary<dim>                     boundary;
    parallel::distributed::Triangulation<dim>   triangulation;
    FESystem<dim>                               fe;
    DoFHandler<dim>                             dof_handler;

    ConstraintMatrix                            constraints;
    std::vector<IndexSet>                       owned_partitioning;
    std::vector<IndexSet>                       relevant_partitioning;

    TrilinosWrappers::BlockSparsityPattern      sparsity_pattern;
    TrilinosWrappers::BlockSparseMatrix         system_matrix;

    TrilinosWrappers::MPI::BlockVector          solution;
    TrilinosWrappers::MPI::BlockVector          system_rhs;

    ConditionalOStream                          pcout;
  };



  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues () : Function<dim>(dim+1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
  };


  template <int dim>
  double
  BoundaryValues<dim>::value (const Point<dim>  &p,
                              const unsigned int component) const
  {
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    return (1-2*(component==0))*p[(component+1)%2]/(p[0]*p[0]+p[1]*p[1]);
  }


  template <int dim>
  void
  BoundaryValues<dim>::vector_value (const Point<dim> &p,
                                     Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = BoundaryValues<dim>::value (p, c);
  }



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>(dim+1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;

  };


  template <int dim>
  double
  RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                             const unsigned int /*component*/) const
  {
    return 0;
  }


  template <int dim>
  void
  RightHandSide<dim>::vector_value (const Point<dim> &p,
                                    Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = RightHandSide<dim>::value (p, c);
  }





  template <class Matrix, class Preconditioner>
  class InverseMatrix : public Preconditioner
  {
  public:
    InverseMatrix (const Matrix         &m,
                   const Preconditioner &preconditioner,
                   const IndexSet       &locally_owned,
                   const MPI_Comm       &mpi_communicator);

    void vmult (TrilinosWrappers::MPI::Vector       &dst,
                const TrilinosWrappers::MPI::Vector &src) const;

  private:
    const SmartPointer<const Matrix> matrix;
    const SmartPointer<const Preconditioner> preconditioner;

    const MPI_Comm *mpi_communicator;
    mutable TrilinosWrappers::MPI::Vector tmp;
  };


  template <class Matrix, class Preconditioner>
  InverseMatrix<Matrix,Preconditioner>::InverseMatrix
  (const Matrix         &m,
   const Preconditioner &preconditioner,
   const IndexSet       &locally_owned,
   const MPI_Comm       &mpi_communicator)
    :
    matrix (&m),
    preconditioner (&preconditioner),
    mpi_communicator (&mpi_communicator),
    tmp(locally_owned, mpi_communicator)
  {}



  template <class Matrix, class Preconditioner>
  void InverseMatrix<Matrix,Preconditioner>::vmult
  (TrilinosWrappers::MPI::Vector       &dst,
   const TrilinosWrappers::MPI::Vector &src) const
  {
    SolverControl solver_control (src.size(), 1e-6*src.l2_norm(), false, false);
    TrilinosWrappers::SolverCG    cg (solver_control,
                                      TrilinosWrappers::SolverCG::AdditionalData());

    tmp = 0.;
    cg.solve (*matrix, tmp, src, *preconditioner);
    dst = tmp;
  }



  template <class Preconditioner>
  class SchurComplement : public TrilinosWrappers::SparseMatrix
  {
  public:
    SchurComplement ( const TrilinosWrappers::BlockSparseMatrix &system_matrix,
                      const InverseMatrix<TrilinosWrappers::SparseMatrix,
                      Preconditioner> &A_inverse,
                      const IndexSet &owned_pres,
                      const IndexSet &relevant_pres,
                      const MPI_Comm &mpi_communicator);

    void vmult (TrilinosWrappers::MPI::Vector       &dst,
                const TrilinosWrappers::MPI::Vector &src) const;

  private:
    const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> system_matrix;
    const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix,
          Preconditioner> > A_inverse;
    mutable TrilinosWrappers::MPI::Vector tmp1, tmp2;
  };



  template <class Preconditioner>
  SchurComplement<Preconditioner>::
  SchurComplement (const TrilinosWrappers::BlockSparseMatrix &system_matrix,
                   const InverseMatrix<TrilinosWrappers::SparseMatrix,
                   Preconditioner> &A_inverse,
                   const IndexSet &owned_vel,
                   const IndexSet &relevant_vel,
                   const MPI_Comm &mpi_communicator)
    :
    system_matrix (&system_matrix),
    A_inverse (&A_inverse),
    tmp1 (owned_vel, mpi_communicator),
    tmp2 (tmp1)
  {}


  template <class Preconditioner>
  void SchurComplement<Preconditioner>::vmult
  (TrilinosWrappers::MPI::Vector       &dst,
   const TrilinosWrappers::MPI::Vector &src) const
  {
    system_matrix->block(0,1).vmult (tmp1, src);
    A_inverse->vmult (tmp2, tmp1);
    system_matrix->block(1,0).vmult (dst, tmp2);
  }




  template <int dim>
  StokesProblem<dim>::StokesProblem (const unsigned int degree)
    :
    degree (degree),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator/*,
                   Triangulation<dim>::maximum_smoothing*/),
    fe (FE_Q<dim>(degree+1), dim,
        FE_Q<dim>(degree), 1),
    dof_handler (triangulation),
    pcout (Utilities::MPI::this_mpi_process(mpi_communicator)
           == 0
           ?
           deallog.get_file_stream()
           :
           std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator)
            == 0))
  {}




  template <int dim>
  void StokesProblem<dim>::setup_dofs ()
  {
    dof_handler.distribute_dofs (fe);

    std::vector<unsigned int> block_component (dim+1,0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise (dof_handler, block_component);

    std::vector<types::global_dof_index> dofs_per_block (2);
    DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
    const unsigned int n_u = dofs_per_block[0],
                       n_p = dofs_per_block[1];

    {
      owned_partitioning.clear();
      IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
      owned_partitioning.push_back(locally_owned_dofs.get_view(0, n_u));
      owned_partitioning.push_back(locally_owned_dofs.get_view(n_u, n_u+n_p));

      relevant_partitioning.clear();
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
      relevant_partitioning.push_back(locally_relevant_dofs.get_view(0, n_u));
      relevant_partitioning.push_back(locally_relevant_dofs.get_view(n_u, n_u+n_p));

      constraints.clear ();
      constraints.reinit(locally_relevant_dofs);

      FEValuesExtractors::Vector velocities(0);
      FEValuesExtractors::Scalar pressure(dim);

      DoFTools::make_hanging_node_constraints (dof_handler,
                                               constraints);
#ifdef PERIODIC
      VectorTools::interpolate_boundary_values (dof_handler,
                                                3,
                                                BoundaryValues<dim>(),
                                                constraints,
                                                fe.component_mask(velocities));
#endif
      VectorTools::interpolate_boundary_values (dof_handler,
                                                0,
                                                BoundaryValues<dim>(),//Functions::ZeroFunction<dim>(dim+1),
                                                constraints,
                                                fe.component_mask(velocities));
      VectorTools::interpolate_boundary_values (dof_handler,
                                                1,
                                                BoundaryValues<dim>(),//Functions::ZeroFunction<dim>(dim+1),
                                                constraints,
                                                fe.component_mask(velocities));

#ifdef PERIODIC
      std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> >
      periodicity_vector;

      FullMatrix<double> matrix(dim);
      matrix[0][1]=1.;
      matrix[1][0]=-1.;

      std::vector<unsigned int> first_vector_components;
      first_vector_components.push_back(0);

      GridTools::collect_periodic_faces(
        dof_handler, 2, 3, 1, periodicity_vector, Tensor<1, dim>(), matrix);

      DoFTools::make_periodicity_constraints<DoFHandler<dim> >(
        periodicity_vector, constraints, fe.component_mask(velocities),
        first_vector_components);
#endif
    }

    constraints.close ();

    {
      TrilinosWrappers::BlockSparsityPattern bsp
      (owned_partitioning, owned_partitioning,
       relevant_partitioning, mpi_communicator);

      DoFTools::make_sparsity_pattern
      (dof_handler, bsp, constraints, false,
       Utilities::MPI::this_mpi_process(mpi_communicator));

      bsp.compress();

      system_matrix.reinit (bsp);
    }

    system_rhs.reinit (owned_partitioning,
                       mpi_communicator);
    solution.reinit (owned_partitioning, relevant_partitioning,
                     mpi_communicator);
  }



  template <int dim>
  void StokesProblem<dim>::assemble_system ()
  {
    system_matrix=0.;
    system_rhs=0.;

    QGauss<dim>   quadrature_formula(degree+2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |
                             update_quadrature_points  |
                             update_JxW_values |
                             update_gradients);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;

    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const RightHandSide<dim>          right_hand_side;
    std::vector<Vector<double> >      rhs_values (n_q_points,
                                                  Vector<double>(dim+1));

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
    std::vector<double>                  div_phi_u   (dofs_per_cell);
    std::vector<double>                  phi_p       (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          local_matrix = 0;
          local_rhs = 0;

          right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                            rhs_values);

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              for (unsigned int k=0; k<dofs_per_cell; ++k)
                {
                  symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
                  div_phi_u[k]     = fe_values[velocities].divergence (k, q);
                  phi_p[k]         = fe_values[pressure].value (k, q);
                }

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<=i; ++j)
                    {
                      local_matrix(i,j) += (symgrad_phi_u[i] * symgrad_phi_u[j]
                                            - div_phi_u[i] * phi_p[j]
                                            - phi_p[i] * div_phi_u[j]
                                            + phi_p[i] * phi_p[j])
                                           * fe_values.JxW(q);

                    }

                  const unsigned int component_i =
                    fe.system_to_component_index(i).first;
                  local_rhs(i) += fe_values.shape_value(i,q) *
                                  rhs_values[q](component_i) *
                                  fe_values.JxW(q);
                }
            }


          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=i+1; j<dofs_per_cell; ++j)
              local_matrix(i,j) = local_matrix(j,i);

          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                  local_dof_indices,
                                                  system_matrix, system_rhs);
        }

    system_matrix.compress (VectorOperation::add);
    system_rhs.compress (VectorOperation::add);
  }




  template <int dim>
  void StokesProblem<dim>::solve ()
  {
    // Compute a reference solution
    TrilinosWrappers::MPI::BlockVector solution_reference (owned_partitioning,
                                                           relevant_partitioning,
                                                           mpi_communicator);
    {
      deallog.push("ManualSchur");

      TrilinosWrappers::PreconditionJacobi A_preconditioner;
      A_preconditioner.initialize(system_matrix.block(0,0));

      const InverseMatrix<TrilinosWrappers::SparseMatrix,
            TrilinosWrappers::PreconditionJacobi>
            A_inverse (system_matrix.block(0,0),
                       A_preconditioner,
                       owned_partitioning[0],
                       mpi_communicator);

      TrilinosWrappers::MPI::BlockVector tmp (owned_partitioning,
                                              mpi_communicator);

      {
        TrilinosWrappers::MPI::Vector schur_rhs (owned_partitioning[1],
                                                 mpi_communicator);
        A_inverse.vmult (tmp.block(0), system_rhs.block(0));
        system_matrix.block(1,0).vmult (schur_rhs, tmp.block(0));
        schur_rhs -= system_rhs.block(1);

        SchurComplement<TrilinosWrappers::PreconditionJacobi>
        schur_complement (system_matrix, A_inverse,
                          owned_partitioning[0], relevant_partitioning[0],
                          mpi_communicator);

        SolverControl solver_control (solution_reference.block(1).size(),
                                      1e-6*schur_rhs.l2_norm());
        //       TrilinosWrappers::SolverCG    cg (solver_control,
        //                                         TrilinosWrappers::SolverCG::AdditionalData(true));
        SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control);

        TrilinosWrappers::PreconditionAMG preconditioner;
        preconditioner.initialize (system_matrix.block(1,1));

        check_solver_within_range(
          cg.solve (schur_complement, tmp.block(1), schur_rhs, preconditioner),
          solver_control.last_step(), 10, 12);

        constraints.distribute (tmp);
        solution_reference.block(1)=tmp.block(1);
      }

      {
        system_matrix.block(0,1).vmult (tmp.block(0), tmp.block(1));
        tmp.block(0) *= -1;
        tmp.block(0) += system_rhs.block(0);

        A_inverse.vmult (tmp.block(0), tmp.block(0));

        constraints.distribute (tmp);
        solution_reference.block(0)=tmp.block(0);
      }

      deallog.pop();
    }

    // Solution using SchurComplement operator
    {
      deallog.push("SchurComplementOp");

      // Near duplicate of solver routine used in lac/schur_complement_03.cc
      typedef TrilinosWrappers::MPI::Vector VectorType;
      typedef TrilinosWrappers::PreconditionJacobi PreconditionerAType;
      typedef TrilinosWrappers::PreconditionJacobi PreconditionerMType; // Solves in ~14 steps
      // typedef TrilinosWrappers::PreconditionAMG    PreconditionerMType; // Solves in 1 step
      typedef TrilinosWrappers::SolverCG SolverType;

      // Linear operators
      // Also note that this call to linear_operator is
      // actually to TrilinosWrappers::linear_operator
      const auto A  = linear_operator<VectorType>(system_matrix.block(0,0));
      const auto B  = linear_operator<VectorType>(system_matrix.block(0,1));
      const auto C  = linear_operator<VectorType>(system_matrix.block(1,0));
      const auto M  = linear_operator<VectorType>(system_matrix.block(1,1)); // Mass matrix stored in this block
      const auto D0 = null_operator<VectorType>(M);

      // Cheap inverse of A
      PreconditionerAType preconditioner_A;
      preconditioner_A.initialize (system_matrix.block(0,0),
                                   PreconditionerAType::AdditionalData());
      ReductionControl solver_control_A (system_matrix.block(0,0).m(),
                                         1e-12, 1e-6);
      SolverType solver_A (solver_control_A);
      const auto A_inv = inverse_operator(A, solver_A, preconditioner_A);

      // Preconditioner for mass matrix
      PreconditionerMType preconditioner_M;
      preconditioner_M.initialize (system_matrix.block(1,1),
                                   PreconditionerMType::AdditionalData());

      // Schur complement
      const auto S = schur_complement(A_inv,B,C,D0);

      // // Proof of concept for more elaborate preconditioners constructed from
      // // an inverse operator.
      // // This particular implementation is rubbish, but it works. It would be
      // // best to construct a more representative approximation S_inv_approx
      // // and use that instead.
      // // Approximate inverse of A
      // typedef TrilinosWrappers::PreconditionAMG PreconditionerAExpType;
      // PreconditionerAExpType preconditioner_A_exp;
      // preconditioner_A_exp.initialize (system_matrix.block(0,0),
      //                                  PreconditionerAExpType::AdditionalData());
      // const auto A_inv_approx = linear_operator<VectorType>(system_matrix.block(0,0), preconditioner_A_exp);
      // // Approximate Schur complement
      // const auto S_approx = schur_complement<VectorType>(A_inv_approx,B,C,D0);
      // IterationNumberControl solver_control_S_approx (1, 1e-12); // system_matrix.block(1,1).m()
      // SolverType solver_S_approx (solver_control_S_approx);
      // const auto S_inv_approx = inverse_operator(S_approx,solver_S_approx,preconditioner_M);

      // Inverse of Schur complement
      ReductionControl solver_control_S (system_matrix.block(1,1).m(),
                                         1e-9, 1e-6);
      SolverType solver_S (solver_control_S);
      const auto S_inv = inverse_operator(S,solver_S,preconditioner_M);

      TrilinosWrappers::MPI::BlockVector solution_tmp;
      solution_tmp.reinit (owned_partitioning,
                           mpi_communicator);
      VectorType &x = solution_tmp.block(0);
      VectorType &y = solution_tmp.block(1);

      const VectorType &f = system_rhs.block(0);
      const VectorType &g = system_rhs.block(1);

      {
        const unsigned int previous_depth = deallog.depth_file(0);

        const VectorType rhs = condense_schur_rhs (A_inv,C,f,g);
        y = S_inv * rhs;
        x = postprocess_schur_solution (A_inv,B,y,f);

        deallog.depth_file(previous_depth);
      }


      constraints.distribute (solution_tmp);
      solution = solution_tmp;

      deallog.pop();
    }

    // Check output of reference solution
    // solution = solution_reference;

    // The solutions will be slightly different since the solvers are
    // configured differently. So we compute the relative error between
    // them and expect them to remain within a reasonable margin of
    // one another
    // Similar goes for the output of the solution at points evaluated
    // in check_periodicity: The solution here will not necessarily coincide
    // exactly with the parent version of this test
    TrilinosWrappers::MPI::BlockVector local_soln_ref(owned_partitioning,
                                                      mpi_communicator);
    local_soln_ref = solution_reference;
    const double ref_local_norm = local_soln_ref.l2_norm();
    const double ref_global_norm = std::sqrt(Utilities::MPI::sum(
                                               ref_local_norm * ref_local_norm,
                                               mpi_communicator));

    TrilinosWrappers::MPI::BlockVector soln_diff(owned_partitioning,
                                                 mpi_communicator);
    soln_diff = solution;
    soln_diff -= local_soln_ref;
    const double local_error = soln_diff.l2_norm();
    const double global_error = std::sqrt(Utilities::MPI::sum(
                                            local_error * local_error,
                                            mpi_communicator));

    if (global_error/ref_global_norm < 1.e-3)
      {
        deallog << "Relative difference of results less than 1.e-3"
                << std::endl;

      }
    else
      {
        deallog
            << "Relative difference of results too large: "
            << global_error/ref_global_norm
            << std::endl;

        AssertThrow(false,
                    ExcMessage("Solution does not match reference result"));
      }
  }

  template <int dim>
  void StokesProblem<dim>::get_point_value
  (const Point<dim> point, const int proc, Vector<double> &value) const
  {
    typename DoFHandler<dim>::active_cell_iterator cell
      = GridTools::find_active_cell_around_point (dof_handler, point);

    if (cell->is_locally_owned())
      VectorTools::point_value (dof_handler, solution,
                                point, value);

    std::vector<double> tmp (value.size());
    for (unsigned int i=0; i<value.size(); ++i)
      tmp[i]=value[i];

    std::vector<double> tmp2 (value.size());
    MPI_Reduce(&(tmp[0]), &(tmp2[0]), value.size(), MPI_DOUBLE,
               MPI_SUM, proc, mpi_communicator);

    for (unsigned int i=0; i<value.size(); ++i)
      value[i]=tmp2[i];
  }

  template <int dim>
  void StokesProblem<dim>::check_periodicity (const unsigned int cycle) const
  {}

  template <>
  void StokesProblem<2>::check_periodicity (const unsigned int cycle) const
  {
    unsigned int n_points = 4;
    for (unsigned int i = 0; i<cycle; i++)
      n_points*=2;

    //don't test exactly at the support points, since point_value is not stable there
    const double eps = 1./(16.*n_points);

    for (unsigned int i=1; i< n_points; i++)
      {
        Vector<double> value1(3);
        Vector<double> value2(3);

        Point<2> point1;
        point1(0)=0;
        point1(1)=.5*(1.+1.*i/n_points+eps);
        Point<2> point2;
        point2(0)=.5*(1.+1.*i/n_points+eps);
        point2(1)=0.;

        get_point_value (point1, 0, value1);
        get_point_value (point2, 0, value2);

        if (Utilities::MPI::this_mpi_process(mpi_communicator)==0)
          {
            pcout << point1 << "\t" << value1[0] << "\t" << value1[1] << std::endl;
            if (std::abs(value2[0]-value1[1])>1e-8)
              {
                std::cout<<point1<< "\t" << value1[1] << std::endl;
                std::cout<<point2<< "\t" << value2[0] << std::endl;
                Assert(false, ExcInternalError());
              }
            if (std::abs(value2[1]+value1[0])>1e-8)
              {
                std::cout<<point1<< "\t" << value1[0] << std::endl;
                std::cout<<point2<< "\t" << value2[1] << std::endl;
                Assert(false, ExcInternalError());
              }
          }
      }
  }


  template <int dim>
  void
  StokesProblem<dim>::output_results (const unsigned int refinement_cycle)  const
  {
    std::vector<std::string> solution_names (dim, "velocity");
    solution_names.push_back ("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, solution_names,
                              DataOut<dim>::type_dof_data,
                              data_component_interpretation);
    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");
    data_out.build_patches (degree+1);

    std::ostringstream filename;
    filename << "solution-"
             << Utilities::int_to_string (refinement_cycle, 2)
             << "."
             << Utilities::int_to_string (triangulation.locally_owned_subdomain(),2)
             << ".vtu";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtu (output);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
          filenames.push_back (std::string("solution-") +
                               Utilities::int_to_string (refinement_cycle, 2) +
                               "." +
                               Utilities::int_to_string(i, 2) +
                               ".vtu");
        const std::string
        pvtu_master_filename = ("solution-" +
                                Utilities::int_to_string (refinement_cycle, 2) +
                                ".pvtu");
        std::ofstream pvtu_master (pvtu_master_filename.c_str());
        data_out.write_pvtu_record (pvtu_master, filenames);
      }
  }



  template <int dim>
  void
  StokesProblem<dim>::refine_mesh ()
  {

    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    FEValuesExtractors::Scalar pressure(dim);
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(degree+1),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell,
                                        fe.component_mask(pressure));

    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_number (triangulation,
                                     estimated_error_per_cell,
                                     0.3, 0.0);
    triangulation.execute_coarsening_and_refinement ();
  }



  template <int dim>
  void StokesProblem<dim>::run ()
  {
    Point<dim> center;
    const double inner_radius = .5;
    const double outer_radius = 1.;

    GridGenerator::quarter_hyper_shell (triangulation,
                                        center,
                                        inner_radius,
                                        outer_radius,
                                        0,
                                        true);

#ifdef PERIODIC
    std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> >
    periodicity_vector;

    FullMatrix<double> matrix(dim);
    matrix[0][1]=1.;
    matrix[1][0]=-1.;

    std::vector<unsigned int> first_vector_components;
    first_vector_components.push_back(0);

    GridTools::collect_periodic_faces(
      triangulation, 2, 3, 1, periodicity_vector, Tensor<1, dim>(), matrix);

    triangulation.add_periodicity(periodicity_vector);
#endif


    triangulation.set_boundary(0, boundary);
    triangulation.set_boundary(1, boundary);

    triangulation.refine_global (4-dim);

    for (unsigned int refinement_cycle = 0; refinement_cycle<3;
         ++refinement_cycle)
      {
        pcout << "Refinement cycle " << refinement_cycle << std::endl;

        if (refinement_cycle > 0)
          refine_mesh ();

        setup_dofs ();

        pcout << "   Assembling..." << std::endl << std::flush;
        assemble_system ();

        pcout << "   Solving..." << std::endl << std::flush;
        solve ();

#ifdef PERIODIC
        check_periodicity(refinement_cycle);
#endif
        // output_results (refinement_cycle);

        pcout << std::endl;
      }
  }
}



int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step22;

      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)==0)
        {
          std::ofstream logfile("output");
          deallog.attach(logfile, false);
          {
            StokesProblem<2> flow_problem(1);
            flow_problem.run ();
          }
        }
      else
        {
          StokesProblem<2> flow_problem(1);
          flow_problem.run ();
        }
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
