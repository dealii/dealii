// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test for correctness of step-37 (without output and only small sizes)


#include "../tests.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/multigrid.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <fstream>
#include <sstream>


namespace Step37
{
  using namespace dealii;

  const unsigned int degree_finite_element = 2;



  template <int dim>
  class Coefficient : public Function<dim>
  {
  public:
    Coefficient ()  : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    template <typename number>
    number value (const Point<dim,number> &p,
                  const unsigned int component = 0) const;

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<double>            &values,
                             const unsigned int              component = 0) const;
  };



  template <int dim>
  template <typename number>
  number Coefficient<dim>::value (const Point<dim,number> &p,
                                  const unsigned int /*component*/) const
  {
    return 1. / (0.05 + 2.*p.square());
  }



  template <int dim>
  double Coefficient<dim>::value (const Point<dim> &p,
                                  const unsigned int component) const
  {
    return value<double>(p,component);
  }



  template <int dim>
  void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
                                     std::vector<double>            &values,
                                     const unsigned int              component) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch (values.size(), points.size()));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    const unsigned int n_points = points.size();
    for (unsigned int i=0; i<n_points; ++i)
      values[i] = value<double>(points[i],component);
  }





  template <int dim, int fe_degree, typename number>
  class LaplaceOperator : public Subscriptor
  {
  public:
    LaplaceOperator ();

    void clear();

    void reinit (const DoFHandler<dim>  &dof_handler,
                 const ConstraintMatrix &constraints,
                 const unsigned int      level = numbers::invalid_unsigned_int);

    unsigned int m () const;
    unsigned int n () const;

    void vmult (Vector<double> &dst,
                const Vector<double> &src) const;
    void Tvmult (Vector<double> &dst,
                 const Vector<double> &src) const;
    void vmult_add (Vector<double> &dst,
                    const Vector<double> &src) const;
    void Tvmult_add (Vector<double> &dst,
                     const Vector<double> &src) const;

    number el (const unsigned int row,
               const unsigned int col) const;
    void set_diagonal (const Vector<number> &diagonal);

    std::size_t memory_consumption () const;

  private:
    void local_apply (const MatrixFree<dim,number>    &data,
                      Vector<double>                      &dst,
                      const Vector<double>                &src,
                      const std::pair<unsigned int,unsigned int> &cell_range) const;

    void evaluate_coefficient(const Coefficient<dim> &function);

    MatrixFree<dim,number>      data;
    AlignedVector<VectorizedArray<number> > coefficient;

    Vector<number>  diagonal_values;
    bool            diagonal_is_available;
  };



  template <int dim, int fe_degree, typename number>
  LaplaceOperator<dim,fe_degree,number>::LaplaceOperator ()
    :
    Subscriptor()
  {}



  template <int dim, int fe_degree, typename number>
  unsigned int
  LaplaceOperator<dim,fe_degree,number>::m () const
  {
    return data.get_vector_partitioner()->size();
  }



  template <int dim, int fe_degree, typename number>
  unsigned int
  LaplaceOperator<dim,fe_degree,number>::n () const
  {
    return data.get_vector_partitioner()->size();
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::clear ()
  {
    data.clear();
    diagonal_is_available = false;
    diagonal_values.reinit(0);
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::reinit (const DoFHandler<dim>  &dof_handler,
                                                 const ConstraintMatrix &constraints,
                                                 const unsigned int      level)
  {
    typename MatrixFree<dim,number>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim,number>::AdditionalData::partition_color;
    additional_data.level_mg_handler = level;
    additional_data.mapping_update_flags = (update_gradients | update_JxW_values |
                                            update_quadrature_points);
    data.reinit (dof_handler, constraints, QGauss<1>(fe_degree+1),
                 additional_data);
    evaluate_coefficient(Coefficient<dim>());
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::
  evaluate_coefficient (const Coefficient<dim> &coefficient_function)
  {
    const unsigned int n_cells = data.n_macro_cells();
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi (data);
    coefficient.resize (n_cells * phi.n_q_points);
    for (unsigned int cell=0; cell<n_cells; ++cell)
      {
        phi.reinit (cell);
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          coefficient[cell*phi.n_q_points+q] =
            coefficient_function.value(phi.quadrature_point(q));
      }
  }




  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::
  local_apply (const MatrixFree<dim,number>         &data,
               Vector<double>                       &dst,
               const Vector<double>                 &src,
               const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi (data);
    AssertDimension (coefficient.size(),
                     data.n_macro_cells() * phi.n_q_points);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);
        phi.read_dof_values(src);
        phi.evaluate (false,true,false);
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          phi.submit_gradient (coefficient[cell*phi.n_q_points+q] *
                               phi.get_gradient(q), q);
        phi.integrate (false,true);
        phi.distribute_local_to_global (dst);
      }
  }




  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::vmult (Vector<double>       &dst,
                                                const Vector<double> &src) const
  {
    dst = 0;
    vmult_add (dst, src);
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::Tvmult (Vector<double>       &dst,
                                                 const Vector<double> &src) const
  {
    dst = 0;
    vmult_add (dst,src);
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::Tvmult_add (Vector<double>       &dst,
                                                     const Vector<double> &src) const
  {
    vmult_add (dst,src);
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::vmult_add (Vector<double>       &dst,
                                                    const Vector<double> &src) const
  {
    data.cell_loop (&LaplaceOperator::local_apply, this, dst, src);

    const std::vector<unsigned int> &
    constrained_dofs = data.get_constrained_dofs();
    for (unsigned int i=0; i<constrained_dofs.size(); ++i)
      dst(constrained_dofs[i]) += src(constrained_dofs[i]);
  }



  template <int dim, int fe_degree, typename number>
  number
  LaplaceOperator<dim,fe_degree,number>::el (const unsigned int row,
                                             const unsigned int col) const
  {
    Assert (row == col, ExcNotImplemented());
    Assert (diagonal_is_available == true, ExcNotInitialized());
    return diagonal_values(row);
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>::set_diagonal(const Vector<number> &diagonal)
  {
    AssertDimension (m(), diagonal.size());

    diagonal_values = diagonal;

    const std::vector<unsigned int> &
    constrained_dofs = data.get_constrained_dofs();
    for (unsigned int i=0; i<constrained_dofs.size(); ++i)
      diagonal_values(constrained_dofs[i]) = 1.0;

    diagonal_is_available = true;
  }



  template <int dim, int fe_degree, typename number>
  std::size_t
  LaplaceOperator<dim,fe_degree,number>::memory_consumption () const
  {
    return (data.memory_consumption () +
            coefficient.memory_consumption() +
            diagonal_values.memory_consumption() +
            MemoryConsumption::memory_consumption(diagonal_is_available));
  }




  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem ();
    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void assemble_multigrid ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    typedef LaplaceOperator<dim,degree_finite_element,double> SystemMatrixType;
    typedef LaplaceOperator<dim,degree_finite_element,float>  LevelMatrixType;

    Triangulation<dim>               triangulation;
    FE_Q<dim>                        fe;
    DoFHandler<dim>                  dof_handler;
    ConstraintMatrix                 constraints;

    SystemMatrixType                 system_matrix;
    MGLevelObject<LevelMatrixType>   mg_matrices;
    FullMatrix<float>                coarse_matrix;
    MGLevelObject<ConstraintMatrix>  mg_constraints;

    Vector<double>                   solution;
    Vector<double>                   system_rhs;
  };



  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem ()
    :
    fe (degree_finite_element),
    dof_handler (triangulation)
  {}




  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    system_matrix.clear();
    mg_matrices.clear();
    mg_constraints.clear();

    dof_handler.distribute_dofs (fe);
    dof_handler.distribute_mg_dofs (fe);

    deallog << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;

    constraints.clear();
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              constraints);
    constraints.close();

    system_matrix.reinit (dof_handler, constraints);
    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    const unsigned int nlevels = triangulation.n_levels();
    mg_matrices.resize(0, nlevels-1);
    mg_constraints.resize (0, nlevels-1);

    typename FunctionMap<dim>::type dirichlet_boundary;
    ZeroFunction<dim>               homogeneous_dirichlet_bc (1);
    dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
    std::vector<std::set<types::global_dof_index> > boundary_indices(triangulation.n_levels());
    MGTools::make_boundary_list (dof_handler,
                                 dirichlet_boundary,
                                 boundary_indices);
    for (unsigned int level=0; level<nlevels; ++level)
      {
        std::set<types::global_dof_index>::iterator bc_it = boundary_indices[level].begin();
        for ( ; bc_it != boundary_indices[level].end(); ++bc_it)
          mg_constraints[level].add_line(*bc_it);

        mg_constraints[level].close();
        mg_matrices[level].reinit(dof_handler,
                                  mg_constraints[level],
                                  level);
      }
    coarse_matrix.reinit (dof_handler.n_dofs(0),
                          dof_handler.n_dofs(0));
  }




  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    const Coefficient<dim> coefficient;
    std::vector<double>    coefficient_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell->get_dof_indices (local_dof_indices);
        fe_values.reinit (cell);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            double rhs_val = 0;
            for (unsigned int q=0; q<n_q_points; ++q)
              rhs_val += (fe_values.shape_value(i,q) * 1.0 *
                          fe_values.JxW(q));
            system_rhs(local_dof_indices[i]) += rhs_val;
          }
      }
    constraints.condense(system_rhs);
  }



  template <int dim>
  void LaplaceProblem<dim>::assemble_multigrid ()
  {
    coarse_matrix = 0;
    QGauss<dim>  quadrature_formula(fe.degree+1);
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_gradients  | update_inverse_jacobians |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    const Coefficient<dim>    coefficient;
    std::vector<double>       coefficient_values (n_q_points);
    FullMatrix<double>        local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>            local_diagonal (dofs_per_cell);

    const unsigned int n_levels = triangulation.n_levels();
    std::vector<Vector<float> > diagonals (n_levels);
    for (unsigned int level=0; level<n_levels; ++level)
      diagonals[level].reinit (dof_handler.n_dofs(level));

    std::vector<unsigned int> cell_no(triangulation.n_levels());
    typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin(),
                                            endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        const unsigned int level = cell->level();
        cell->get_mg_dof_indices (local_dof_indices);
        fe_values.reinit (cell);
        coefficient.value_list (fe_values.get_quadrature_points(),
                                coefficient_values);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            double local_diag = 0;
            for (unsigned int q=0; q<n_q_points; ++q)
              local_diag += ((fe_values.shape_grad(i,q) *
                              fe_values.shape_grad(i,q)) *
                             coefficient_values[q] * fe_values.JxW(q));
            local_diagonal(i) = local_diag;
          }
        mg_constraints[level].distribute_local_to_global(local_diagonal,
                                                         local_dof_indices,
                                                         diagonals[level]);

        if (level == 0)
          {
            local_matrix = 0;

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                {
                  double add_value = 0;
                  for (unsigned int q=0; q<n_q_points; ++q)
                    add_value += (fe_values.shape_grad(i,q) *
                                  fe_values.shape_grad(j,q) *
                                  coefficient_values[q] *
                                  fe_values.JxW(q));
                  local_matrix(i,j) = add_value;
                }
            mg_constraints[0].distribute_local_to_global (local_matrix,
                                                          local_dof_indices,
                                                          coarse_matrix);
          }
      }

    for (unsigned int level=0; level<n_levels; ++level)
      mg_matrices[level].set_diagonal (diagonals[level]);
  }



  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    MGTransferPrebuilt<Vector<double> > mg_transfer;
    mg_transfer.build_matrices(dof_handler);

    MGCoarseGridHouseholder<float, Vector<double> > mg_coarse;
    mg_coarse.initialize(coarse_matrix);

    typedef PreconditionChebyshev<LevelMatrixType,Vector<double> > SMOOTHER;
    MGSmootherPrecondition<LevelMatrixType, SMOOTHER, Vector<double> >
    mg_smoother;

    typename SMOOTHER::AdditionalData smoother_data;
    smoother_data.smoothing_range = 10.;
    smoother_data.degree = 6;
    smoother_data.eig_cg_n_iterations = 10;
    mg_smoother.initialize(mg_matrices, smoother_data);

    MGMatrix<LevelMatrixType, Vector<double> >
    mg_matrix(&mg_matrices);

    Multigrid<Vector<double> > mg(dof_handler,
                                  mg_matrix,
                                  mg_coarse,
                                  mg_transfer,
                                  mg_smoother,
                                  mg_smoother);
    PreconditionMG<dim, Vector<double>,
                   MGTransferPrebuilt<Vector<double> > >
                   preconditioner(dof_handler, mg, mg_transfer);

    const std::size_t multigrid_memory
      = (mg_matrices.memory_consumption() +
         mg_transfer.memory_consumption() +
         coarse_matrix.memory_consumption());

    SolverControl           solver_control (1000, 1e-12*system_rhs.l2_norm());
    SolverCG<>              cg (solver_control);

    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);
  }



  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<3; ++cycle)
      {
        deallog << "Cycle " << cycle << std::endl;

        if (cycle == 0)
          {
            GridGenerator::hyper_cube (triangulation, 0., 1.);
            triangulation.refine_global (3-dim);
          }
        triangulation.refine_global (1);
        setup_system ();
        assemble_system ();
        assemble_multigrid ();
        solve ();
        deallog << std::endl;
      };
  }
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision (3);
  deallog.threshold_double(1e-10);
  deallog.depth_console(0);

  {
    deallog.push("2d");
    Step37::LaplaceProblem<2> laplace_problem;
    laplace_problem.run();
    deallog.pop();
  }
  {
    deallog.push("3d");
    Step37::LaplaceProblem<3> laplace_problem;
    laplace_problem.run();
    deallog.pop();
  }


  return 0;
}

