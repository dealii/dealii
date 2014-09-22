// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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



// a hp-ified version of step-13


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
std::ofstream logfile("output");


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <iostream>
#include <fstream>
#include <list>
#include <sstream>



namespace Evaluation
{

  template <int dim>
  class EvaluationBase
  {
  public:
    virtual ~EvaluationBase ();

    void set_refinement_cycle (const unsigned int refinement_cycle);

    virtual void operator () (const hp::DoFHandler<dim> &dof_handler,
                              const Vector<double>  &solution) const = 0;
  protected:
    unsigned int refinement_cycle;
  };


  template <int dim>
  EvaluationBase<dim>::~EvaluationBase ()
  {}



  template <int dim>
  void
  EvaluationBase<dim>::set_refinement_cycle (const unsigned int step)
  {
    refinement_cycle = step;
  }



  template <int dim>
  class PointValueEvaluation : public EvaluationBase<dim>
  {
  public:
    PointValueEvaluation (const Point<dim>   &evaluation_point,
                          TableHandler       &results_table);

    virtual void operator () (const hp::DoFHandler<dim> &dof_handler,
                              const Vector<double>  &solution) const;

    DeclException1 (ExcEvaluationPointNotFound,
                    Point<dim>,
                    << "The evaluation point " << arg1
                    << " was not found among the vertices of the present grid.");
  private:
    const Point<dim>  evaluation_point;
    TableHandler     &results_table;
  };


  template <int dim>
  PointValueEvaluation<dim>::
  PointValueEvaluation (const Point<dim>   &evaluation_point,
                        TableHandler       &results_table)
    :
    evaluation_point (evaluation_point),
    results_table (results_table)
  {}



  template <int dim>
  void
  PointValueEvaluation<dim>::
  operator () (const hp::DoFHandler<dim> &dof_handler,
               const Vector<double>  &solution) const
  {
    double point_value = 1e20;

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    bool evaluation_point_found = false;
    for (; (cell!=endc) && !evaluation_point_found; ++cell)
      for (unsigned int vertex=0;
           vertex<GeometryInfo<dim>::vertices_per_cell;
           ++vertex)
        if (cell->vertex(vertex) == evaluation_point)
          {
            point_value = solution(cell->vertex_dof_index(vertex,0,
                                                          cell->active_fe_index()));

            evaluation_point_found = true;
            break;
          };

    AssertThrow (evaluation_point_found,
                 ExcEvaluationPointNotFound(evaluation_point));

    results_table.add_value ("DoFs", dof_handler.n_dofs());
    results_table.add_value ("u(x_0)", point_value);
  }





  template <int dim>
  class SolutionOutput : public EvaluationBase<dim>
  {
  public:
    SolutionOutput (const std::string                         &output_name_base,
                    const DataOutBase::OutputFormat  output_format);

    virtual void operator () (const hp::DoFHandler<dim> &dof_handler,
                              const Vector<double>  &solution) const;
  private:
    const std::string                         output_name_base;
    const DataOutBase::OutputFormat output_format;
  };


  template <int dim>
  SolutionOutput<dim>::
  SolutionOutput (const std::string                         &output_name_base,
                  const DataOutBase::OutputFormat  output_format)
    :
    output_name_base (output_name_base),
    output_format (output_format)
  {}


  template <int dim>
  void
  SolutionOutput<dim>::operator () (const hp::DoFHandler<dim> &dof_handler,
                                    const Vector<double>  &solution) const
  {
    DataOut<dim,hp::DoFHandler<dim> > data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.build_patches ();

    std::ostringstream filename;
    filename << output_name_base << "-"
             << this->refinement_cycle
             << data_out.default_suffix (output_format)
             << std::ends;

    data_out.write (deallog.get_file_stream(), output_format);
  }




}



namespace LaplaceSolver
{

  template <int dim>
  class Base
  {
  public:
    Base (Triangulation<dim> &coarse_grid);
    virtual ~Base ();

    virtual void solve_problem () = 0;
    virtual void postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const = 0;
    virtual void refine_grid () = 0;
    virtual unsigned int n_dofs () const = 0;

  protected:
    const SmartPointer<Triangulation<dim> > triangulation;
  };


  template <int dim>
  Base<dim>::Base (Triangulation<dim> &coarse_grid)
    :
    triangulation (&coarse_grid)
  {}


  template <int dim>
  Base<dim>::~Base ()
  {}



  template <int dim>
  class Solver : public virtual Base<dim>
  {
  public:
    Solver (Triangulation<dim>       &triangulation,
            const hp::FECollection<dim> &fe,
            const hp::QCollection<dim>    &quadrature,
            const Function<dim>      &boundary_values);
    virtual
    ~Solver ();

    virtual
    void
    solve_problem ();

    virtual
    void
    postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

    virtual
    unsigned int
    n_dofs () const;

  protected:
    const SmartPointer<const hp::FECollection<dim> >  fe;
    const SmartPointer<const hp::QCollection<dim> >     quadrature;
    hp::DoFHandler<dim>                                dof_handler;
    Vector<double>                                 solution;
    const SmartPointer<const Function<dim> >       boundary_values;

    virtual void assemble_rhs (Vector<double> &rhs) const = 0;

  private:
    struct LinearSystem
    {
      LinearSystem (const hp::DoFHandler<dim> &dof_handler);

      void solve (Vector<double> &solution) const;

      ConstraintMatrix     hanging_node_constraints;
      SparsityPattern      sparsity_pattern;
      SparseMatrix<double> matrix;
      Vector<double>       rhs;
    };

    void
    assemble_linear_system (LinearSystem &linear_system);

    void
    assemble_matrix (LinearSystem                                         &linear_system,
                     const typename hp::DoFHandler<dim>::active_cell_iterator &begin_cell,
                     const typename hp::DoFHandler<dim>::active_cell_iterator &end_cell,
                     Threads::ThreadMutex                                 &mutex) const;
  };



  template <int dim>
  Solver<dim>::Solver (Triangulation<dim>       &triangulation,
                       const hp::FECollection<dim> &fe,
                       const hp::QCollection<dim>    &quadrature,
                       const Function<dim>      &boundary_values)
    :
    Base<dim> (triangulation),
    fe (&fe),
    quadrature (&quadrature),
    dof_handler (triangulation),
    boundary_values (&boundary_values)
  {}


  template <int dim>
  Solver<dim>::~Solver ()
  {
    dof_handler.clear ();
  }


  template <int dim>
  void
  Solver<dim>::solve_problem ()
  {
    dof_handler.distribute_dofs (*fe);
    solution.reinit (dof_handler.n_dofs());

    LinearSystem linear_system (dof_handler);
    assemble_linear_system (linear_system);
    linear_system.solve (solution);
  }


  template <int dim>
  void
  Solver<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    postprocessor (dof_handler, solution);
  }


  template <int dim>
  unsigned int
  Solver<dim>::n_dofs () const
  {
    return dof_handler.n_dofs();
  }


  template <int dim>
  void
  Solver<dim>::assemble_linear_system (LinearSystem &linear_system)
  {
    typedef
    typename hp::DoFHandler<dim>::active_cell_iterator
    active_cell_iterator;

    const unsigned int n_threads = multithread_info.n_default_threads;
    std::vector<std::pair<active_cell_iterator,active_cell_iterator> >
    thread_ranges
      = Threads::split_range<active_cell_iterator> (dof_handler.begin_active (),
                                                    dof_handler.end (),
                                                    n_threads);

    Threads::ThreadMutex mutex;
    Threads::ThreadGroup<> threads;
    for (unsigned int thread=0; thread<n_threads; ++thread)
      threads += Threads::spawn (*this, &Solver<dim>::assemble_matrix)
                 (linear_system,
                  thread_ranges[thread].first,
                  thread_ranges[thread].second,
                  mutex);

    assemble_rhs (linear_system.rhs);
    linear_system.hanging_node_constraints.condense (linear_system.rhs);

    std::map<types::global_dof_index,double> boundary_value_map;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              *boundary_values,
                                              boundary_value_map);


    threads.join_all ();
    linear_system.hanging_node_constraints.condense (linear_system.matrix);

    MatrixTools::apply_boundary_values (boundary_value_map,
                                        linear_system.matrix,
                                        solution,
                                        linear_system.rhs);

  }


  template <int dim>
  void
  Solver<dim>::assemble_matrix (LinearSystem                                         &linear_system,
                                const typename hp::DoFHandler<dim>::active_cell_iterator &begin_cell,
                                const typename hp::DoFHandler<dim>::active_cell_iterator &end_cell,
                                Threads::ThreadMutex                                 &mutex) const
  {
    hp::FEValues<dim> fe_values (*fe, *quadrature,
                                 update_gradients | update_JxW_values);

    const unsigned int   dofs_per_cell = (*fe)[0].dofs_per_cell;
    const unsigned int   n_q_points    = (*quadrature)[0].size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    for (typename hp::DoFHandler<dim>::active_cell_iterator cell=begin_cell;
         cell!=end_cell; ++cell)
      {
        cell_matrix = 0;

        fe_values.reinit (cell);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.get_present_fe_values().shape_grad(i,q_point) *
                                   fe_values.get_present_fe_values().shape_grad(j,q_point) *
                                   fe_values.get_present_fe_values().JxW(q_point));


        cell->get_dof_indices (local_dof_indices);

        Threads::ThreadMutex::ScopedLock lock (mutex);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            linear_system.matrix.add (local_dof_indices[i],
                                      local_dof_indices[j],
                                      cell_matrix(i,j));
      };
  }


  template <int dim>
  Solver<dim>::LinearSystem::
  LinearSystem (const hp::DoFHandler<dim> &dof_handler)
  {
    hanging_node_constraints.clear ();

    void (*mhnc_p) (const hp::DoFHandler<dim> &,
                    ConstraintMatrix &)
      = &DoFTools::make_hanging_node_constraints;

    Threads::Thread<>
    mhnc_thread = Threads::spawn (mhnc_p)(dof_handler,
                                          hanging_node_constraints);

    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

    mhnc_thread.join ();
    hanging_node_constraints.close ();
    hanging_node_constraints.condense (sparsity_pattern);

    sparsity_pattern.compress();
    matrix.reinit (sparsity_pattern);
    rhs.reinit (dof_handler.n_dofs());
  }



  template <int dim>
  void
  Solver<dim>::LinearSystem::solve (Vector<double> &solution) const
  {
    SolverControl           solver_control (1000, 1e-12);
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(matrix, 1.2);

    cg.solve (matrix, solution, rhs, preconditioner);

    hanging_node_constraints.distribute (solution);
  }





  template <int dim>
  class PrimalSolver : public Solver<dim>
  {
  public:
    PrimalSolver (Triangulation<dim>       &triangulation,
                  const hp::FECollection<dim> &fe,
                  const hp::QCollection<dim>    &quadrature,
                  const Function<dim>      &rhs_function,
                  const Function<dim>      &boundary_values);
  protected:
    const SmartPointer<const Function<dim> > rhs_function;
    virtual void assemble_rhs (Vector<double> &rhs) const;
  };


  template <int dim>
  PrimalSolver<dim>::
  PrimalSolver (Triangulation<dim>       &triangulation,
                const hp::FECollection<dim> &fe,
                const hp::QCollection<dim>    &quadrature,
                const Function<dim>      &rhs_function,
                const Function<dim>      &boundary_values)
    :
    Base<dim> (triangulation),
    Solver<dim> (triangulation, fe,
                 quadrature, boundary_values),
    rhs_function (&rhs_function)
  {}



  template <int dim>
  void
  PrimalSolver<dim>::
  assemble_rhs (Vector<double> &rhs) const
  {
    hp::FEValues<dim> fe_values (*this->fe, *this->quadrature,
                                 update_values | update_q_points  |
                                 update_JxW_values);

    const unsigned int   dofs_per_cell = (*this->fe)[0].dofs_per_cell;
    const unsigned int   n_q_points    = (*this->quadrature)[0].size();

    Vector<double>       cell_rhs (dofs_per_cell);
    std::vector<double>  rhs_values (n_q_points);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = this->dof_handler.begin_active(),
    endc = this->dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell_rhs = 0;
        fe_values.reinit (cell);
        rhs_function->value_list (fe_values.get_present_fe_values().get_quadrature_points(),
                                  rhs_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_rhs(i) += (fe_values.get_present_fe_values().shape_value(i,q_point) *
                            rhs_values[q_point] *
                            fe_values.get_present_fe_values().JxW(q_point));

        cell->get_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          rhs(local_dof_indices[i]) += cell_rhs(i);
      };
  }



  template <int dim>
  class RefinementGlobal : public PrimalSolver<dim>
  {
  public:
    RefinementGlobal (Triangulation<dim>       &coarse_grid,
                      const hp::FECollection<dim> &fe,
                      const hp::QCollection<dim>    &quadrature,
                      const Function<dim>      &rhs_function,
                      const Function<dim>      &boundary_values);

    virtual void refine_grid ();
  };



  template <int dim>
  RefinementGlobal<dim>::
  RefinementGlobal (Triangulation<dim>       &coarse_grid,
                    const hp::FECollection<dim> &fe,
                    const hp::QCollection<dim>    &quadrature,
                    const Function<dim>      &rhs_function,
                    const Function<dim>      &boundary_values)
    :
    Base<dim> (coarse_grid),
    PrimalSolver<dim> (coarse_grid, fe, quadrature,
                       rhs_function, boundary_values)
  {}



  template <int dim>
  void
  RefinementGlobal<dim>::refine_grid ()
  {
    this->triangulation->refine_global (1);
  }



  template <int dim>
  class RefinementKelly : public PrimalSolver<dim>
  {
  public:
    RefinementKelly (Triangulation<dim>       &coarse_grid,
                     const hp::FECollection<dim> &fe,
                     const hp::QCollection<dim>    &quadrature,
                     const Function<dim>      &rhs_function,
                     const Function<dim>      &boundary_values);

    virtual void refine_grid ();
  };



  template <int dim>
  RefinementKelly<dim>::
  RefinementKelly (Triangulation<dim>       &coarse_grid,
                   const hp::FECollection<dim> &fe,
                   const hp::QCollection<dim>    &quadrature,
                   const Function<dim>      &rhs_function,
                   const Function<dim>      &boundary_values)
    :
    Base<dim> (coarse_grid),
    PrimalSolver<dim> (coarse_grid, fe, quadrature,
                       rhs_function, boundary_values)
  {}



  template <int dim>
  void
  RefinementKelly<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (this->triangulation->n_active_cells());
    KellyErrorEstimator<dim>::estimate (this->dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        this->solution,
                                        estimated_error_per_cell);
    GridRefinement::refine_and_coarsen_fixed_number (*this->triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.03);
    this->triangulation->execute_coarsening_and_refinement ();
  }

}





template <int dim>
class Solution : public Function<dim>
{
public:
  Solution () : Function<dim> () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const;
};


template <int dim>
double
Solution<dim>::value (const Point<dim>   &p,
                      const unsigned int  /*component*/) const
{
  double q = p(0);
  for (unsigned int i=1; i<dim; ++i)
    q += std::sin(10*p(i)+5*p(0)*p(0));
  const double exponential = std::exp(q);
  return exponential;
}



template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide () : Function<dim> () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>   &p,
                           const unsigned int  /*component*/) const
{
  double q = p(0);
  for (unsigned int i=1; i<dim; ++i)
    q += std::sin(10*p(i)+5*p(0)*p(0));
  const double u = std::exp(q);
  double t1 = 1,
         t2 = 0,
         t3 = 0;
  for (unsigned int i=1; i<dim; ++i)
    {
      t1 += std::cos(10*p(i)+5*p(0)*p(0)) * 10 * p(0);
      t2 += 10*std::cos(10*p(i)+5*p(0)*p(0)) -
            100*std::sin(10*p(i)+5*p(0)*p(0)) * p(0)*p(0);
      t3 += 100*std::cos(10*p(i)+5*p(0)*p(0))*std::cos(10*p(i)+5*p(0)*p(0)) -
            100*std::sin(10*p(i)+5*p(0)*p(0));
    };
  t1 = t1*t1;

  return -u*(t1+t2+t3);
}




template <int dim>
void
run_simulation (LaplaceSolver::Base<dim>                     &solver,
                const std::list<Evaluation::EvaluationBase<dim> *> &postprocessor_list)
{
  deallog << "Refinement cycle: ";

  for (unsigned int step=0; true; ++step)
    {
      deallog << step << " " << std::flush;

      solver.solve_problem ();

      for (typename std::list<Evaluation::EvaluationBase<dim> *>::const_iterator
           i = postprocessor_list.begin();
           i != postprocessor_list.end(); ++i)
        {
          (*i)->set_refinement_cycle (step);
          solver.postprocess (**i);
        };


      if (solver.n_dofs() < 2000)
        solver.refine_grid ();
      else
        break;
    };

  deallog << std::endl;
}



template <int dim>
void solve_problem (const std::string &solver_name)
{
  const std::string header = "Running tests with \"" + solver_name +
                             "\" refinement criterion:";
  deallog << header << std::endl
          << std::string (header.size(), '-') << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (2);
  const hp::FECollection<dim>          fe(FE_Q<dim>(1));
  const hp::QCollection<dim>       quadrature(QGauss<dim>(4));
  const RightHandSide<dim> rhs_function;
  const Solution<dim>      boundary_values;

  LaplaceSolver::Base<dim> *solver = 0;
  if (solver_name == "global")
    solver = new LaplaceSolver::RefinementGlobal<dim> (triangulation, fe,
                                                       quadrature,
                                                       rhs_function,
                                                       boundary_values);
  else if (solver_name == "kelly")
    solver = new LaplaceSolver::RefinementKelly<dim> (triangulation, fe,
                                                      quadrature,
                                                      rhs_function,
                                                      boundary_values);
  else
    AssertThrow (false, ExcNotImplemented());

  TableHandler results_table;
  Evaluation::PointValueEvaluation<dim>
  postprocessor1 (Point<dim>(0.5,0.5), results_table);

  Evaluation::SolutionOutput<dim>
  postprocessor2 (std::string("solution-")+solver_name,
                  DataOutBase::gnuplot);

  std::list<Evaluation::EvaluationBase<dim> *> postprocessor_list;
  postprocessor_list.push_back (&postprocessor1);
  postprocessor_list.push_back (&postprocessor2);

  run_simulation (*solver, postprocessor_list);

  results_table.write_text (deallog.get_file_stream());
  delete solver;

  deallog << std::endl;
}



int main ()
{
  try
    {
      logfile.precision(2);
      deallog << std::setprecision(2);

      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      solve_problem<2> ("global");
      solve_problem<2> ("kelly");
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
    };

  return 0;
}
