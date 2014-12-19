// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// a hp-ified version of step-14


#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <fstream>
std::ofstream logfile("step-14/output");


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <numeric>
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
    PointValueEvaluation (const Point<dim>   &evaluation_point);

    virtual void operator () (const hp::DoFHandler<dim> &dof_handler,
                              const Vector<double>  &solution) const;

    DeclException1 (ExcEvaluationPointNotFound,
                    Point<dim>,
                    << "The evaluation point " << arg1
                    << " was not found among the vertices of the present grid.");
  private:
    const Point<dim>  evaluation_point;
  };


  template <int dim>
  PointValueEvaluation<dim>::
  PointValueEvaluation (const Point<dim>   &evaluation_point)
    :
    evaluation_point (evaluation_point)
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
        if (cell->vertex(vertex).distance (evaluation_point)
            <
            cell->diameter() * 1e-8)
          {
            point_value = solution(cell->vertex_dof_index(vertex,0));

            evaluation_point_found = true;
            break;
          };

    AssertThrow (evaluation_point_found,
                 ExcEvaluationPointNotFound(evaluation_point));

    deallog << "   Point value=" << point_value
            << std::endl;
  }



  template <int dim>
  class PointXDerivativeEvaluation : public EvaluationBase<dim>
  {
  public:
    PointXDerivativeEvaluation (const Point<dim>   &evaluation_point);

    virtual void operator () (const hp::DoFHandler<dim> &dof_handler,
                              const Vector<double>  &solution) const;

    DeclException1 (ExcEvaluationPointNotFound,
                    Point<dim>,
                    << "The evaluation point " << arg1
                    << " was not found among the vertices of the present grid.");
  private:
    const Point<dim>  evaluation_point;
  };


  template <int dim>
  PointXDerivativeEvaluation<dim>::
  PointXDerivativeEvaluation (const Point<dim>   &evaluation_point)
    :
    evaluation_point (evaluation_point)
  {}


  template <int dim>
  void
  PointXDerivativeEvaluation<dim>::
  operator () (const hp::DoFHandler<dim> &dof_handler,
               const Vector<double>  &solution) const
  {
    double point_derivative = 0;

    QTrapez<dim>  vertex_quadrature;
    hp::FEValues<dim> fe_values (dof_handler.get_fe(),
                                 vertex_quadrature,
                                 update_gradients | update_q_points);
    std::vector<Tensor<1,dim> >
    solution_gradients (vertex_quadrature.size());

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    unsigned int evaluation_point_hits = 0;
    for (; cell!=endc; ++cell)
      for (unsigned int vertex=0;
           vertex<GeometryInfo<dim>::vertices_per_cell;
           ++vertex)
        if (cell->vertex(vertex) == evaluation_point)
          {
            fe_values.reinit (cell);
            fe_values.get_present_fe_values().get_function_grads (solution,
                                                                  solution_gradients);

            unsigned int q_point = 0;
            for (; q_point<solution_gradients.size(); ++q_point)
              if (fe_values.get_present_fe_values().quadrature_point(q_point) ==
                  evaluation_point)
                break;

            Assert (q_point < solution_gradients.size(),
                    ExcInternalError());
            point_derivative += solution_gradients[q_point][0];
            ++evaluation_point_hits;

            break;
          };

    AssertThrow (evaluation_point_hits > 0,
                 ExcEvaluationPointNotFound(evaluation_point));

    point_derivative /= evaluation_point_hits;
    deallog << "   Point x-derivative=" << point_derivative
            << std::endl;
  }




  template <int dim>
  class GridOutput : public EvaluationBase<dim>
  {
  public:
    GridOutput (const std::string &output_name_base);

    virtual void operator () (const hp::DoFHandler<dim> &dof_handler,
                              const Vector<double>  &solution) const;
  private:
    const std::string output_name_base;
  };


  template <int dim>
  GridOutput<dim>::
  GridOutput (const std::string &output_name_base)
    :
    output_name_base (output_name_base)
  {}


  template <int dim>
  void
  GridOutput<dim>::operator () (const hp::DoFHandler<dim> &dof_handler,
                                const Vector<double>  &/*solution*/) const
  {
    std::ostringstream filename;
    filename << output_name_base << "-"
             << this->refinement_cycle
             << ".eps"
             << std::ends;

    GridOut().write_eps (dof_handler.get_tria(),
                         deallog.get_file_stream());
  }
}



namespace LaplaceSolver
{
  template <int dim> class WeightedResidual;



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

    virtual void set_refinement_cycle (const unsigned int cycle);

    virtual void output_solution () const = 0;

  protected:
    const SmartPointer<Triangulation<dim> > triangulation;

    unsigned int refinement_cycle;
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
  void
  Base<dim>::set_refinement_cycle (const unsigned int cycle)
  {
    refinement_cycle = cycle;
  }



  template <int dim>
  class Solver : public virtual Base<dim>
  {
  public:
    Solver (Triangulation<dim>       &triangulation,
            const hp::FECollection<dim> &fe,
            const hp::QCollection<dim>    &quadrature,
            const hp::QCollection<dim-1>  &face_quadrature,
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
    const SmartPointer<const hp::QCollection<dim-1> >   face_quadrature;
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
                       const hp::QCollection<dim-1>  &face_quadrature,
                       const Function<dim>      &boundary_values)
    :
    Base<dim> (triangulation),
    fe (&fe),
    quadrature (&quadrature),
    face_quadrature (&face_quadrature),
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
    mhnc_thread = Threads::spawn (mhnc_p)(dof_handler, hanging_node_constraints);

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
    SolverControl           solver_control (5000, 1e-12);
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
                  const hp::QCollection<dim-1>  &face_quadrature,
                  const Function<dim>      &rhs_function,
                  const Function<dim>      &boundary_values);

    virtual
    void solve_problem ();

    virtual
    unsigned int n_dofs () const;

    virtual
    void postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

    virtual
    void output_solution () const;

  protected:
    const SmartPointer<const Function<dim> > rhs_function;
    virtual void assemble_rhs (Vector<double> &rhs) const;

    friend class WeightedResidual<dim>;
  };


  template <int dim>
  PrimalSolver<dim>::
  PrimalSolver (Triangulation<dim>       &triangulation,
                const hp::FECollection<dim> &fe,
                const hp::QCollection<dim>    &quadrature,
                const hp::QCollection<dim-1>  &face_quadrature,
                const Function<dim>      &rhs_function,
                const Function<dim>      &boundary_values)
    :
    Base<dim> (triangulation),
    Solver<dim> (triangulation, fe,
                 quadrature, face_quadrature,
                 boundary_values),
    rhs_function (&rhs_function)
  {}


  template <int dim>
  void
  PrimalSolver<dim>::solve_problem ()
  {
    Solver<dim>::solve_problem ();
  }



  template <int dim>
  unsigned int
  PrimalSolver<dim>::n_dofs() const
  {
    return Solver<dim>::n_dofs();
  }


  template <int dim>
  void
  PrimalSolver<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    Solver<dim>::postprocess(postprocessor);
  }


  template <int dim>
  void
  PrimalSolver<dim>::output_solution () const
  {
    DataOut<dim,hp::DoFHandler<dim> > data_out;
    data_out.attach_dof_handler (this->dof_handler);
    data_out.add_data_vector (this->solution, "solution");
    data_out.build_patches ();

    std::ostringstream filename;
    filename << "solution-"
             << this->refinement_cycle
             << ".gnuplot"
             << std::ends;
    data_out.write (deallog.get_file_stream(), DataOut<dim>::gnuplot);
  }



  template <int dim>
  void
  PrimalSolver<dim>::
  assemble_rhs (Vector<double> &rhs) const
  {
    hp::FEValues<dim> fe_values (*this->fe, *this->quadrature,
                                 update_values  | update_q_points  |
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
      }
  }



  template <int dim>
  class RefinementGlobal : public PrimalSolver<dim>
  {
  public:
    RefinementGlobal (Triangulation<dim>       &coarse_grid,
                      const hp::FECollection<dim> &fe,
                      const hp::QCollection<dim>    &quadrature,
                      const hp::QCollection<dim-1>  &face_quadrature,
                      const Function<dim>      &rhs_function,
                      const Function<dim>      &boundary_values);

    virtual void refine_grid ();
  };



  template <int dim>
  RefinementGlobal<dim>::
  RefinementGlobal (Triangulation<dim>       &coarse_grid,
                    const hp::FECollection<dim> &fe,
                    const hp::QCollection<dim>    &quadrature,
                    const hp::QCollection<dim-1>  &face_quadrature,
                    const Function<dim>      &rhs_function,
                    const Function<dim>      &boundary_values)
    :
    Base<dim> (coarse_grid),
    PrimalSolver<dim> (coarse_grid, fe, quadrature,
                       face_quadrature, rhs_function,
                       boundary_values)
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
                     const hp::QCollection<dim-1>  &face_quadrature,
                     const Function<dim>      &rhs_function,
                     const Function<dim>      &boundary_values);

    virtual void refine_grid ();
  };



  template <int dim>
  RefinementKelly<dim>::
  RefinementKelly (Triangulation<dim>       &coarse_grid,
                   const hp::FECollection<dim> &fe,
                   const hp::QCollection<dim>    &quadrature,
                   const hp::QCollection<dim-1>  &face_quadrature,
                   const Function<dim>      &rhs_function,
                   const Function<dim>      &boundary_values)
    :
    Base<dim> (coarse_grid),
    PrimalSolver<dim> (coarse_grid, fe, quadrature,
                       face_quadrature,
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




  template <int dim>
  class RefinementWeightedKelly : public PrimalSolver<dim>
  {
  public:
    RefinementWeightedKelly (Triangulation<dim>       &coarse_grid,
                             const hp::FECollection<dim> &fe,
                             const hp::QCollection<dim>    &quadrature,
                             const hp::QCollection<dim-1>  &face_quadrature,
                             const Function<dim>      &rhs_function,
                             const Function<dim>      &boundary_values,
                             const Function<dim>      &weighting_function);

    virtual void refine_grid ();

  private:
    const SmartPointer<const Function<dim> > weighting_function;
  };



  template <int dim>
  RefinementWeightedKelly<dim>::
  RefinementWeightedKelly (Triangulation<dim>       &coarse_grid,
                           const hp::FECollection<dim> &fe,
                           const hp::QCollection<dim>    &quadrature,
                           const hp::QCollection<dim-1>  &face_quadrature,
                           const Function<dim>      &rhs_function,
                           const Function<dim>      &boundary_values,
                           const Function<dim>      &weighting_function)
    :
    Base<dim> (coarse_grid),
    PrimalSolver<dim> (coarse_grid, fe, quadrature,
                       face_quadrature,
                       rhs_function, boundary_values),
    weighting_function (&weighting_function)
  {}



  template <int dim>
  void
  RefinementWeightedKelly<dim>::refine_grid ()
  {
    Vector<float> estimated_error (this->triangulation->n_active_cells());
    KellyErrorEstimator<dim>::estimate (this->dof_handler,
                                        (*this->face_quadrature)[0],
                                        typename FunctionMap<dim>::type(),
                                        this->solution,
                                        estimated_error);

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = this->dof_handler.begin_active(),
    endc = this->dof_handler.end();
    for (unsigned int cell_index=0; cell!=endc; ++cell, ++cell_index)
      estimated_error(cell_index)
      *= weighting_function->value (cell->center());

    GridRefinement::refine_and_coarsen_fixed_number (*this->triangulation,
                                                     estimated_error,
                                                     0.3, 0.03);
    this->triangulation->execute_coarsening_and_refinement ();
  }

}


namespace Data
{

  template <int dim>
  struct SetUpBase : public Subscriptor
  {
    virtual
    const Function<dim>   &get_boundary_values () const = 0;

    virtual
    const Function<dim>   &get_right_hand_side () const = 0;

    virtual
    void create_coarse_grid (Triangulation<dim> &coarse_grid) const = 0;
  };


  template <class Traits, int dim>
  struct SetUp : public SetUpBase<dim>
  {
    SetUp () {}

    virtual
    const Function<dim>   &get_boundary_values () const;

    virtual
    const Function<dim>   &get_right_hand_side () const;


    virtual
    void create_coarse_grid (Triangulation<dim> &coarse_grid) const;

  private:
    static const typename Traits::BoundaryValues boundary_values;
    static const typename Traits::RightHandSide  right_hand_side;
  };

  template <class Traits, int dim>
  const typename Traits::BoundaryValues  SetUp<Traits,dim>::boundary_values;
  template <class Traits, int dim>
  const typename Traits::RightHandSide   SetUp<Traits,dim>::right_hand_side;

  template <class Traits, int dim>
  const Function<dim> &
  SetUp<Traits,dim>::get_boundary_values () const
  {
    return boundary_values;
  }


  template <class Traits, int dim>
  const Function<dim> &
  SetUp<Traits,dim>::get_right_hand_side () const
  {
    return right_hand_side;
  }


  template <class Traits, int dim>
  void
  SetUp<Traits,dim>::
  create_coarse_grid (Triangulation<dim> &coarse_grid) const
  {
    Traits::create_coarse_grid (coarse_grid);
  }



  template <int dim>
  struct CurvedRidges
  {
    class BoundaryValues : public Function<dim>
    {
    public:
      BoundaryValues () : Function<dim> () {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component) const;
    };


    class RightHandSide : public Function<dim>
    {
    public:
      RightHandSide () : Function<dim> () {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component) const;
    };

    static
    void
    create_coarse_grid (Triangulation<dim> &coarse_grid);
  };


  template <int dim>
  double
  CurvedRidges<dim>::BoundaryValues::
  value (const Point<dim>   &p,
         const unsigned int  /*component*/) const
  {
    double q = p(0);
    for (unsigned int i=1; i<dim; ++i)
      q += std::sin(10*p(i)+5*p(0)*p(0));
    const double exponential = std::exp(q);
    return exponential;
  }



  template <int dim>
  double
  CurvedRidges<dim>::RightHandSide::value (const Point<dim>   &p,
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
  CurvedRidges<dim>::
  create_coarse_grid (Triangulation<dim> &coarse_grid)
  {
    GridGenerator::hyper_cube (coarse_grid, -1, 1);
    coarse_grid.refine_global (2);
  }



  template <int dim>
  struct Exercise_2_3
  {
    typedef ZeroFunction<dim> BoundaryValues;

    class RightHandSide : public ConstantFunction<dim>
    {
    public:
      RightHandSide () : ConstantFunction<dim> (1.) {}
    };

    static
    void
    create_coarse_grid (Triangulation<dim> &coarse_grid);
  };


  template <>
  void
  Exercise_2_3<2>::
  create_coarse_grid (Triangulation<2> &coarse_grid)
  {
    const unsigned int dim = 2;

    static const Point<2> vertices_1[]
      = {  Point<2> (-1.,   -1.),
           Point<2> (-1./2, -1.),
           Point<2> (0.,    -1.),
           Point<2> (+1./2, -1.),
           Point<2> (+1,    -1.),

           Point<2> (-1.,   -1./2.),
           Point<2> (-1./2, -1./2.),
           Point<2> (0.,    -1./2.),
           Point<2> (+1./2, -1./2.),
           Point<2> (+1,    -1./2.),

           Point<2> (-1.,   0.),
           Point<2> (-1./2, 0.),
           Point<2> (+1./2, 0.),
           Point<2> (+1,    0.),

           Point<2> (-1.,   1./2.),
           Point<2> (-1./2, 1./2.),
           Point<2> (0.,    1./2.),
           Point<2> (+1./2, 1./2.),
           Point<2> (+1,    1./2.),

           Point<2> (-1.,   1.),
           Point<2> (-1./2, 1.),
           Point<2> (0.,    1.),
           Point<2> (+1./2, 1.),
           Point<2> (+1,    1.)
        };
    const unsigned int
    n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);

    const std::vector<Point<dim> > vertices (&vertices_1[0],
                                             &vertices_1[n_vertices]);

    static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
    = {{0, 1, 5, 6},
      {1, 2, 6, 7},
      {2, 3, 7, 8},
      {3, 4, 8, 9},
      {5, 6, 10, 11},
      {8, 9, 12, 13},
      {10, 11, 14, 15},
      {12, 13, 17, 18},
      {14, 15, 19, 20},
      {15, 16, 20, 21},
      {16, 17, 21, 22},
      {17, 18, 22, 23}
    };
    const unsigned int
    n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

    std::vector<CellData<dim> > cells (n_cells, CellData<dim>());
    for (unsigned int i=0; i<n_cells; ++i)
      {
        for (unsigned int j=0;
             j<GeometryInfo<dim>::vertices_per_cell;
             ++j)
          cells[i].vertices[j] = cell_vertices[i][j];
        cells[i].material_id = 0;
      };

    coarse_grid.create_triangulation (vertices,
                                      cells,
                                      SubCellData());

    coarse_grid.refine_global (1);
  }
}




namespace DualFunctional
{

  template <int dim>
  class DualFunctionalBase : public Subscriptor
  {
  public:
    virtual
    void
    assemble_rhs (const hp::DoFHandler<dim> &dof_handler,
                  Vector<double>        &rhs) const = 0;
  };



  template <int dim>
  class PointValueEvaluation : public DualFunctionalBase<dim>
  {
  public:
    PointValueEvaluation (const Point<dim> &evaluation_point);

    virtual
    void
    assemble_rhs (const hp::DoFHandler<dim> &dof_handler,
                  Vector<double>        &rhs) const;

    DeclException1 (ExcEvaluationPointNotFound,
                    Point<dim>,
                    << "The evaluation point " << arg1
                    << " was not found among the vertices of the present grid.");

  protected:
    const Point<dim> evaluation_point;
  };


  template <int dim>
  PointValueEvaluation<dim>::
  PointValueEvaluation (const Point<dim> &evaluation_point)
    :
    evaluation_point (evaluation_point)
  {}


  template <int dim>
  void
  PointValueEvaluation<dim>::
  assemble_rhs (const hp::DoFHandler<dim> &dof_handler,
                Vector<double>        &rhs) const
  {
    rhs.reinit (dof_handler.n_dofs());

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int vertex=0;
           vertex<GeometryInfo<dim>::vertices_per_cell;
           ++vertex)
        if (cell->vertex(vertex).distance(evaluation_point)
            < cell->diameter()*1e-8)
          {
            rhs(cell->vertex_dof_index(vertex,0)) = 1;
            return;
          };

    AssertThrow (false, ExcEvaluationPointNotFound(evaluation_point));
  }



  template <int dim>
  class PointXDerivativeEvaluation : public DualFunctionalBase<dim>
  {
  public:
    PointXDerivativeEvaluation (const Point<dim> &evaluation_point);

    virtual
    void
    assemble_rhs (const hp::DoFHandler<dim> &dof_handler,
                  Vector<double>        &rhs) const;

    DeclException1 (ExcEvaluationPointNotFound,
                    Point<dim>,
                    << "The evaluation point " << arg1
                    << " was not found among the vertices of the present grid.");

  protected:
    const Point<dim> evaluation_point;
  };


  template <int dim>
  PointXDerivativeEvaluation<dim>::
  PointXDerivativeEvaluation (const Point<dim> &evaluation_point)
    :
    evaluation_point (evaluation_point)
  {}


  template <int dim>
  void
  PointXDerivativeEvaluation<dim>::
  assemble_rhs (const hp::DoFHandler<dim> &dof_handler,
                Vector<double>        &rhs) const
  {
    rhs.reinit (dof_handler.n_dofs());

    QGauss<dim> quadrature(4);
    hp::FEValues<dim>  fe_values (dof_handler.get_fe(), quadrature,
                                  update_gradients |
                                  update_q_points  |
                                  update_JxW_values);
    const unsigned int n_q_points = fe_values.get_present_fe_values().n_quadrature_points;
    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    Vector<double> cell_rhs (dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    double total_volume = 0;

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->center().distance(evaluation_point) <=
          cell->diameter())
        {
          fe_values.reinit (cell);
          cell_rhs = 0;

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                cell_rhs(i) += fe_values.get_present_fe_values().shape_grad(i,q)[0] *
                               fe_values.get_present_fe_values().JxW (q);
              total_volume += fe_values.get_present_fe_values().JxW (q);
            };

          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            rhs(local_dof_indices[i]) += cell_rhs(i);
        };

    AssertThrow (total_volume > 0,
                 ExcEvaluationPointNotFound(evaluation_point));

    rhs /= total_volume;
  }


}


namespace LaplaceSolver
{


  template <int dim>
  class DualSolver : public Solver<dim>
  {
  public:
    DualSolver (Triangulation<dim>       &triangulation,
                const hp::FECollection<dim> &fe,
                const hp::QCollection<dim>    &quadrature,
                const hp::QCollection<dim-1>  &face_quadrature,
                const DualFunctional::DualFunctionalBase<dim> &dual_functional);

    virtual
    void
    solve_problem ();

    virtual
    unsigned int
    n_dofs () const;

    virtual
    void
    postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

  protected:
    const SmartPointer<const DualFunctional::DualFunctionalBase<dim> > dual_functional;
    virtual void assemble_rhs (Vector<double> &rhs) const;

    static const ZeroFunction<dim> boundary_values;

    friend class WeightedResidual<dim>;
  };

  template <int dim>
  const ZeroFunction<dim> DualSolver<dim>::boundary_values;

  template <int dim>
  DualSolver<dim>::
  DualSolver (Triangulation<dim>       &triangulation,
              const hp::FECollection<dim> &fe,
              const hp::QCollection<dim>    &quadrature,
              const hp::QCollection<dim-1>  &face_quadrature,
              const DualFunctional::DualFunctionalBase<dim> &dual_functional)
    :
    Base<dim> (triangulation),
    Solver<dim> (triangulation, fe,
                 quadrature, face_quadrature,
                 boundary_values),
    dual_functional (&dual_functional)
  {}


  template <int dim>
  void
  DualSolver<dim>::solve_problem ()
  {
    Solver<dim>::solve_problem ();
  }



  template <int dim>
  unsigned int
  DualSolver<dim>::n_dofs() const
  {
    return Solver<dim>::n_dofs();
  }


  template <int dim>
  void
  DualSolver<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    Solver<dim>::postprocess(postprocessor);
  }



  template <int dim>
  void
  DualSolver<dim>::
  assemble_rhs (Vector<double> &rhs) const
  {
    dual_functional->assemble_rhs (this->dof_handler, rhs);
  }



  template <int dim>
  class WeightedResidual : public PrimalSolver<dim>,
    public DualSolver<dim>
  {
  public:
    WeightedResidual (Triangulation<dim>       &coarse_grid,
                      const hp::FECollection<dim> &primal_fe,
                      const hp::FECollection<dim> &dual_fe,
                      const hp::QCollection<dim>    &quadrature,
                      const hp::QCollection<dim-1>  &face_quadrature,
                      const Function<dim>      &rhs_function,
                      const Function<dim>      &boundary_values,
                      const DualFunctional::DualFunctionalBase<dim> &dual_functional);

    virtual
    void
    solve_problem ();

    virtual
    void
    postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const;

    virtual
    unsigned int
    n_dofs () const;

    virtual void refine_grid ();

    virtual
    void
    output_solution () const;

  private:
    void solve_primal_problem ();
    void solve_dual_problem ();

    typedef
    typename hp::DoFHandler<dim>::active_cell_iterator
    active_cell_iterator;

    typedef
    typename std::map<typename hp::DoFHandler<dim>::face_iterator,double>
    FaceIntegrals;

    struct CellData
    {
      hp::FEValues<dim>    fe_values;
      const SmartPointer<const Function<dim> > right_hand_side;

      std::vector<double> cell_residual;
      std::vector<double> rhs_values;
      std::vector<double> dual_weights;
      typename std::vector<Tensor<2,dim> > cell_grad_grads;
      CellData (const hp::FECollection<dim> &fe,
                const hp::QCollection<dim>    &quadrature,
                const Function<dim>      &right_hand_side);
    };

    struct FaceData
    {
      hp::FEFaceValues<dim>    fe_face_values_cell;
      hp::FEFaceValues<dim>    fe_face_values_neighbor;
      hp::FESubfaceValues<dim> fe_subface_values_cell;

      std::vector<double> jump_residual;
      std::vector<double> dual_weights;
      typename std::vector<Tensor<1,dim> > cell_grads;
      typename std::vector<Tensor<1,dim> > neighbor_grads;
      FaceData (const hp::FECollection<dim> &fe,
                const hp::QCollection<dim-1>  &face_quadrature);
    };



    void estimate_error (Vector<float> &error_indicators) const;

    void estimate_some (const Vector<double> &primal_solution,
                        const Vector<double> &dual_weights,
                        const unsigned int    n_threads,
                        const unsigned int    this_thread,
                        Vector<float>        &error_indicators,
                        FaceIntegrals        &face_integrals) const;

    void
    integrate_over_cell (const active_cell_iterator &cell,
                         const unsigned int          cell_index,
                         const Vector<double>       &primal_solution,
                         const Vector<double>       &dual_weights,
                         CellData                   &cell_data,
                         Vector<float>              &error_indicators) const;

    void
    integrate_over_regular_face (const active_cell_iterator &cell,
                                 const unsigned int          face_no,
                                 const Vector<double>       &primal_solution,
                                 const Vector<double>       &dual_weights,
                                 FaceData                   &face_data,
                                 FaceIntegrals              &face_integrals) const;
    void
    integrate_over_irregular_face (const active_cell_iterator &cell,
                                   const unsigned int          face_no,
                                   const Vector<double>       &primal_solution,
                                   const Vector<double>       &dual_weights,
                                   FaceData                   &face_data,
                                   FaceIntegrals              &face_integrals) const;
  };



  template <int dim>
  WeightedResidual<dim>::CellData::
  CellData (const hp::FECollection<dim> &fe,
            const hp::QCollection<dim>    &quadrature,
            const Function<dim>      &right_hand_side)
    :
    fe_values (fe, quadrature,
               update_values             |
               update_second_derivatives |
               update_q_points           |
               update_JxW_values),
    right_hand_side (&right_hand_side)
  {
    const unsigned int n_q_points
      = quadrature[0].size();

    cell_residual.resize(n_q_points);
    rhs_values.resize(n_q_points);
    dual_weights.resize(n_q_points);
    cell_grad_grads.resize(n_q_points);
  }



  template <int dim>
  WeightedResidual<dim>::FaceData::
  FaceData (const hp::FECollection<dim> &fe,
            const hp::QCollection<dim-1>  &face_quadrature)
    :
    fe_face_values_cell (fe, face_quadrature,
                         update_values        |
                         update_gradients     |
                         update_JxW_values    |
                         update_normal_vectors),
    fe_face_values_neighbor (fe, face_quadrature,
                             update_values     |
                             update_gradients  |
                             update_JxW_values |
                             update_normal_vectors),
    fe_subface_values_cell (fe, face_quadrature,
                            update_gradients)
  {
    const unsigned int n_face_q_points
      = face_quadrature[0].size();

    jump_residual.resize(n_face_q_points);
    dual_weights.resize(n_face_q_points);
    cell_grads.resize(n_face_q_points);
    neighbor_grads.resize(n_face_q_points);
  }




  template <int dim>
  WeightedResidual<dim>::
  WeightedResidual (Triangulation<dim>       &coarse_grid,
                    const hp::FECollection<dim> &primal_fe,
                    const hp::FECollection<dim> &dual_fe,
                    const hp::QCollection<dim>    &quadrature,
                    const hp::QCollection<dim-1>  &face_quadrature,
                    const Function<dim>      &rhs_function,
                    const Function<dim>      &bv,
                    const DualFunctional::DualFunctionalBase<dim> &dual_functional)
    :
    Base<dim> (coarse_grid),
    PrimalSolver<dim> (coarse_grid, primal_fe,
                       quadrature, face_quadrature,
                       rhs_function, bv),
    DualSolver<dim> (coarse_grid, dual_fe,
                     quadrature, face_quadrature,
                     dual_functional)
  {}


  template <int dim>
  void
  WeightedResidual<dim>::solve_problem ()
  {
    Threads::ThreadGroup<> threads;
    threads += Threads::spawn (*this, &WeightedResidual<dim>::solve_primal_problem)();
    threads += Threads::spawn (*this, &WeightedResidual<dim>::solve_dual_problem)();
    threads.join_all ();
  }


  template <int dim>
  void
  WeightedResidual<dim>::solve_primal_problem ()
  {
    PrimalSolver<dim>::solve_problem ();
  }

  template <int dim>
  void
  WeightedResidual<dim>::solve_dual_problem ()
  {
    DualSolver<dim>::solve_problem ();
  }


  template <int dim>
  void
  WeightedResidual<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    PrimalSolver<dim>::postprocess (postprocessor);
  }


  template <int dim>
  unsigned int
  WeightedResidual<dim>::n_dofs () const
  {
    return PrimalSolver<dim>::n_dofs();
  }



  template <int dim>
  void
  WeightedResidual<dim>::refine_grid ()
  {
    Vector<float> error_indicators (this->triangulation->n_active_cells());
    estimate_error (error_indicators);

    for (Vector<float>::iterator i=error_indicators.begin();
         i != error_indicators.end(); ++i)
      *i = std::fabs (*i);

    GridRefinement::refine_and_coarsen_fixed_fraction (*this->triangulation,
                                                       error_indicators,
                                                       0.8, 0.02);
    this->triangulation->execute_coarsening_and_refinement ();
  }


  template <int dim>
  void
  WeightedResidual<dim>::output_solution () const
  {
    const PrimalSolver<dim> &primal_solver = *this;
    const DualSolver<dim>   &dual_solver   = *this;

    ConstraintMatrix primal_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints (primal_solver.dof_handler,
                                             primal_hanging_node_constraints);
    primal_hanging_node_constraints.close();
    Vector<double> dual_solution (primal_solver.dof_handler.n_dofs());
    FETools::interpolate (dual_solver.dof_handler,
                          dual_solver.solution,
                          primal_solver.dof_handler,
                          primal_hanging_node_constraints,
                          dual_solution);

    DataOut<dim,hp::DoFHandler<dim> > data_out;
    data_out.attach_dof_handler (primal_solver.dof_handler);

    data_out.add_data_vector (primal_solver.solution,
                              "primal_solution");
    data_out.add_data_vector (dual_solution,
                              "dual_solution");

    data_out.build_patches ();

    std::ostringstream filename;
    filename << "solution-"
             << this->refinement_cycle
             << ".gnuplot"
             << std::ends;
    data_out.write (deallog.get_file_stream(), DataOut<dim>::gnuplot);
  }




  template <int dim>
  void
  WeightedResidual<dim>::
  estimate_error (Vector<float> &error_indicators) const
  {
    const PrimalSolver<dim> &primal_solver = *this;
    const DualSolver<dim>   &dual_solver   = *this;

    ConstraintMatrix dual_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints (dual_solver.dof_handler,
                                             dual_hanging_node_constraints);
    dual_hanging_node_constraints.close();
    Vector<double> primal_solution (dual_solver.dof_handler.n_dofs());
    FETools::interpolate (primal_solver.dof_handler,
                          primal_solver.solution,
                          dual_solver.dof_handler,
                          dual_hanging_node_constraints,
                          primal_solution);

    ConstraintMatrix primal_hanging_node_constraints;
    DoFTools::make_hanging_node_constraints (primal_solver.dof_handler,
                                             primal_hanging_node_constraints);
    primal_hanging_node_constraints.close();
    Vector<double> dual_weights (dual_solver.dof_handler.n_dofs());
    FETools::interpolation_difference (dual_solver.dof_handler,
                                       dual_hanging_node_constraints,
                                       dual_solver.solution,
                                       primal_solver.dof_handler,
                                       primal_hanging_node_constraints,
                                       dual_weights);


    FaceIntegrals face_integrals;
    for (active_cell_iterator cell=dual_solver.dof_handler.begin_active();
         cell!=dual_solver.dof_handler.end();
         ++cell)
      for (unsigned int face_no=0;
           face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        face_integrals[cell->face(face_no)] = -1e20;

    error_indicators.reinit (dual_solver.dof_handler
                             .get_tria().n_active_cells());

    const unsigned int n_threads = multithread_info.n_default_threads;
    Threads::ThreadGroup<> threads;
    for (unsigned int i=0; i<n_threads; ++i)
      threads += Threads::spawn (*this, &WeightedResidual<dim>::estimate_some)
                 (primal_solution,
                  dual_weights,
                  n_threads, i,
                  error_indicators,
                  face_integrals);
    threads.join_all();

    unsigned int present_cell=0;
    for (active_cell_iterator cell=dual_solver.dof_handler.begin_active();
         cell!=dual_solver.dof_handler.end();
         ++cell, ++present_cell)
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        {
          Assert(face_integrals.find(cell->face(face_no)) !=
                 face_integrals.end(),
                 ExcInternalError());
          error_indicators(present_cell)
          -= 0.5*face_integrals[cell->face(face_no)];
        };
    deallog << "   Estimated error="
            << std::accumulate (error_indicators.begin(),
                                error_indicators.end(), 0.)
            << std::endl;
  }



  template <int dim>
  void
  WeightedResidual<dim>::
  estimate_some (const Vector<double> &primal_solution,
                 const Vector<double> &dual_weights,
                 const unsigned int    n_threads,
                 const unsigned int    this_thread,
                 Vector<float>        &error_indicators,
                 FaceIntegrals        &face_integrals) const
  {
    const PrimalSolver<dim> &primal_solver = *this;
    const DualSolver<dim>   &dual_solver   = *this;

    CellData cell_data (*dual_solver.fe,
                        *dual_solver.quadrature,
                        *primal_solver.rhs_function);
    FaceData face_data (*dual_solver.fe,
                        *dual_solver.face_quadrature);

    active_cell_iterator cell=dual_solver.dof_handler.begin_active();
    for (unsigned int t=0;
         (t<this_thread) && (cell!=dual_solver.dof_handler.end());
         ++t, ++cell);

    if (cell == dual_solver.dof_handler.end())
      return;

    for (unsigned int cell_index=this_thread; true; )
      {
        integrate_over_cell (cell, cell_index,
                             primal_solution,
                             dual_weights,
                             cell_data,
                             error_indicators);

        for (unsigned int face_no=0;
             face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          {
            if (cell->face(face_no)->at_boundary())
              {
                face_integrals[cell->face(face_no)] = 0;
                continue;
              };

            if ((cell->neighbor(face_no)->has_children() == false) &&
                (cell->neighbor(face_no)->level() == cell->level()) &&
                (cell->neighbor(face_no)->index() < cell->index()))
              continue;

            if (cell->at_boundary(face_no) == false)
              if (cell->neighbor(face_no)->level() < cell->level())
                continue;


            if (cell->face(face_no)->has_children() == false)
              integrate_over_regular_face (cell, face_no,
                                           primal_solution,
                                           dual_weights,
                                           face_data,
                                           face_integrals);
            else
              integrate_over_irregular_face (cell, face_no,
                                             primal_solution,
                                             dual_weights,
                                             face_data,
                                             face_integrals);
          };

        for (unsigned int t=0;
             ((t<n_threads) && (cell!=dual_solver.dof_handler.end()));
             ++t, ++cell, ++cell_index);
        if (cell == dual_solver.dof_handler.end())
          break;
      };
  }



  template <int dim>
  void WeightedResidual<dim>::
  integrate_over_cell (const active_cell_iterator &cell,
                       const unsigned int          cell_index,
                       const Vector<double>       &primal_solution,
                       const Vector<double>       &dual_weights,
                       CellData                   &cell_data,
                       Vector<float>              &error_indicators) const
  {
    cell_data.fe_values.reinit (cell);
    cell_data.right_hand_side
    ->value_list (cell_data.fe_values.get_present_fe_values().get_quadrature_points(),
                  cell_data.rhs_values);
    cell_data.fe_values.get_present_fe_values().get_function_2nd_derivatives (primal_solution,
        cell_data.cell_grad_grads);

    cell_data.fe_values.get_present_fe_values().get_function_values (dual_weights,
        cell_data.dual_weights);

    double sum = 0;
    for (unsigned int p=0; p<cell_data.fe_values.get_present_fe_values().n_quadrature_points; ++p)
      sum += ((cell_data.rhs_values[p]+trace(cell_data.cell_grad_grads[p])) *
              cell_data.dual_weights[p] *
              cell_data.fe_values.get_present_fe_values().JxW (p));
    error_indicators(cell_index) += sum;
  }



  template <int dim>
  void WeightedResidual<dim>::
  integrate_over_regular_face (const active_cell_iterator &cell,
                               const unsigned int          face_no,
                               const Vector<double>       &primal_solution,
                               const Vector<double>       &dual_weights,
                               FaceData                   &face_data,
                               FaceIntegrals              &face_integrals) const
  {
    const unsigned int
    n_q_points = face_data.fe_face_values_cell.get_present_fe_values().n_quadrature_points;

    face_data.fe_face_values_cell.reinit (cell, face_no);
    face_data.fe_face_values_cell.get_present_fe_values().get_function_grads (primal_solution,
        face_data.cell_grads);

    Assert (cell->neighbor(face_no).state() == IteratorState::valid,
            ExcInternalError());
    const unsigned int
    neighbor_neighbor = cell->neighbor_of_neighbor (face_no);
    const active_cell_iterator neighbor = cell->neighbor(face_no);
    face_data.fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor);
    face_data.fe_face_values_neighbor.get_present_fe_values().get_function_grads (primal_solution,
        face_data.neighbor_grads);

    for (unsigned int p=0; p<n_q_points; ++p)
      face_data.jump_residual[p]
        = ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) *
           face_data.fe_face_values_cell.get_present_fe_values().normal_vector(p));

    face_data.fe_face_values_cell.get_present_fe_values().get_function_values (dual_weights,
        face_data.dual_weights);

    double face_integral = 0;
    for (unsigned int p=0; p<n_q_points; ++p)
      face_integral += (face_data.jump_residual[p] *
                        face_data.dual_weights[p]  *
                        face_data.fe_face_values_cell.get_present_fe_values().JxW(p));

    Assert (face_integrals.find (cell->face(face_no)) != face_integrals.end(),
            ExcInternalError());
    Assert (face_integrals[cell->face(face_no)] == -1e20,
            ExcInternalError());

    face_integrals[cell->face(face_no)] = face_integral;
  }



  template <int dim>
  void WeightedResidual<dim>::
  integrate_over_irregular_face (const active_cell_iterator &cell,
                                 const unsigned int          face_no,
                                 const Vector<double>       &primal_solution,
                                 const Vector<double>       &dual_weights,
                                 FaceData                   &face_data,
                                 FaceIntegrals              &face_integrals) const
  {
    const unsigned int
    n_q_points = face_data.fe_face_values_cell.get_present_fe_values().n_quadrature_points;

    const typename hp::DoFHandler<dim>::face_iterator
    face = cell->face(face_no);
    const typename hp::DoFHandler<dim>::cell_iterator
    neighbor = cell->neighbor(face_no);
    Assert (neighbor.state() == IteratorState::valid,
            ExcInternalError());
    Assert (neighbor->has_children(),
            ExcInternalError());

    const unsigned int
    neighbor_neighbor = cell->neighbor_of_neighbor (face_no);

    for (unsigned int subface_no=0;
         subface_no<face->n_children(); ++subface_no)
      {
        const active_cell_iterator neighbor_child
          = cell->neighbor_child_on_subface (face_no, subface_no);
        Assert (neighbor_child->face(neighbor_neighbor) ==
                cell->face(face_no)->child(subface_no),
                ExcInternalError());

        face_data.fe_subface_values_cell.reinit (cell, face_no, subface_no);
        face_data.fe_subface_values_cell.get_present_fe_values().get_function_grads (primal_solution,
            face_data.cell_grads);
        face_data.fe_face_values_neighbor.reinit (neighbor_child,
                                                  neighbor_neighbor);
        face_data.fe_face_values_neighbor.get_present_fe_values().get_function_grads (primal_solution,
            face_data.neighbor_grads);

        for (unsigned int p=0; p<n_q_points; ++p)
          face_data.jump_residual[p]
            = ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) *
               face_data.fe_face_values_neighbor.get_present_fe_values().normal_vector(p));

        face_data.fe_face_values_neighbor.get_present_fe_values().get_function_values (dual_weights,
            face_data.dual_weights);

        double face_integral = 0;
        for (unsigned int p=0; p<n_q_points; ++p)
          face_integral += (face_data.jump_residual[p] *
                            face_data.dual_weights[p] *
                            face_data.fe_face_values_neighbor.get_present_fe_values().JxW(p));
        face_integrals[neighbor_child->face(neighbor_neighbor)]
          = face_integral;
      };

    double sum = 0;
    for (unsigned int subface_no=0;
         subface_no<face->n_children(); ++subface_no)
      {
        Assert (face_integrals.find(face->child(subface_no)) !=
                face_integrals.end(),
                ExcInternalError());
        Assert (face_integrals[face->child(subface_no)] != -1e20,
                ExcInternalError());

        sum += face_integrals[face->child(subface_no)];
      };
    face_integrals[face] = sum;
  }

}



template <int dim>
struct Framework
{
public:
  typedef Evaluation::EvaluationBase<dim> Evaluator;
  typedef std::list<Evaluator *>           EvaluatorList;


  struct ProblemDescription
  {
    unsigned int primal_fe_degree;
    unsigned int dual_fe_degree;

    SmartPointer<const Data::SetUpBase<dim> > data;

    enum RefinementCriterion
    {
      dual_weighted_error_estimator,
      global_refinement,
      kelly_indicator,
      weighted_kelly_indicator
    };

    RefinementCriterion refinement_criterion;

    SmartPointer<const DualFunctional::DualFunctionalBase<dim> > dual_functional;

    EvaluatorList evaluator_list;

    SmartPointer<const Function<dim> > kelly_weight;

    unsigned int max_degrees_of_freedom;

    ProblemDescription ();
  };

  static void run (const ProblemDescription &descriptor);
};


template <int dim>
Framework<dim>::ProblemDescription::ProblemDescription ()
  :
  primal_fe_degree (1),
  dual_fe_degree (2),
  refinement_criterion (dual_weighted_error_estimator),
  max_degrees_of_freedom (5000)
{}



template <int dim>
void Framework<dim>::run (const ProblemDescription &descriptor)
{
  Triangulation<dim>
  triangulation (Triangulation<dim>::smoothing_on_refinement);
  descriptor.data->create_coarse_grid (triangulation);

  const hp::FECollection<dim>     primal_fe(FE_Q<dim>(descriptor.primal_fe_degree));
  const hp::FECollection<dim>     dual_fe(FE_Q<dim>(descriptor.dual_fe_degree));
  const hp::QCollection<dim>   quadrature(QGauss<dim>(descriptor.dual_fe_degree+1));
  const hp::QCollection<dim-1> face_quadrature(QGauss<dim-1>(descriptor.dual_fe_degree+1));

  LaplaceSolver::Base<dim> *solver = 0;
  switch (descriptor.refinement_criterion)
    {
    case ProblemDescription::dual_weighted_error_estimator:
    {
      solver
        = new LaplaceSolver::WeightedResidual<dim> (triangulation,
                                                    primal_fe,
                                                    dual_fe,
                                                    quadrature,
                                                    face_quadrature,
                                                    descriptor.data->get_right_hand_side(),
                                                    descriptor.data->get_boundary_values(),
                                                    *descriptor.dual_functional);
      break;
    };

    case ProblemDescription::global_refinement:
    {
      solver
        = new LaplaceSolver::RefinementGlobal<dim> (triangulation,
                                                    primal_fe,
                                                    quadrature,
                                                    face_quadrature,
                                                    descriptor.data->get_right_hand_side(),
                                                    descriptor.data->get_boundary_values());
      break;
    };

    case ProblemDescription::kelly_indicator:
    {
      solver
        = new LaplaceSolver::RefinementKelly<dim> (triangulation,
                                                   primal_fe,
                                                   quadrature,
                                                   face_quadrature,
                                                   descriptor.data->get_right_hand_side(),
                                                   descriptor.data->get_boundary_values());
      break;
    };

    case ProblemDescription::weighted_kelly_indicator:
    {
      solver
        = new LaplaceSolver::RefinementWeightedKelly<dim> (triangulation,
                                                           primal_fe,
                                                           quadrature,
                                                           face_quadrature,
                                                           descriptor.data->get_right_hand_side(),
                                                           descriptor.data->get_boundary_values(),
                                                           *descriptor.kelly_weight);
      break;
    };

    default:
      AssertThrow (false, ExcInternalError());
    };

  for (unsigned int step=0; true; ++step)
    {
      deallog << "Refinement cycle: " << step
              << std::endl;

      solver->set_refinement_cycle (step);
      solver->solve_problem ();
      solver->output_solution ();

      deallog << "   Number of degrees of freedom="
              << solver->n_dofs() << std::endl;

      for (typename EvaluatorList::const_iterator
           e = descriptor.evaluator_list.begin();
           e != descriptor.evaluator_list.end(); ++e)
        {
          (*e)->set_refinement_cycle (step);
          solver->postprocess (**e);
        };


      if (solver->n_dofs() < descriptor.max_degrees_of_freedom)
        solver->refine_grid ();
      else
        break;
    };

  deallog << std::endl;
  delete solver;
  solver = 0;
}





int main ()
{
  deallog.depth_console (0);
  try
    {
      logfile.precision(2);

      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      const unsigned int dim = 2;
      Framework<dim>::ProblemDescription descriptor;

      descriptor.refinement_criterion
        = Framework<dim>::ProblemDescription::dual_weighted_error_estimator;

      descriptor.primal_fe_degree = 1;
      descriptor.dual_fe_degree   = 2;

      descriptor.data = new Data::SetUp<Data::Exercise_2_3<dim>,dim> ();

      const Point<dim> evaluation_point (0.75, 0.75);
      descriptor.dual_functional
        = new DualFunctional::PointValueEvaluation<dim> (evaluation_point);

      Evaluation::PointValueEvaluation<dim>
      postprocessor1 (evaluation_point);
      Evaluation::GridOutput<dim>
      postprocessor2 ("grid");

      descriptor.evaluator_list.push_back (&postprocessor1);
      descriptor.evaluator_list.push_back (&postprocessor2);

      descriptor.max_degrees_of_freedom = 20000;

      Framework<dim>::run (descriptor);
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
