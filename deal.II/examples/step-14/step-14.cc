/* $Id$ */
/* Author: Wolfgang Bangerth, ETH Zurich, 2002 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2002 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/table_handler.h>
#include <base/thread_management.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_refinement.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <numeric>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif



namespace Evaluation
{

  template <int dim>
  class EvaluationBase 
  {
    public:
      virtual ~EvaluationBase ();

      void set_refinement_cycle (const unsigned int refinement_cycle);
      
      virtual void operator () (const DoFHandler<dim> &dof_handler,
				const Vector<double>  &solution) const = 0;
    protected:
      unsigned int refinement_cycle;
  };


  template <int dim>
  EvaluationBase<dim>::~EvaluationBase ()
  {};
  

  
  template <int dim>
  void
  EvaluationBase<dim>::set_refinement_cycle (const unsigned int step)
  {
    refinement_cycle = step;
  };



  template <int dim>
  class PointValueEvaluation : public EvaluationBase<dim>
  {
    public:
      PointValueEvaluation (const Point<dim>   &evaluation_point,
			    TableHandler       &results_table);
      
      virtual void operator () (const DoFHandler<dim> &dof_handler,
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
  {};
  


  template <int dim>
  void
  PointValueEvaluation<dim>::
  operator () (const DoFHandler<dim> &dof_handler,
	       const Vector<double>  &solution) const 
  {
    double point_value = 1e20;

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    bool evaluation_point_found = false;
    for (; (cell!=endc) && !evaluation_point_found; ++cell)
      for (unsigned int vertex=0;
	   vertex<GeometryInfo<dim>::vertices_per_cell;
	   ++vertex)
	if (cell->vertex(vertex) == evaluation_point)
	  {
	    point_value = 1.-solution(cell->vertex_dof_index(vertex,0));

	    evaluation_point_found = true;
	    break;
	  };

    AssertThrow (evaluation_point_found,
		 ExcEvaluationPointNotFound(evaluation_point));

    results_table.add_value ("DoFs", dof_handler.n_dofs());
    results_table.add_value ("u(x_0)", point_value);
  };

};

  

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
  {};


  template <int dim>
  Base<dim>::~Base () 
  {};



  template <int dim>
  void
  Base<dim>::set_refinement_cycle (const unsigned int cycle)
  {
    refinement_cycle = cycle;
  };
  


  template <int dim>
  class Solver : public virtual Base<dim>
  {
    public:
      Solver (Triangulation<dim>       &triangulation,
	      const FiniteElement<dim> &fe,
	      const Quadrature<dim>    &quadrature,
	      const Quadrature<dim-1>  &face_quadrature,	      
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
      const SmartPointer<const FiniteElement<dim> >  fe;
      const SmartPointer<const Quadrature<dim> >     quadrature;
      const SmartPointer<const Quadrature<dim-1> >   face_quadrature;      
      DoFHandler<dim>                                dof_handler;
      Vector<double>                                 solution;
      const SmartPointer<const Function<dim> >       boundary_values;

      virtual void assemble_rhs (Vector<double> &rhs) const = 0;
    
    private:
      struct LinearSystem
      {
	  LinearSystem (const DoFHandler<dim> &dof_handler);

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
		       const typename DoFHandler<dim>::active_cell_iterator &begin_cell,
		       const typename DoFHandler<dim>::active_cell_iterator &end_cell,
		       Threads::ThreadMutex                                 &mutex) const      ;
  };



  template <int dim>
  Solver<dim>::Solver (Triangulation<dim>       &triangulation,
		       const FiniteElement<dim> &fe,
		       const Quadrature<dim>    &quadrature,
		       const Quadrature<dim-1>  &face_quadrature,
		       const Function<dim>      &boundary_values)
		  :
		  Base<dim> (triangulation),
		  fe (&fe),
                  quadrature (&quadrature),
                  face_quadrature (&face_quadrature),    
		  dof_handler (triangulation),
		  boundary_values (&boundary_values)
  {};


  template <int dim>
  Solver<dim>::~Solver () 
  {
    dof_handler.clear ();
  };


  template <int dim>
  void
  Solver<dim>::solve_problem ()
  {
    dof_handler.distribute_dofs (*fe);
    solution.reinit (dof_handler.n_dofs());

    LinearSystem linear_system (dof_handler);
    assemble_linear_system (linear_system);
    linear_system.solve (solution);
  };


  template <int dim>
  void
  Solver<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    postprocessor (dof_handler, solution);
  };


  template <int dim>
  unsigned int
  Solver<dim>::n_dofs () const
  {
    return dof_handler.n_dofs();
  };
  

  template <int dim>
  void
  Solver<dim>::assemble_linear_system (LinearSystem &linear_system)
  {
    typedef
      typename DoFHandler<dim>::active_cell_iterator
      active_cell_iterator;

    const unsigned int n_threads = multithread_info.n_default_threads;
    std::vector<std::pair<active_cell_iterator,active_cell_iterator> >
      thread_ranges 
      = Threads::split_range<active_cell_iterator> (dof_handler.begin_active (),
						    dof_handler.end (),
						    n_threads);

    Threads::ThreadMutex mutex;
    Threads::ThreadManager thread_manager;
    for (unsigned int thread=0; thread<n_threads; ++thread)
      Threads::spawn (thread_manager,
		      Threads::encapsulate(&Solver<dim>::assemble_matrix)
		      .collect_args (this,
				     linear_system,
				     thread_ranges[thread].first,
				     thread_ranges[thread].second,
				     mutex));

    assemble_rhs (linear_system.rhs);
    linear_system.hanging_node_constraints.condense (linear_system.rhs);

    std::map<unsigned int,double> boundary_value_map;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      *boundary_values,
					      boundary_value_map);
    
    
    thread_manager.wait ();
    linear_system.hanging_node_constraints.condense (linear_system.matrix);

    MatrixTools::apply_boundary_values (boundary_value_map,
					linear_system.matrix,
					solution,
					linear_system.rhs);

  };


  template <int dim>
  void
  Solver<dim>::assemble_matrix (LinearSystem                                         &linear_system,
				const typename DoFHandler<dim>::active_cell_iterator &begin_cell,
				const typename DoFHandler<dim>::active_cell_iterator &end_cell,
				Threads::ThreadMutex                                 &mutex) const
  {
    FEValues<dim> fe_values (*fe, *quadrature, 
			     UpdateFlags(update_gradients |
					 update_JxW_values));

    const unsigned int   dofs_per_cell = fe->dofs_per_cell;
    const unsigned int   n_q_points    = quadrature->n_quadrature_points;

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    for (typename DoFHandler<dim>::active_cell_iterator cell=begin_cell;
	 cell!=end_cell; ++cell)
      {
	cell_matrix.clear ();

	fe_values.reinit (cell);
	const std::vector<std::vector<Tensor<1,dim> > >
	  & shape_grads  = fe_values.get_shape_grads();
	const std::vector<double>
	  & JxW_values   = fe_values.get_JxW_values();

	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (shape_grads[i][q_point] *
				   shape_grads[j][q_point] *
				   JxW_values[q_point]);


	cell->get_dof_indices (local_dof_indices);
	mutex.acquire ();
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    linear_system.matrix.add (local_dof_indices[i],
				      local_dof_indices[j],
				      cell_matrix(i,j));
	mutex.release ();
      };
  };


  template <int dim>
  Solver<dim>::LinearSystem::
  LinearSystem (const DoFHandler<dim> &dof_handler)
  {
    hanging_node_constraints.clear ();

    void (*mhnc_p) (const DoFHandler<dim> &,
		    ConstraintMatrix      &)
      = &DoFTools::make_hanging_node_constraints;
    
    Threads::ThreadManager thread_manager;
    Threads::spawn (thread_manager,
		    Threads::encapsulate (mhnc_p)
		    .collect_args (dof_handler,
				   hanging_node_constraints));

    sparsity_pattern.reinit (dof_handler.n_dofs(),
			     dof_handler.n_dofs(),
			     dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

    thread_manager.wait ();
    hanging_node_constraints.close ();
    hanging_node_constraints.condense (sparsity_pattern);

    sparsity_pattern.compress();
    matrix.reinit (sparsity_pattern);
    rhs.reinit (dof_handler.n_dofs());
  };



  template <int dim>
  void
  Solver<dim>::LinearSystem::solve (Vector<double> &solution) const
  {
    SolverControl           solver_control (1000, 1e-12);
    PrimitiveVectorMemory<> vector_memory;
    SolverCG<>              cg (solver_control, vector_memory);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(matrix, 1.2);

    cg.solve (matrix, solution, rhs, preconditioner);

    hanging_node_constraints.distribute (solution);
  };





  template <int dim>
  class PrimalSolver : public Solver<dim>
  {
    public:
      PrimalSolver (Triangulation<dim>       &triangulation,
		    const FiniteElement<dim> &fe,
		    const Quadrature<dim>    &quadrature,
		    const Quadrature<dim-1>  &face_quadrature,
		    const Function<dim>      &rhs_function,
		    const Function<dim>      &boundary_values);

				       // XXX
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
      const SmartPointer<const Function<dim> > rhs_function;
      virtual void assemble_rhs (Vector<double> &rhs) const;
  };


  template <int dim>
  PrimalSolver<dim>::
  PrimalSolver (Triangulation<dim>       &triangulation,
		const FiniteElement<dim> &fe,
		const Quadrature<dim>    &quadrature,
		const Quadrature<dim-1>  &face_quadrature,
		const Function<dim>      &rhs_function,
		const Function<dim>      &boundary_values)
		  :
		  Base<dim> (triangulation),
		  Solver<dim> (triangulation, fe,
			       quadrature, face_quadrature,
			       boundary_values),
                  rhs_function (&rhs_function)
  {};


  template <int dim>
  void
  PrimalSolver<dim>::solve_problem ()
  {
    Solver<dim>::solve_problem ();
  };



  template <int dim>
  unsigned int
  PrimalSolver<dim>::n_dofs() const
  {
    return Solver<dim>::n_dofs();
  };


  template <int dim>
  void
  PrimalSolver<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    return Solver<dim>::postprocess(postprocessor);
  };
  


  template <int dim>
  void
  PrimalSolver<dim>::
  assemble_rhs (Vector<double> &rhs) const 
  {
    FEValues<dim> fe_values (*fe, *quadrature, 
			     UpdateFlags(update_values    |
					 update_q_points  |
					 update_JxW_values));

    const unsigned int   dofs_per_cell = fe->dofs_per_cell;
    const unsigned int   n_q_points    = quadrature->n_quadrature_points;

    Vector<double>       cell_rhs (dofs_per_cell);
    std::vector<double>  rhs_values (n_q_points);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
	cell_rhs.clear ();

	fe_values.reinit (cell);
	const FullMatrix<double> 
	  & shape_values = fe_values.get_shape_values();
	const std::vector<double>
	  & JxW_values   = fe_values.get_JxW_values();
	const std::vector<Point<dim> >
	  & q_points     = fe_values.get_quadrature_points();

	rhs_function->value_list (q_points, rhs_values);
      
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    cell_rhs(i) += (shape_values (i,q_point) *
			    rhs_values[q_point] *
			    JxW_values[q_point]);

	cell->get_dof_indices (local_dof_indices);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  rhs(local_dof_indices[i]) += cell_rhs(i);
      };
  };



  template <int dim>
  class RefinementGlobal : public PrimalSolver<dim>
  {
    public:
      RefinementGlobal (Triangulation<dim>       &coarse_grid,
			const FiniteElement<dim> &fe,
			const Quadrature<dim>    &quadrature,
			const Function<dim>      &rhs_function,
			const Function<dim>      &boundary_values);

      virtual void refine_grid ();
  };



  template <int dim>
  RefinementGlobal<dim>::
  RefinementGlobal (Triangulation<dim>       &coarse_grid,
		    const FiniteElement<dim> &fe,
		    const Quadrature<dim>    &quadrature,
		    const Function<dim>      &rhs_function,
		    const Function<dim>      &boundary_values)
		  :
		  Base<dim> (coarse_grid),
                  PrimalSolver<dim> (coarse_grid, fe, quadrature,
				     rhs_function, boundary_values)
  {};



  template <int dim>
  void
  RefinementGlobal<dim>::refine_grid ()
  {
    triangulation->refine_global (1);
  };



  template <int dim>
  class RefinementKelly : public PrimalSolver<dim>
  {
    public:
      RefinementKelly (Triangulation<dim>       &coarse_grid,
		       const FiniteElement<dim> &fe,
		       const Quadrature<dim>    &quadrature,
		       const Quadrature<dim-1>  &face_quadrature,
		       const Function<dim>      &rhs_function,
		       const Function<dim>      &boundary_values);

      virtual void refine_grid ();
  };



  template <int dim>
  RefinementKelly<dim>::
  RefinementKelly (Triangulation<dim>       &coarse_grid,
		   const FiniteElement<dim> &fe,
		   const Quadrature<dim>    &quadrature,
		   const Quadrature<dim-1>  &face_quadrature,
		   const Function<dim>      &rhs_function,
		   const Function<dim>      &boundary_values)
		  :
		  Base<dim> (coarse_grid),
                  PrimalSolver<dim> (coarse_grid, fe, quadrature,
				     face_quadrature,
				     rhs_function, boundary_values)
  {};



  template <int dim>
  void
  RefinementKelly<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (triangulation->n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
					QGauss3<dim-1>(),
					typename FunctionMap<dim>::type(),
					solution,
					estimated_error_per_cell);
    GridRefinement::refine_and_coarsen_fixed_number (*triangulation,
						     estimated_error_per_cell,
						     0.3, 0.03);
    triangulation->execute_coarsening_and_refinement ();
  };

};


template <int dim>
class Solution : public Function<dim>
{
  public:
    Solution () : Function<dim> () {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const;
};


template <int dim>
double
Solution<dim>::value (const Point<dim>   &p,
		      const unsigned int  /*component*/) const
{
//    double q = p(0);
//    for (unsigned int i=1; i<dim; ++i)
//      q += sin(10*p(i)+5*p(0)*p(0));
//    const double exponential = exp(q);
//    return exponential;
  double s = 1;
  for (unsigned int i=0; i<dim; ++i)
    s *= sin(3.1415926536*p(i));
  return s;
};



template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>   &p,
			   const unsigned int  /*component*/) const
{
//    double q = p(0);
//    for (unsigned int i=1; i<dim; ++i)
//      q += sin(10*p(i)+5*p(0)*p(0));
//    const double u = exp(q);
//    double t1 = 1,
//  	 t2 = 0,
//  	 t3 = 0;
//    for (unsigned int i=1; i<dim; ++i)
//      {
//        t1 += cos(10*p(i)+5*p(0)*p(0)) * 10 * p(0);
//        t2 += 10*cos(10*p(i)+5*p(0)*p(0)) -
//  	    100*sin(10*p(i)+5*p(0)*p(0)) * p(0)*p(0);
//        t3 += 100*cos(10*p(i)+5*p(0)*p(0))*cos(10*p(i)+5*p(0)*p(0)) -
//  	    100*sin(10*p(i)+5*p(0)*p(0));
//      };
//    t1 = t1*t1;
  
//    return -u*(t1+t2+t3);
  double s = 1;
  for (unsigned int i=0; i<dim; ++i)
    s *= sin(3.1415926536*p(i));
  return dim*3.1415926536*3.1415926536*s;  
};




namespace DualFunctional
{
  template <int dim>
  class DualFunctionalBase : public Subscriptor
  {
    public:
      virtual
      void
      assemble_rhs (const DoFHandler<dim> &dof_handler,
		    Vector<double>        &rhs) const = 0;
  };


  template <int dim>
  class PointValueEvaluation : public DualFunctionalBase<dim>
  {
    public:
      PointValueEvaluation (const Point<dim> &evaluation_point);

      virtual
      void
      assemble_rhs (const DoFHandler<dim> &dof_handler,
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
  {};
  

  template <int dim>
  void
  PointValueEvaluation<dim>::
  assemble_rhs (const DoFHandler<dim> &dof_handler,
		Vector<double>        &rhs) const
  {
    rhs.reinit (dof_handler.n_dofs());
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    bool evaluation_point_found = false;
    for (; (cell!=endc) && !evaluation_point_found; ++cell)
      for (unsigned int vertex=0;
	   vertex<GeometryInfo<dim>::vertices_per_cell;
	   ++vertex)
	if (cell->vertex(vertex) == evaluation_point)
	  {
	    rhs(cell->vertex_dof_index(vertex,0)) = 1;

	    evaluation_point_found = true;
	    break;
	  };

    AssertThrow (evaluation_point_found,
		 ExcEvaluationPointNotFound(evaluation_point));
  };
  

};


namespace LaplaceSolver
{
  template <int dim>
  class DualSolver : public Solver<dim>
  {
    public:
      DualSolver (Triangulation<dim>       &triangulation,
		  const FiniteElement<dim> &fe,
		  const Quadrature<dim>    &quadrature,
		  const Quadrature<dim-1>  &face_quadrature,
		  const DualFunctional::DualFunctionalBase<dim> &dual_functional);

				       // XXX
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
  };

  template <int dim>
  const ZeroFunction<dim> DualSolver<dim>::boundary_values;

  template <int dim>
  DualSolver<dim>::
  DualSolver (Triangulation<dim>       &triangulation,
	      const FiniteElement<dim> &fe,
	      const Quadrature<dim>    &quadrature,
	      const Quadrature<dim-1>  &face_quadrature,
	      const DualFunctional::DualFunctionalBase<dim> &dual_functional)
		  :
		  Base<dim> (triangulation),
		  Solver<dim> (triangulation, fe,
			       quadrature, face_quadrature,
			       boundary_values),
                  dual_functional (&dual_functional)
  {};


  template <int dim>
  void
  DualSolver<dim>::solve_problem ()
  {
    Solver<dim>::solve_problem ();
  };



  template <int dim>
  unsigned int
  DualSolver<dim>::n_dofs() const
  {
    return Solver<dim>::n_dofs();
  };


  template <int dim>
  void
  DualSolver<dim>::
  postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    return Solver<dim>::postprocess(postprocessor);
  };
  


  template <int dim>
  void
  DualSolver<dim>::
  assemble_rhs (Vector<double> &rhs) const 
  {
    dual_functional->assemble_rhs (dof_handler, rhs);
  };


  template <int dim>
  class WeightedResidual : public PrimalSolver<dim>,
			   public DualSolver<dim>
  {
    public:
      WeightedResidual (Triangulation<dim>       &coarse_grid,
			const FiniteElement<dim> &primal_fe,
			const FiniteElement<dim> &dual_fe,
			const Quadrature<dim>    &quadrature,
			const Quadrature<dim-1>  &face_quadrature,
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

				       /**
					* Declare a data type to
					* represent the mapping between
					* faces and integrated jumps of
					* gradients of each of the
					* solution vectors. See the
					* general documentation of this
					* class for more information.
					*/
      typedef typename std::pair<double,typename DoFHandler<dim>::active_cell_iterator> FaceEntry;
      typedef typename std::map<typename DoFHandler<dim>::face_iterator,FaceEntry> FaceIntegrals;


				       /**
					* Redeclare an active cell iterator.
					* This is simply for convenience.
					*/
      typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;

				       /**
				      * All data needed by the several
				      * functions of the error
				      * estimator is gathered in this
				      * struct. It is passed as a
				      * reference to the separate
				      * functions in this class.
				      *
				      * The reason for invention of
				      * this object is two-fold:
				      * first, class member data is
				      * not possible because no real
				      * object is created (all
				      * functions are @p{static}),
				      * which is a historical
				      * reason. Second, if we don't
				      * collect the data the various
				      * functions need somewhere at a
				      * central place, that would mean
				      * that the functions would have
				      * to allocate them upon
				      * need. However, then some
				      * variables would be allocated
				      * over and over again, which can
				      * take a significant amount of
				      * time (10-20 per cent) and most
				      * importantly, memory allocation
				      * requires synchronisation in
				      * multithreaded mode. While that
				      * is done by the C++ library and
				      * has not to be handcoded, it
				      * nevertheless seriously damages
				      * the ability to efficiently run
				      * the functions of this class in
				      * parallel, since they are quite
				      * often blocked by these
				      * synchronisation points.
				      *
				      * Thus, every thread gets an
				      * instance of this class to work
				      * with and needs not allocate
				      * memory itself, or synchronise
				      * with other threads.
				      */
      struct FaceData
      {
	  FEFaceValues<dim>    fe_face_values_cell;
	  FEFaceValues<dim>    fe_face_values_neighbor;
	  FESubfaceValues<dim> fe_subface_values_cell;

	  std::vector<double> jump_residual;
	  std::vector<double> dual_weights;	  
	  typename std::vector<Tensor<1,dim> > cell_grads;
	  typename std::vector<Tensor<1,dim> > neighbor_grads;
	  FaceData (const FiniteElement<dim> &dof_handler,
		    const Quadrature<dim-1>  &face_quadrature);
      };

      struct CellData
      {
	  FEValues<dim>    fe_values;
	  const SmartPointer<const Function<dim> > right_hand_side;

	  std::vector<double> cell_residual;
	  std::vector<double> rhs_values;	  
	  std::vector<double> dual_weights;	  
	  typename std::vector<Tensor<2,dim> > cell_grad_grads;
	  CellData (const FiniteElement<dim> &dof_handler,
		    const Quadrature<dim>    &quadrature,
		    const Function<dim>      &right_hand_side);
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
      
				       /**
					* Actually do the computation on
					* a face which has no hanging
					* nodes (it is regular), i.e.
					* either on the other side there
					* is nirvana (face is at
					* boundary), or the other side's
					* refinement level is the same
					* as that of this side, then
					* handle the integration of
					* these both cases together.
					*
					* The meaning of the parameters
					* becomes clear when looking at
					* the source code. This function
					* is only externalized from
					* @p{estimate_error} to avoid
					* ending up with a function of
					* 500 lines of code.
					*/
      void
      integrate_over_regular_face (const active_cell_iterator &cell,
				   const unsigned int          face_no,
				   const Vector<double>       &primal_solution,
				   const Vector<double>       &dual_weights,
				   FaceData                   &face_data,
				   FaceIntegrals              &face_integrals) const;
      

				       /**
					* The same applies as for the
					* function above, except that
					* integration is over face
					* @p{face_no} of @p{cell}, where
					* the respective neighbor is
					* refined, so that the
					* integration is a bit more
					* complex.
					*/
      void
      integrate_over_irregular_face (const active_cell_iterator &cell,
				     const unsigned int          face_no,
				     const Vector<double>       &primal_solution,
				     const Vector<double>       &dual_weights,
				     FaceData                   &face_data,
				     FaceIntegrals              &face_integrals) const;
  };





  template <int dim>
  WeightedResidual<dim>::FaceData::
  FaceData (const FiniteElement<dim> &fe,
	    const Quadrature<dim-1>  &face_quadrature)
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
      = face_quadrature.n_quadrature_points;
  
    jump_residual.resize(n_face_q_points);
    dual_weights.resize(n_face_q_points);    
    cell_grads.resize(n_face_q_points);
    neighbor_grads.resize(n_face_q_points);
  };
  


  template <int dim>
  WeightedResidual<dim>::CellData::
  CellData (const FiniteElement<dim> &fe,
	    const Quadrature<dim>    &quadrature,
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
      = quadrature.n_quadrature_points;
  
    cell_residual.resize(n_q_points);
    rhs_values.resize(n_q_points);    
    dual_weights.resize(n_q_points);    
    cell_grad_grads.resize(n_q_points);
  };
  
  


  template <int dim>
  WeightedResidual<dim>::
  WeightedResidual (Triangulation<dim>       &coarse_grid,
		    const FiniteElement<dim> &primal_fe,
		    const FiniteElement<dim> &dual_fe,
		    const Quadrature<dim>    &quadrature,
		    const Quadrature<dim-1>  &face_quadrature,
		    const Function<dim>      &rhs_function,
		    const Function<dim>      &boundary_values,
		    const DualFunctional::DualFunctionalBase<dim> &dual_functional)
		  :
		  Base<dim> (coarse_grid),
                  PrimalSolver<dim> (coarse_grid, primal_fe,
				     quadrature, face_quadrature,
				     rhs_function, boundary_values),
                  DualSolver<dim> (coarse_grid, dual_fe,
				   quadrature, face_quadrature,
				   dual_functional)
  {};


  template <int dim>
  void
  WeightedResidual<dim>::solve_problem ()
  {
    PrimalSolver<dim>::solve_problem ();
    DualSolver<dim>::solve_problem ();
  };


  template <int dim>
  void
  WeightedResidual<dim>::postprocess (const Evaluation::EvaluationBase<dim> &postprocessor) const
  {
    PrimalSolver<dim>::postprocess (postprocessor);
  };
  
  
  template <int dim>
  unsigned int
  WeightedResidual<dim>::n_dofs () const
  {
    return PrimalSolver<dim>::n_dofs();
  };



  template <int dim>
  void
  WeightedResidual<dim>::refine_grid ()
  {
    Vector<float> error_indicators (triangulation->n_active_cells());
    estimate_error (error_indicators);
    std::cout << "Estimated error="
	      << std::accumulate (error_indicators.begin(),
				  error_indicators.end(), 0.)
	      << std::endl;
    DataOut<dim> data_out;
    ofstream x("x");
    Vector<double> xe (error_indicators.begin(),
		       error_indicators.end());
    data_out.attach_dof_handler (DualSolver<dim>::dof_handler);
    data_out.add_data_vector (xe, "e");
    data_out.build_patches ();
    data_out.write_gnuplot (x);
    
    std::transform (error_indicators.begin(),
		    error_indicators.end(),
		    error_indicators.begin(),
		    &fabs);
    GridRefinement::refine_and_coarsen_fixed_number (*triangulation,
						     error_indicators,
						     0.3, 0.03);
    triangulation->execute_coarsening_and_refinement ();
  };
  

  template <int dim>
  void
  WeightedResidual<dim>::output_solution () const
  {
    for (unsigned int solution_type=0; solution_type<2; ++solution_type)
      {
	DataOut<dim> data_out;

	switch (solution_type)
	  {
	    case 0:
		  data_out.attach_dof_handler (PrimalSolver<dim>::dof_handler);
		  data_out.add_data_vector (PrimalSolver<dim>::solution,
					    "primal_solution");
		  break;
	    case 1:
		  data_out.attach_dof_handler (DualSolver<dim>::dof_handler);
		  data_out.add_data_vector (DualSolver<dim>::solution,
					    "dual_solution");
		  break;
	    default:
		  Assert (false, ExcInternalError());
	  };
	data_out.build_patches ();
  
#ifdef HAVE_STD_STRINGSTREAM
	std::ostringstream filename;
#else
	std::ostrstream filename;
#endif
	filename << "solution-"
		 << (solution_type == 0 ?
		     "primal-" : "dual-")
		 << refinement_cycle
		 << ".gnuplot"
		 << std::ends;
#ifdef HAVE_STD_STRINGSTREAM
	std::ofstream out (filename.str().c_str());
#else
	std::ofstream out (filename.str());
#endif
    
	data_out.write (out, DataOut<dim>::gnuplot);
      };
  };


				   // @sect3{Estimating errors}

				   // @sect4{Error estimation driver functions}
				   //
				   // As for the actual computation of
				   // error estimates, let's start
				   // with the function that drives
				   // all this, i.e. calls those
				   // functions that actually do the
				   // work, and finally collects the
				   // results.
  
  template <int dim>
  void
  WeightedResidual<dim>::
  estimate_error (Vector<float> &error_indicators) const
  {
				     // The first task in computing
				     // the error is to set up vectors
				     // that denote the primal
				     // solution, and the weights
				     // (z-z_h)=(z-I_hz), both in the
				     // finite element space for which
				     // we have computed the dual
				     // solution. For this, we have to
				     // interpolate the primal
				     // solution to the dual finite
				     // element space, and to subtract
				     // the interpolation of the
				     // computed dual solution to the
				     // primal finite element
				     // space. Fortunately, the
				     // library provides functions for
				     // these two actions.
    Vector<double> primal_solution (DualSolver<dim>::dof_handler.n_dofs());
    FETools::interpolate (PrimalSolver<dim>::dof_handler,
			  PrimalSolver<dim>::solution,
			  DualSolver<dim>::dof_handler,
			  primal_solution);    
    Vector<double> dual_weights (DualSolver<dim>::dof_handler.n_dofs());
    FETools::interpolation_difference (DualSolver<dim>::dof_handler,
				       DualSolver<dim>::solution,
				       *PrimalSolver<dim>::fe,
				       dual_weights);
    
				     // Then we set up a map between
				     // face iterators and their jump
				     // term contributions of faces to
				     // the error estimator. The
				     // reason is that we compute the
				     // jump terms only once, from one
				     // side of the face, and want to
				     // collect them only afterwards
				     // when looping over all cells a
				     // second time.
				     //
				     // We initialize this map already
				     // with a value of -1e20 for all
				     // faces, since this value will
				     // strike in the results if
				     // something should go wrong and
				     // we fail to compute the value
				     // for a face for some
				     // reason. Secondly, we
				     // initialize the map once before
				     // we branch to different threads
				     // since this way the map's
				     // structure is no more modified
				     // by the individual threads,
				     // only existing entries are set
				     // to new values. This relieves
				     // us from the necessity to
				     // synchronise the threads
				     // through a mutex each time they
				     // write to this map.
    FaceIntegrals face_integrals;
    for (active_cell_iterator cell=DualSolver<dim>::dof_handler.begin_active();
	 cell!=DualSolver<dim>::dof_handler.end();
	 ++cell)
      for (unsigned int face_no=0;
	   face_no<GeometryInfo<dim>::faces_per_cell;
	   ++face_no)
	face_integrals[cell->face(face_no)].first = -1e20;

				     // Then set up a vector with
				     // error indicators.  Reserve one
				     // slot for each cell and set it
				     // to zero.
    error_indicators.reinit (DualSolver<dim>::dof_handler
			     .get_tria().n_active_cells());

				     // Now start a number of threads
				     // which compute the error
				     // formula on parts of all the
				     // cells, and once they are all
				     // started wait until they have
				     // all finished:
    const unsigned int n_threads = multithread_info.n_default_threads;
    Threads::ThreadManager thread_manager;
    for (unsigned int i=0; i<n_threads; ++i)
      Threads::spawn (thread_manager,
		      Threads::encapsulate (&WeightedResidual<dim>::
					    estimate_some)
		      .collect_args (this,
				     primal_solution,
				     dual_weights,
				     n_threads, i,
				     error_indicators,
				     face_integrals));
    thread_manager.wait();

				     // Once the error contributions
				     // are computed, sum them up. For
				     // this, note that the cell terms
				     // are already set, and that only
				     // the edge terms need to be
				     // collected. For this, loop over
				     // all cells and their faces,
				     // make sure that the
				     // contributions of each of the
				     // faces are there, and add them
				     // up.
    unsigned int present_cell=0;  
    for (active_cell_iterator cell=DualSolver<dim>::dof_handler.begin_active();
	 cell!=DualSolver<dim>::dof_handler.end();
	 ++cell, ++present_cell)
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
	   ++face_no)
	{
	  Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(),
		 ExcInternalError());
	  if (true || (face_integrals[cell->face(face_no)].second
		       ==
		       cell))
	    error_indicators(present_cell)
	      += -0.5*face_integrals[cell->face(face_no)].first;
	  else
	    error_indicators(present_cell)
	      -= -0.5*face_integrals[cell->face(face_no)].first;
	};
  };


				   // @sect4{Estimating on a subset of cells}

				   // Next we have the function that
				   // is called to estimate the error
				   // on a subset of cells. The
				   // function may be called multiply
				   // if the library was configured to
				   // use multi-threading. Here it
				   // goes:
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
				     // At the beginning, we
				     // initialize two variables for
				     // each thread which may be
				     // running this function. The
				     // reason for these functions was
				     // discussed above, when the
				     // respective classes were
				     // discussed, so we here only
				     // point out that since they are
				     // local to the function that is
				     // spawned when running more than
				     // one thread, the data of these
				     // objects exists actually once
				     // per thread, so we don't have
				     // to take care about
				     // synchronising access to them.
    CellData cell_data (*DualSolver<dim>::fe,
			*DualSolver<dim>::quadrature,
			*PrimalSolver<dim>::rhs_function);
    FaceData face_data (*DualSolver<dim>::fe,
			*DualSolver<dim>::face_quadrature);    

				     // Then calculate the start cell
				     // for this thread. We let the
				     // different threads run on
				     // interleaved cells, i.e. for
				     // example if we have 4 threads,
				     // then the first thread treates
				     // cells 0, 4, 8, etc, while the
				     // second threads works on cells 1,
				     // 5, 9, and so on. The reason is
				     // that it takes vastly more time
				     // to work on cells with hanging
				     // nodes than on regular cells, but
				     // such cells are not evenly
				     // distributed across the range of
				     // cell iterators, so in order to
				     // have the different threads do
				     // approximately the same amount of
				     // work, we have to let them work
				     // interleaved to the effect of a
				     // pseudorandom distribution of the
				     // `hard' cells to the different
				     // threads.
    active_cell_iterator cell=DualSolver<dim>::dof_handler.begin_active();
    for (unsigned int t=0;
	 (t<this_thread) && (cell!=DualSolver<dim>::dof_handler.end());
	 ++t, ++cell);
  
				     // Next loop over all cells. The
				     // check for loop end is done at
				     // the end of the loop, along
				     // with incrementing the loop
				     // index.
    for (unsigned int cell_index=this_thread; true; )
      {
					 // First task on each cell is
					 // to compute the cell
					 // residual contributions of
					 // this cell, and put them
					 // into the
					 // ``error_indicators''
					 // variable:
	integrate_over_cell (cell, cell_index,
			     primal_solution,
			     dual_weights,
			     cell_data,
			     error_indicators);
	
					 // After computing the cell
					 // terms, turn to the face
					 // terms. For this, loop over
					 // all faces of the present
					 // cell, and see whether
					 // something needs to be
					 // computed on it:
	for (unsigned int face_no=0;
	     face_no<GeometryInfo<dim>::faces_per_cell;
	     ++face_no)
	  {
					     // First, if this face is
					     // part of the boundary,
					     // then there is nothing
					     // to do. However, to
					     // make things easier
					     // when summing up the
					     // contributions of the
					     // faces of cells, we
					     // enter this face into
					     // the list of faces with
					     // a zero contribution to
					     // the error, and also
					     // mark the cell on which
					     // we computed this
					     // value.
	    if (cell->face(face_no)->at_boundary()) 
	      {
		face_integrals[cell->face(face_no)].first = 0;
		face_integrals[cell->face(face_no)].second = cell;
		continue;
	      };
	    
					     // Next, note that since
					     // we want to compute the
					     // jump terms on each
					     // face only once,
					     // although we access it
					     // twice if it is not at
					     // the boundary, we have
					     // to define some rules
					     // who is responsible for
					     // computing on a face:
					     //
					     // First, if the
					     // neighboring cell is on
					     // the same level as this
					     // one, i.e. neither
					     // further refined not
					     // coarser, then the one
					     // with the lower index
					     // within this level does
					     // the work. In other
					     // words: if the other
					     // one has a lower index,
					     // then skip work on this
					     // face:
	    if ((cell->neighbor(face_no)->has_children() == false) &&
		(cell->neighbor(face_no)->level() == cell->level()) &&
		(cell->neighbor(face_no)->index() < cell->index()))
	      continue;

					     // Likewise, we always
					     // work from the coarser
					     // cell if this and its
					     // neighbor differ in
					     // refinement. Thus, if
					     // the neighboring cell
					     // is less refined than
					     // the present one, then
					     // do nothing since we
					     // integrate over the
					     // subfaces when we visit
					     // the coarse cell.
	    if (cell->at_boundary(face_no) == false)
	      if (cell->neighbor(face_no)->level() < cell->level())
		continue;	  


					     // Now we know that we
					     // are in charge here, so
					     // actually compute the
					     // face jump terms. If
					     // the face is a regular
					     // one, i.e.  the other
					     // side's cell is neither
					     // coarser not finer than
					     // this cell, then call
					     // one function, and if
					     // the cell on the other
					     // side is further
					     // refined, then use
					     // another function. Note
					     // that the case that the
					     // cell on the other side
					     // is coarser cannot
					     // happen since we have
					     // decided above that we
					     // handle this case when
					     // we pass over that
					     // other cell.
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

					 // After computing the cell
					 // contributions and looping
					 // over the faces, go to the
					 // next cell for this
					 // thread. Note again that
					 // the cells for each of the
					 // threads are
					 // interleaved. If we are at
					 // the end of our workload,
					 // jump out of the loop.
	for (unsigned int t=0;
	     ((t<n_threads) && (cell!=DualSolver<dim>::dof_handler.end()));
	     ++t, ++cell, ++cell_index);
	if (cell == DualSolver<dim>::dof_handler.end())
	  break;
      };    
  };


				   // @sect4{Computing cell term error contributions}

				   // As for the actual computation of
				   // the error contributions, first
				   // turn to the cell terms:
  template <int dim>
  void WeightedResidual<dim>::
  integrate_over_cell (const active_cell_iterator &cell,
		       const unsigned int          cell_index,
		       const Vector<double>       &primal_solution,
		       const Vector<double>       &dual_weights,
		       CellData                   &cell_data,
		       Vector<float>              &error_indicators) const
  {
				     // The tasks to be done are what
				     // appears natural from looking
				     // at the error estimation
				     // formula: first compute the the
				     // right hand side and the
				     // Laplacian of the numerical
				     // solution at the quadrature
				     // points for the cell residual,
    cell_data.fe_values.reinit (cell);
    cell_data.right_hand_side
      ->value_list (cell_data.fe_values.get_quadrature_points(),
		    cell_data.rhs_values);
    cell_data.fe_values.get_function_2nd_derivatives (primal_solution,
						      cell_data.cell_grad_grads);

				     // ...then get the dual weights...
    cell_data.fe_values.get_function_values (dual_weights,
					     cell_data.dual_weights);

				     // ...and finally build the sum
				     // over all quadrature points:
    double sum = 0;
    for (unsigned int p=0; p<cell_data.fe_values.n_quadrature_points; ++p)
      sum += ((cell_data.rhs_values[p]+trace(cell_data.cell_grad_grads[p])) *
	      cell_data.dual_weights[p] *
	      cell_data.fe_values.JxW (p));
    error_indicators(cell_index) += sum;
  };


				   // @sect4{Computing edge term error contributions - 1}
  
				   // On the other hand, computation
				   // of the edge terms for the error
				   // estimate is not so
				   // simple. First, we have to
				   // distinguish between faces with
				   // and without hanging
				   // nodes. Because it is the simple
				   // case, we first consider the case
				   // without hanging nodes on a face
				   // (let's call this the `regular'
				   // case):
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
      n_q_points = face_data.fe_face_values_cell.n_quadrature_points;

				     // The first step is to get the
				     // values of the gradients at the
				     // quadrature points of the
				     // finite element field on the
				     // present cell. For this,
				     // initialize the
				     // ``FEFaceValues'' object
				     // corresponding to this side of
				     // the face, and extract the
				     // gradients using that
				     // object.
    face_data.fe_face_values_cell.reinit (cell, face_no);
    face_data.fe_face_values_cell.get_function_grads (primal_solution,
						      face_data.cell_grads);

				     // The second step is then to
				     // extract the gradients of the
				     // finite element solution at the
				     // quadrature points on the other
				     // side of the face, i.e. from
				     // the neighboring cell.
				     //
				     // For this, do a sanity check
				     // before: make sure that the
				     // neigbor actually exists (yes,
				     // we should not have come here
				     // if the neighbor did not exist,
				     // but in complicated software
				     // there are bugs, so better
				     // check this), and if this is
				     // not the case throw an error.
    Assert (cell->neighbor(face_no).state() == IteratorState::valid,
	    ExcInternalError());
				     // If we have that, then we need
				     // to find out with which face of
				     // the neighboring cell we have
				     // to work, i.e. the
				     // ``home-many''the neighbor the
				     // present cell is of the cell
				     // behind the present face. For
				     // this, there is a function, and
				     // we put the result into a
				     // variable with the name
				     // ``neighbor_neighbor'':
    const unsigned int
      neighbor_neighbor = cell->neighbor_of_neighbor (face_no);
				     // Then define an abbreviation
				     // for the neigbor cell,
				     // initialize the
				     // ``FEFaceValues'' object on
				     // that cell, and extract the
				     // gradients on that cell:
    const active_cell_iterator neighbor = cell->neighbor(face_no);
    face_data.fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor);      
    face_data.fe_face_values_neighbor.get_function_grads (primal_solution,
							  face_data.neighbor_grads);

				     // Now that we have the gradients
				     // on this and the neighboring
				     // cell, compute the jump
				     // residual by multiplying the
				     // jump in the gradient with the
				     // normal vector:
    for (unsigned int p=0; p<n_q_points; ++p)
      face_data.jump_residual[p]
	= ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) *
	   face_data.fe_face_values_cell.normal_vector(p));

				     // Next get the dual weights for
				     // this face:
    face_data.fe_face_values_cell.get_function_values (dual_weights,
						       face_data.dual_weights);
    
				     // Finally, we have to compute
				     // the sum over jump residuals,
				     // dual weights, and quadrature
				     // weights, to get the result for
				     // this face:
    double face_integral = 0;
    for (unsigned int p=0; p<n_q_points; ++p)
      face_integral += (face_data.jump_residual[p] *
			face_data.dual_weights[p]  *
			face_data.fe_face_values_cell.JxW(p));

				     // Double check that the element
				     // already exists and that it was
				     // not already written to...
    Assert (face_integrals.find (cell->face(face_no)) != face_integrals.end(),
	    ExcInternalError());
    Assert (face_integrals[cell->face(face_no)].first == -1e20,
	    ExcInternalError());

				     // ...then store computed value
				     // at assigned location:
    face_integrals[cell->face(face_no)].first = face_integral;
    face_integrals[cell->face(face_no)].second = cell;
  };


  				   // @sect4{Computing edge term error contributions - 2}
  
				   // We are still missing the case of
				   // faces with hanging nodes. This
				   // is what is covered in this
				   // function:
  template <int dim>
  void WeightedResidual<dim>::
  integrate_over_irregular_face (const active_cell_iterator &cell,
				 const unsigned int          face_no,
				 const Vector<double>       &primal_solution,
				 const Vector<double>       &dual_weights,
				 FaceData                   &face_data,
				 FaceIntegrals              &face_integrals) const
  {
				     // First again two abbreviations,
				     // and some consistency checks
				     // whether the function is called
				     // only on faces for which it is
				     // supposed to be called:
    const unsigned int
      n_q_points = face_data.fe_face_values_cell.n_quadrature_points;

    const typename DoFHandler<dim>::cell_iterator
      neighbor = cell->neighbor(face_no);    
    Assert (neighbor.state() == IteratorState::valid,
	    ExcInternalError());
    Assert (neighbor->has_children(),
	    ExcInternalError());

				     // Then find out which neighbor
				     // the present cell is of the
				     // adjacent cell. Note that we
				     // will operator on the children
				     // of this adjacent cell, but
				     // that their orientation is the
				     // same as that of their mother,
				     // i.e. the neigbor direction is
				     // the same.
    const unsigned int
      neighbor_neighbor = cell->neighbor_of_neighbor (face_no);
  
				     // Then simply do everything we
				     // did in the previous function
				     // for one face for all the
				     // sub-faces now:
    for (unsigned int subface_no=0;
	 subface_no<GeometryInfo<dim>::subfaces_per_face;
	 ++subface_no)
      {
					 // Start with some checks
					 // again: get an iterator
					 // pointing to the cell
					 // behind the present subface
					 // and check whether its face
					 // is a subface of the one we
					 // are considering. If that
					 // were not the case, then
					 // there would be either a
					 // bug in the
					 // ``neighbor_neighbor''
					 // function called above, or
					 // -- worse -- some function
					 // in the library did not
					 // keep to some underlying
					 // assumptions about cells,
					 // their children, and their
					 // faces. In any case, even
					 // though this assertion
					 // should not be triggered,
					 // it does not harm to be
					 // cautious, and in optimized
					 // mode computations the
					 // assertion will be removed
					 // anyway.
	const active_cell_iterator neighbor_child
	  = neighbor->child(GeometryInfo<dim>::
			    child_cell_on_face(neighbor_neighbor,
					       subface_no));
	Assert (neighbor_child->face(neighbor_neighbor) ==
		cell->face(face_no)->child(subface_no),
		ExcInternalError());

					 // Now start the work by
					 // again getting the gradient
					 // of the solution first at
					 // this side of the
					 // interface,
	face_data.fe_subface_values_cell.reinit (cell, face_no, subface_no);
	face_data.fe_subface_values_cell.get_function_grads (primal_solution,
							     face_data.cell_grads);
					 // then at the other side,
	face_data.fe_face_values_neighbor.reinit (neighbor_child,
					     neighbor_neighbor);
	face_data.fe_face_values_neighbor.get_function_grads (primal_solution,
							      face_data.neighbor_grads);
      
					 // and finally building the
					 // jump residuals. Since we
					 // take the normal vector
					 // from the other cell this
					 // time, revert the sign of
					 // the first term compared to
					 // the other function:
	for (unsigned int p=0; p<n_q_points; ++p)
	  face_data.jump_residual[p]
	     = ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) *
		face_data.fe_face_values_neighbor.normal_vector(p));

					 // Then get dual weights:
	face_data.fe_face_values_neighbor.get_function_values (dual_weights,
							       face_data.dual_weights);
	
					 // At last, sum up the
					 // contribution of this
					 // sub-face, and set it in
					 // the global map:
	double face_integral = 0;
	for (unsigned int p=0; p<n_q_points; ++p)
	  face_integral += (face_data.jump_residual[p] *
			    face_data.dual_weights[p] *
			    face_data.fe_face_values_neighbor.JxW(p));
	face_integrals[neighbor_child->face(neighbor_neighbor)].first
	  = face_integral;
	face_integrals[neighbor_child->face(neighbor_neighbor)].second
	  = cell;
      };

				     // Once the contributions of all
				     // sub-faces are computed, loop
				     // over all sub-faces to collect
				     // and store them with the mother
				     // face for simple use when later
				     // collecting the error terms of
				     // cells. Again make safety
				     // checks that the entries for
				     // the sub-faces have been
				     // computed and do not carry an
				     // invalid value.
    double sum = 0;
    typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
    for (unsigned int subface_no=0;
	 subface_no<GeometryInfo<dim>::subfaces_per_face;
	 ++subface_no) 
      {
	Assert (face_integrals.find(face->child(subface_no)) !=
		face_integrals.end(),
		ExcInternalError());
	Assert (face_integrals[face->child(subface_no)].first != -1e20,
		ExcInternalError());
      
	sum += face_integrals[face->child(subface_no)].first;
      };
				     // Finally store the value with
				     // the parent face.
    face_integrals[face].first = sum;
    face_integrals[face].second = cell;    
  };
  
};




template <int dim>
void
run_simulation (LaplaceSolver::Base<dim>                     &solver,
		const std::list<Evaluation::EvaluationBase<dim> *> &postprocessor_list)
{
  std::cout << "Refinement cycle: ";

  for (unsigned int step=0; true; ++step)
    {
      std::cout << step << " Solving "
		<< solver.n_dofs()
		<< std::endl;

      solver.set_refinement_cycle (step);
      solver.solve_problem ();
      solver.output_solution ();

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

  std::cout << std::endl;
};



void
create_triangulation (Triangulation<2> &tria)
{
  const Point<2>
    vertices[16] = {  Point<2> (-1.,   -1.),
		      Point<2> (-1./3, -1.),
		      Point<2> (+1./3, -1.),
		      Point<2> (+1,    -1.),
		      Point<2> (-1.,   -1./3.),
		      Point<2> (-1./3, -1./3.),
		      Point<2> (+1./3, -1./3.),
		      Point<2> (+1,    -1./3.),
		      Point<2> (-1.,   1./3.),
		      Point<2> (-1./3, 1./3.),
		      Point<2> (+1./3, 1./3.),
		      Point<2> (+1,    1./3.),
		      Point<2> (-1.,   1.),
		      Point<2> (-1./3, 1.),
		      Point<2> (+1./3, 1.),
		      Point<2> (+1,    1.)    };
  
  const int cell_vertices[8][4] = {{0, 1, 5, 4},
				   {1, 2, 6, 5},
				   {2, 3, 7, 6},
				   {4, 5, 9, 8},
				   {6, 7, 11, 10},
				   {8,9,13,12},
				   {9,10,14,13},
				   {10,11,15,14}};

  std::vector<CellData<2> > cells (8, CellData<2>());
  
  for (unsigned int i=0; i<8; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (std::vector<Point<2> >(&vertices[0], &vertices[16]),
			     cells,
			     SubCellData());       // no boundary information
};



template <int dim>
void solve_problem (const std::string &solver_name) 
{
  const std::string header = "Running tests with \"" + solver_name +
			     "\" refinement criterion:";
  std::cout << header << std::endl
	    << std::string (header.size(), '-') << std::endl;

  Triangulation<dim> triangulation;
//  create_triangulation (triangulation);
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (2);
  const FE_Q<dim>          primal_fe(1);
  const FE_Q<dim>          dual_fe(2);
  const QGauss4<dim>       quadrature;
  const QGauss4<dim-1>     face_quadrature;  
  const RightHandSide<dim> rhs_function;
  const Solution<dim>      boundary_values;

  const Point<dim> evaluation_point(0.5,0.5);
  const DualFunctional::PointValueEvaluation<dim>
    dual_functional (evaluation_point);
  
  LaplaceSolver::Base<dim> * solver = 0;
  solver = new LaplaceSolver::WeightedResidual<dim> (triangulation,
						     primal_fe,
						     dual_fe,
						     quadrature,
						     face_quadrature,
						     rhs_function,
						     boundary_values,
						     dual_functional);

  TableHandler results_table;
  Evaluation::PointValueEvaluation<dim>
    postprocessor1 (Point<dim>(0.5,0.5), results_table);

  std::list<Evaluation::EvaluationBase<dim> *> postprocessor_list;
  postprocessor_list.push_back (&postprocessor1);

  run_simulation (*solver, postprocessor_list);

  results_table.write_text (std::cout);
  delete solver;

  std::cout << std::endl;
};



int main () 
{
  try
    {
      deallog.depth_console (0);

      solve_problem<2> ("global");
//      solve_problem<2> ("kelly");      
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
};
