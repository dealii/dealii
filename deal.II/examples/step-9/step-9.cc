/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_bicgstab.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <fe/fe_lib.lagrange.h>
#include <grid/grid_out.h>

#ifdef DEAL_II_USE_MT
#  include <base/thread_management.h>
#  include <base/multithread_info.h>
#endif

#include <dofs/dof_constraints.h>

#include <numerics/error_estimator.h>

#include <fstream>



// in strict ANSI C mode, the following constants are not defined by
// default, so we do it ourselves
#ifndef M_PI
#  define	M_PI		3.14159265358979323846
#endif



template <int dim>
class AdvectionProblem 
{
  public:
    AdvectionProblem ();
    ~AdvectionProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void assemble_system_interval (const DoFHandler<dim>::active_cell_iterator &begin,
				   const DoFHandler<dim>::active_cell_iterator &begin);
    
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FEQ1<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;

#ifdef DEAL_II_USE_MT
    ACE_Thread_Mutex     assembler_lock;
#endif
};



template <int dim>
class AdvectionField 
{
  public:
    Point<dim> value (const Point<dim> &p) const;
    
    void value_list (const vector<Point<dim> > &points,
		     vector<Point<dim> >       &values) const;

				     /**
				      * Exception
				      */
    DeclException2 (ExcDimensionMismatch,
		    int, int,
		    << "The vector has size " << arg1 << " but should have "
		    << arg2 << " elements.");
};



template <int dim>
Point<dim>
AdvectionField<dim>::value (const Point<dim> &p) const 
{
  Point<dim> value;
  value[0] = 2;
  for (unsigned int i=1; i<dim; ++i)
    value[i] = 1+0.8*sin(8*M_PI*p[0]);

  return value;
};



template <int dim>
void
AdvectionField<dim>::value_list (const vector<Point<dim> > &points,
				 vector<Point<dim> >       &values) const 
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch (values.size(), points.size()));
  
  for (unsigned int i=0; i<points.size(); ++i)
    values[i] = AdvectionField<dim>::value (points[i]);
};




template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;
    
  private:
    static const Point<dim> center_point;
};


template <>
const Point<1> RightHandSide<1>::center_point = Point<1> (-0.75);

template <>
const Point<2> RightHandSide<2>::center_point = Point<2> (-0.75, -0.75);

template <>
const Point<3> RightHandSide<3>::center_point = Point<3> (-0.75, -0.75, -0.75);



template <int dim>
double
RightHandSide<dim>::value (const Point<dim>   &p,
			   const unsigned int  component) const 
{
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  const double diameter = 0.1;
  return ( (p-center_point).square() < diameter*diameter ?
	   .1/pow(diameter,dim) :
	   0);
};



template <int dim>
void
RightHandSide<dim>::value_list (const vector<Point<dim> > &points,
				vector<double>            &values,
				const unsigned int         component) const 
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch (values.size(), points.size()));
  
  for (unsigned int i=0; i<points.size(); ++i)
    values[i] = RightHandSide<dim>::value (points[i], component);
};



template <int dim>
class BoundaryValues : public Function<dim>
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    virtual void value_list (const vector<Point<dim> > &points,
			     vector<double>            &values,
			     const unsigned int         component = 0) const;
};



template <int dim>
double
BoundaryValues<dim>::value (const Point<dim>   &p,
			    const unsigned int  component) const 
{
  Assert (component == 0, ExcIndexRange (component, 0, 1));

  const double sine_term = sin(16*M_PI*sqrt(p.square()));
  const double weight    = exp(-5*p.square()) / exp(-5.);
  return sine_term * weight;
};



template <int dim>
void
BoundaryValues<dim>::value_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
				 const unsigned int         component) const 
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch (values.size(), points.size()));
  
  for (unsigned int i=0; i<points.size(); ++i)
    values[i] = BoundaryValues<dim>::value (points[i], component);
};



template <int dim>
AdvectionProblem<dim>::AdvectionProblem () :
		dof_handler (triangulation)
{};



template <int dim>
AdvectionProblem<dim>::~AdvectionProblem () 
{
  dof_handler.clear ();
};



template <int dim>
void AdvectionProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  hanging_node_constraints.condense (sparsity_pattern);

  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



template <int dim>
void AdvectionProblem<dim>::assemble_system () 
{
#ifdef DEAL_II_USE_MT
  const unsigned int n_threads = multithread_info.n_default_threads;
  ACE_Thread_Manager thread_manager;

				   // define starting and end point
				   // for each thread
  typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;
  vector<pair<active_cell_iterator,active_cell_iterator> >
    thread_ranges = Threads::split_range<active_cell_iterator> (dof_handler.begin_active (),
								dof_handler.end (),
								n_threads);

  for (unsigned int thread=0; thread<n_threads; ++thread)
    Threads::spawn (thread_manager,
		    Threads::encapsulate(&AdvectionProblem<2>::assemble_system_interval)
		    .collect_args (this,
				   thread_ranges[thread].first,
				   thread_ranges[thread].second));
  thread_manager.wait ();  
#else
  assemble_system_interval (dof_handler.begin_active(),
			    dof_handler.end());
#endif
};



template <int dim>
void
AdvectionProblem<dim>::
assemble_system_interval (const DoFHandler<dim>::active_cell_iterator &begin,
			  const DoFHandler<dim>::active_cell_iterator &end) 
{
  AdvectionField<dim> advection_field;
  RightHandSide<dim>  right_hand_side;
  BoundaryValues<dim> boundary_values;
  
  QGauss3<dim>   quadrature_formula;
  QGauss3<dim-1> face_quadrature_formula;
  
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
				    UpdateFlags (update_values     |
						 update_q_points   |
						 update_JxW_values |
						 update_normal_vectors));

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.n_quadrature_points;
  const unsigned int   n_face_q_points = face_quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  vector<unsigned int> local_dof_indices (dofs_per_cell);

  vector<double>       rhs_values (n_q_points);
  vector<Point<dim> >  advection_directions (n_q_points);
  vector<double>       face_boundary_values (n_face_q_points);
  vector<Point<dim> >  face_advection_directions (n_face_q_points);

  DoFHandler<dim>::active_cell_iterator cell;
  for (cell=begin; cell!=end; ++cell)
    {
      cell_matrix.clear ();
      cell_rhs.clear ();

      fe_values.reinit (cell);
      const FullMatrix<double> 
	& shape_values = fe_values.get_shape_values();
      const vector<vector<Tensor<1,dim> > >
	& shape_grads  = fe_values.get_shape_grads();
      const vector<double>
	& JxW_values   = fe_values.get_JxW_values();
      const vector<Point<dim> >
	& q_points     = fe_values.get_quadrature_points();

      advection_field.value_list (q_points, advection_directions);
      right_hand_side.value_list (q_points, rhs_values);

      const double delta = 0.1 * cell->diameter ();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += ((advection_directions[q_point] *
				    shape_grads[j][q_point]    *
				    (shape_values(i,q_point) +
				     delta *
				     advection_directions[q_point] *
				     shape_grads[i][q_point])) *
				   JxW_values[q_point]);

	    cell_rhs(i) += ((shape_values (i,q_point) +
			     delta *
			     advection_directions[q_point] *
			     shape_grads[i][q_point]        ) *
			    rhs_values[i] *
			    fe_values.JxW (q_point));
	  };

      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	{
	  fe_face_values.reinit (cell, face);
	  
	  const FullMatrix<double> 
	    & face_shape_values = fe_face_values.get_shape_values();
	  const vector<double>
	    & face_JxW_values   = fe_face_values.get_JxW_values();
	  const vector<Point<dim> >
	    & face_q_points     = fe_face_values.get_quadrature_points();
	  const vector<Point<dim> >
	    & normal_vectors    = fe_face_values.get_normal_vectors();

	  boundary_values.value_list (face_q_points, face_boundary_values);
	  advection_field.value_list (face_q_points, face_advection_directions);

	  for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
	    if (cell->face(face)->at_boundary () &&
		(normal_vectors[q_point] * face_advection_directions[q_point] < 0))
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    cell_matrix(i,j) -= (face_advection_directions[q_point] *
					 normal_vectors[q_point] *
					 face_shape_values(i,q_point) *
					 face_shape_values(j,q_point) *
					 face_JxW_values[q_point]);
		  
		  cell_rhs(i) -= (face_advection_directions[q_point] *
				  normal_vectors[q_point] *
				  face_boundary_values[q_point] *
				  face_shape_values(i,q_point) *
				  face_JxW_values[q_point]);
		};
	};
	  

      cell->get_dof_indices (local_dof_indices);
#ifdef DEAL_II_USE_MT
      assembler_lock.acquire ();
#endif
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	};
#ifdef DEAL_II_USE_MT
      assembler_lock.release ();
#endif
    };

  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
				   // no bdr val
};



template <int dim>
void AdvectionProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverBicgstab<>        bicgstab (solver_control, vector_memory);

  PreconditionJacobi<> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);

  bicgstab.solve (system_matrix, solution, system_rhs,
		  preconditioner);

  hanging_node_constraints.distribute (solution);
};


template <int dim>
void AdvectionProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::FunctionMap neumann_boundary;

  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss3<dim-1>(),
				      neumann_boundary,
				      solution,
				      estimated_error_per_cell);

  triangulation.refine_and_coarsen_fixed_number (estimated_error_per_cell,
						 0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
};



template <int dim>
void AdvectionProblem<dim>::output_results (const unsigned int cycle) const
{
  string filename = "grid-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".eps";
  ofstream output (filename.c_str());

  GridOut grid_out;
  grid_out.write_eps (triangulation, output);
};



template <int dim>
void AdvectionProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      cout << "Cycle " << cycle << ':' << endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation, -1, 1);
	  triangulation.refine_global (4);
	}
      else
	{
	  refine_grid ();
	};
      

      cout << "   Number of active cells:       "
	   << triangulation.n_active_cells()
	   << endl;

      setup_system ();

      cout << "   Number of degrees of freedom: "
	   << dof_handler.n_dofs()
	   << endl;
      
      assemble_system ();
      solve ();
      output_results (cycle);
    };

  typename DataOut<dim>::EpsFlags eps_flags;
  eps_flags.z_scaling = 4;
  
  DataOut<dim> data_out;
  data_out.set_flags (eps_flags);

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  
  ofstream output ("final-solution.gmv");
  data_out.write_gmv (output);
};



				 // The ``main'' function is exactly
				 // like in previous examples, with
				 // the only difference in the name of
				 // the main class that actually does
				 // the computation.
int main () 
{
  try
    {
      deallog.depth_console (0);

      AdvectionProblem<2> advection_problem_2d;
      advection_problem_2d.run ();
    }
  catch (exception &exc)
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Exception on processing: " << endl
	   << exc.what() << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
      return 1;
    }
  catch (...) 
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Unknown exception!" << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
      return 1;
    };

  return 0;
};
