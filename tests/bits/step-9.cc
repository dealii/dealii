//----------------------------  step-9.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  step-9.cc  ---------------------------


// a un-hp-ified version of hp/step-9


#include "../tests.h"
#include <base/logstream.h>
#include <fstream>
std::ofstream logfile("step-9/output");


#include <base/quadrature_lib.h>
#include <base/function.h>
#include "../tests.h"
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_bicgstab.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <hp/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <lac/constraint_matrix.h>
#include <hp/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <fe/fe_q.h>
#include <grid/grid_out.h>

#include <base/thread_management.h>
#include <base/multithread_info.h>

#include <base/tensor_function.h>

#include <numerics/error_estimator.h>

#include <fstream>
#include <iomanip>



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
    void assemble_system_interval (const typename DoFHandler<dim>::active_cell_iterator &begin,
			 	   const typename DoFHandler<dim>::active_cell_iterator &end);
    
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FE_Q<dim>            fe;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    Threads::ThreadMutex     assembler_lock;
};




template <int dim>
class AdvectionField : public TensorFunction<1,dim>
{
  public:
    AdvectionField () : TensorFunction<1,dim> () {}
    
    virtual Tensor<1,dim> value (const Point<dim> &p) const;
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<Tensor<1,dim> >    &values) const;

    DeclException2 (ExcDimensionMismatch,
		    unsigned int, unsigned int,
		    << "The vector has size " << arg1 << " but should have "
		    << arg2 << " elements.");
};



template <int dim>
Tensor<1,dim> 
AdvectionField<dim>::value (const Point<dim> &p) const 
{
  Point<dim> value;
  value[0] = 2;
  for (unsigned int i=1; i<dim; ++i)
    value[i] = 1+0.8*std::sin(8*M_PI*p[0]);

  return value;
}



template <int dim>
void
AdvectionField<dim>::value_list (const std::vector<Point<dim> > &points,
				 std::vector<Tensor<1,dim> >    &values) const 
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch (values.size(), points.size()));
  
  for (unsigned int i=0; i<points.size(); ++i)
    values[i] = AdvectionField<dim>::value (points[i]);
}




template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;
    
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
	   .1/std::pow(diameter,dim) :
	   0);
}



template <int dim>
void
RightHandSide<dim>::value_list (const std::vector<Point<dim> > &points,
				std::vector<double>            &values,
				const unsigned int              component) const 
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch (values.size(), points.size()));
  
  for (unsigned int i=0; i<points.size(); ++i)
    values[i] = RightHandSide<dim>::value (points[i], component);
}



template <int dim>
class BoundaryValues : public Function<dim>
{
  public:
    BoundaryValues () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;
};



template <int dim>
double
BoundaryValues<dim>::value (const Point<dim>   &p,
			    const unsigned int  component) const 
{
  Assert (component == 0, ExcIndexRange (component, 0, 1));

  const double sine_term = std::sin(16*M_PI*std::sqrt(p.square()));
  const double weight    = std::exp(-5*p.square()) / std::exp(-5.);
  return sine_term * weight;
}



template <int dim>
void
BoundaryValues<dim>::value_list (const std::vector<Point<dim> > &points,
				 std::vector<double>            &values,
				 const unsigned int              component) const 
{
  Assert (values.size() == points.size(),
	  ExcDimensionMismatch (values.size(), points.size()));
  
  for (unsigned int i=0; i<points.size(); ++i)
    values[i] = BoundaryValues<dim>::value (points[i], component);
}




class GradientEstimation
{
  public:
    template <int dim>
    static void estimate (const DoFHandler<dim> &dof,
			  const Vector<double>  &solution,
			  Vector<float>         &error_per_cell);

    DeclException2 (ExcInvalidVectorLength,
		    int, int,
		    << "Vector has length " << arg1 << ", but should have "
		    << arg2);
    DeclException0 (ExcInsufficientDirections);

  private:
    typedef std::pair<unsigned int,unsigned int> IndexInterval;

    template <int dim>
    static void estimate_interval (const DoFHandler<dim> &dof,
				   const Vector<double>  &solution,
				   const IndexInterval   &index_interval,
				   Vector<float>         &error_per_cell);    
};





template <int dim>
AdvectionProblem<dim>::AdvectionProblem () :
		dof_handler (triangulation),
		fe(1)
{}



template <int dim>
AdvectionProblem<dim>::~AdvectionProblem () 
{
  dof_handler.clear ();
}



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
}



template <int dim>
void AdvectionProblem<dim>::assemble_system () 
{
  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadGroup<> threads;

  typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> >
    thread_ranges 
    = Threads::split_range<active_cell_iterator> (dof_handler.begin_active (),
						  dof_handler.end (),
						  n_threads);

  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (*this, &AdvectionProblem<dim>::assemble_system_interval)
               (thread_ranges[thread].first,
                thread_ranges[thread].second);

  threads.join_all ();  


  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
}


 
template <int dim>
void
AdvectionProblem<dim>::
assemble_system_interval (const typename DoFHandler<dim>::active_cell_iterator &begin,
			  const typename DoFHandler<dim>::active_cell_iterator &end) 
{
  const AdvectionField<dim> advection_field;
  const RightHandSide<dim>  right_hand_side;
  const BoundaryValues<dim> boundary_values;
  
  QGauss<dim>   quadrature_formula (2);
  QGauss<dim-1> face_quadrature_formula (2);
  
  FEValues<dim> x_fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients |
                           update_q_points | update_JxW_values);
  FEFaceValues<dim> x_fe_face_values (fe, face_quadrature_formula,
				    update_values     | update_q_points   |
                                    update_JxW_values | update_normal_vectors);

  const unsigned int   dofs_per_cell   = fe[0].dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.n_quadrature_points;
  const unsigned int   n_face_q_points = face_quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  std::vector<double>         rhs_values (n_q_points);
  std::vector<Tensor<1,dim> > advection_directions (n_q_points);
  std::vector<double>         face_boundary_values (n_face_q_points);
  std::vector<Tensor<1,dim> > face_advection_directions (n_face_q_points);

  typename DoFHandler<dim>::active_cell_iterator cell;
  for (cell=begin; cell!=end; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      x_fe_values.reinit (cell);
      const FEValues<dim> &fe_values = x_fe_values.get_present_fe_values();

      advection_field.value_list (fe_values.get_quadrature_points(),
				  advection_directions);
      right_hand_side.value_list (fe_values.get_quadrature_points(),
				  rhs_values);

      const double delta = 0.1 * cell->diameter ();

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += ((advection_directions[q_point] *
				    fe_values.shape_grad(j,q_point)   *
				    (fe_values.shape_value(i,q_point) +
				     delta *
				     (advection_directions[q_point] *
				      fe_values.shape_grad(i,q_point)))) *
				   fe_values.JxW(q_point));

	    cell_rhs(i) += ((fe_values.shape_value(i,q_point) +
			     delta *
			     (advection_directions[q_point] *
			      fe_values.shape_grad(i,q_point))        ) *
			    rhs_values[i] *
			    fe_values.JxW (q_point));
	  };

      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	if (cell->face(face)->at_boundary())
	  {
	    x_fe_face_values.reinit (cell, face);
	    const FEFaceValues<dim> &fe_face_values
	      = x_fe_face_values.get_present_fe_values ();
	    
	    boundary_values.value_list (fe_face_values.get_quadrature_points(),
					face_boundary_values);
	    advection_field.value_list (fe_face_values.get_quadrature_points(),
					face_advection_directions);
	    
	    for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
	      if (fe_face_values.normal_vector(q_point) *
		  face_advection_directions[q_point]
		  < 0)
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  {
		    for (unsigned int j=0; j<dofs_per_cell; ++j)
		      cell_matrix(i,j) -= (face_advection_directions[q_point] *
					   fe_face_values.normal_vector(q_point) *
					   fe_face_values.shape_value(i,q_point) *
					   fe_face_values.shape_value(j,q_point) *
					   fe_face_values.JxW(q_point));
		    
		    cell_rhs(i) -= (face_advection_directions[q_point] *
				    fe_face_values.normal_vector(q_point) *
				    face_boundary_values[q_point]         *
				    fe_face_values.shape_value(i,q_point) *
				    fe_face_values.JxW(q_point));
		  };
	  };
      

      cell->get_dof_indices (local_dof_indices);

      assembler_lock.acquire ();
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	};
      assembler_lock.release ();
    };
}



template <int dim>
void AdvectionProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  SolverBicgstab<>        bicgstab (solver_control);

  PreconditionJacobi<> preconditioner;
  preconditioner.initialize(system_matrix, 1.0);

  bicgstab.solve (system_matrix, solution, system_rhs,
		  preconditioner);

  hanging_node_constraints.distribute (solution);
}


template <int dim>
void AdvectionProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  GradientEstimation::estimate (dof_handler,
				solution,
				estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.5, 0.03);

  triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void AdvectionProblem<dim>::output_results (const unsigned int) const
{
  GridOut grid_out;
  grid_out.write_eps (triangulation, deallog.get_file_stream());
}


template <int dim>
void AdvectionProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<3; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation, -1, 1);
	  triangulation.refine_global (2);
	}
      else
	{
	  refine_grid ();
	};
      

      deallog << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      setup_system ();

      deallog << "   Number of degrees of freedom: "
		<< dof_handler.n_dofs()
		<< std::endl;
      
      assemble_system ();
      solve ();
      output_results (cycle);
    };

  DataOut<dim,DoFHandler<dim> > data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();
  
  data_out.write_gmv (deallog.get_file_stream());
}




template <int dim>
void 
GradientEstimation::estimate (const DoFHandler<dim> &dof_handler,
			      const Vector<double>  &solution,
			      Vector<float>         &error_per_cell)
{
  Assert (error_per_cell.size() == dof_handler.get_tria().n_active_cells(),
	  ExcInvalidVectorLength (error_per_cell.size(),
				  dof_handler.get_tria().n_active_cells()));

  const unsigned int n_threads = multithread_info.n_default_threads;
  std::vector<IndexInterval> index_intervals
    = Threads::split_interval (0, dof_handler.get_tria().n_active_cells(),
			       n_threads);

  Threads::ThreadGroup<> threads;
  void (*estimate_interval_ptr) (const DoFHandler<dim> &,
				 const Vector<double> &,
				 const IndexInterval &,
				 Vector<float> &)
    = &GradientEstimation::template estimate_interval<dim>;
  for (unsigned int i=0; i<n_threads; ++i)
    threads += Threads::spawn (estimate_interval_ptr)(dof_handler, solution,
                                                      index_intervals[i],
                                                      error_per_cell);
  threads.join_all ();
}


template <int dim>
void 
GradientEstimation::estimate_interval (const DoFHandler<dim> &dof_handler,
				       const Vector<double>  &solution,
				       const IndexInterval   &index_interval,
				       Vector<float>         &error_per_cell)
{
  QMidpoint<dim> midpoint_rule;
  FEValues<dim>  x_fe_midpoint_value (dof_handler.get_fe(),
				    midpoint_rule,
				    update_values | update_q_points);
  
  Tensor<2,dim> Y;

  typename DoFHandler<dim>::active_cell_iterator cell, endc;

  cell = dof_handler.begin_active();
  advance (cell, static_cast<signed int>(index_interval.first));

  endc = dof_handler.begin_active();
  advance (endc, static_cast<signed int>(index_interval.second));

  Vector<float>::iterator
    error_on_this_cell = error_per_cell.begin() + index_interval.first;
  

  std::vector<typename DoFHandler<dim>::active_cell_iterator> active_neighbors;
  active_neighbors.reserve (GeometryInfo<dim>::faces_per_cell *
			    GeometryInfo<dim>::max_children_per_face);

  for (; cell!=endc; ++cell, ++error_on_this_cell)
    {
      x_fe_midpoint_value.reinit (cell);
      const FEValues<dim> &fe_midpoint_value = x_fe_midpoint_value.get_present_fe_values ();
      
      Y.clear ();

      Tensor<1,dim> projected_gradient;


      active_neighbors.clear ();
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	if (! cell->at_boundary(face_no))
	  {
	    const typename DoFHandler<dim>::face_iterator 
	      face = cell->face(face_no);
	    const typename DoFHandler<dim>::cell_iterator 
	      neighbor = cell->neighbor(face_no);

	    if (neighbor->active())
	      active_neighbors.push_back (neighbor);
	    else
	      {
		if (dim == 1)
		  {
		    typename DoFHandler<dim>::cell_iterator
		      neighbor_child = neighbor;
		    while (neighbor_child->has_children())
		      neighbor_child = neighbor_child->child (face_no==0 ? 1 : 0);
		    
		    Assert (neighbor_child->neighbor(face_no==0 ? 1 : 0)==cell,
			    ExcInternalError());
		    
		    active_neighbors.push_back (neighbor_child);
		  }
		else
		  for (unsigned int subface_no=0; subface_no<face->n_children(); ++subface_no)
		    active_neighbors.push_back (
		      cell->neighbor_child_on_subface(face_no, subface_no));
	      };
	  };

      const Point<dim> this_center = fe_midpoint_value.quadrature_point(0);

      std::vector<double> this_midpoint_value(1);
      fe_midpoint_value.get_function_values (solution, this_midpoint_value);
		

      std::vector<double> neighbor_midpoint_value(1);
      typename std::vector<typename DoFHandler<dim>::active_cell_iterator>::const_iterator
	neighbor_ptr = active_neighbors.begin();
      for (; neighbor_ptr!=active_neighbors.end(); ++neighbor_ptr)
	{
	  const typename DoFHandler<dim>::active_cell_iterator
	    neighbor = *neighbor_ptr;
	    
	  x_fe_midpoint_value.reinit (neighbor);
	  const FEValues<dim> &fe_midpoint_value = x_fe_midpoint_value.get_present_fe_values ();
	  const Point<dim> neighbor_center = fe_midpoint_value.quadrature_point(0);

	  fe_midpoint_value.get_function_values (solution,
                                                 neighbor_midpoint_value);

	  Point<dim>   y        = neighbor_center - this_center;
	  const double distance = std::sqrt(y.square());
	  y /= distance;
	  
	  for (unsigned int i=0; i<dim; ++i)
	    for (unsigned int j=0; j<dim; ++j)
	      Y[i][j] += y[i] * y[j];
	  
	  projected_gradient += (neighbor_midpoint_value[0] -
				 this_midpoint_value[0]) /
				distance *
				y;
	};

      AssertThrow (determinant(Y) != 0,
		   ExcInsufficientDirections());

      const Tensor<2,dim> Y_inverse = invert(Y);
      
      Point<dim> gradient;
      contract (gradient, Y_inverse, projected_gradient);
      
      *error_on_this_cell = (std::pow(cell->diameter(),
				      1+1.0*dim/2) *
			     std::sqrt(gradient.square()));
    };
}



int main () 
{
  deallog << std::setprecision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  try
    {
      deallog.depth_console (0);

      AdvectionProblem<2> advection_problem_2d;
      advection_problem_2d.run ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...) 
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };

  return 0;
}
