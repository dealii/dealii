/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */

				 // As usual, the program starts with
				 // a rather long list of include
				 // files which you are probably
				 // already used to by now:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/mapping_q.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

				 // Just this one is new: it declares
				 // a class
				 // ``CompressedSparsityPattern'',
				 // which we will use and explain
				 // further down below.
#include <lac/compressed_sparsity_pattern.h>

				 // We will make use of the std::find
				 // algorithm of the C++ standard
				 // library, so we have to include the
				 // following file for its
				 // declaration:
#include <algorithm>


template <int dim>
double measure (const DoFHandler<dim> &dof_handler,
		const Mapping<dim>    &mapping)
{
  QGauss4<dim> quadrature;
  FEValues<dim> fe_values (mapping, dof_handler.get_fe(), quadrature,
			   update_JxW_values);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  double measure = 0;
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
	measure += fe_values.JxW (i);
    };
  return measure;
};


template <int dim>
double measure (const Triangulation<dim> &triangulation,
		const Mapping<dim>       &mapping)
{
  FE_Q<dim> dummy_fe(1);
  DoFHandler<dim> dof_handler (const_cast<Triangulation<dim>&>(triangulation));
  dof_handler.distribute_dofs(dummy_fe);
  return measure (dof_handler, mapping);
};


template <int dim>
double surface (const DoFHandler<dim> &dof_handler,
		const Mapping<dim>    &mapping)
{
  QGauss4<dim-1> quadrature;
  FEFaceValues<dim> fe_values (mapping, dof_handler.get_fe(), quadrature,
			       update_JxW_values);
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  double surface = 0;
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->at_boundary())
	{
	  fe_values.reinit (cell, face);
	  for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
	    surface += fe_values.JxW (i);
	};
  return surface;
};


template <int dim>
double surface (const Triangulation<dim> &triangulation,
		const Mapping<dim>       &mapping)
{
  FE_Q<dim> dummy_fe(1);
  DoFHandler<dim> dof_handler (const_cast<Triangulation<dim>&>(triangulation));
  dof_handler.distribute_dofs(dummy_fe);
  return surface (dof_handler, mapping);
};


template double surface (const Triangulation<2> &, const Mapping<2> &);
template double measure (const Triangulation<2> &, const Mapping<2> &);




				 // Then we declare a class which
				 // represents the solution of a
				 // Laplace problem. As this example
				 // program is based on step-5, the
				 // class looks rather the same, with
				 // the sole structural difference
				 // that we have merged the functions
				 // ``assemble_system'' and ``solve'',
				 // and the output function was
				 // dropped since the solution
				 // function is so boring that it is
				 // not worth being viewed.
				 //
				 // The only other noteworthy change
				 // is that the constructor takes a
				 // value representing the polynomial
				 // degree of the mapping to be used
				 // later on, and that it has another
				 // member variable representing
				 // exactly this mapping. In general,
				 // this variable will occur in real
				 // applications at the same places
				 // where the finite element is
				 // declared or used.
template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem (const unsigned int mapping_degree);
    void run ();
    
  private:
    void setup_system ();
    void assemble_and_solve ();
    void solve ();

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    MappingQ<dim>        mapping;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    ConstraintMatrix     mean_value_constraints;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};



				 // Construct such an object, by
				 // initializing the variables. Here,
				 // we use linear finite elements (the
				 // argument to the ``fe'' variable
				 // denotes the polynomial degree),
				 // and mappings of given order. Print
				 // to screen what we are about to do.
template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int mapping_degree) :
                fe (1),
		dof_handler (triangulation),
		mapping (mapping_degree)
{
  std::cout << "Using mapping with degree " << mapping_degree << ":"
	    << std::endl
	    << "============================"
	    << std::endl;
};



				 // The first task is to set up the
				 // variables for this problem. This
				 // includes generating a valid
				 // ``DoFHandler'' object, as well as
				 // the sparsity patterns for the
				 // matrix, and the object
				 // representing the constraints that
				 // the mean value of the degrees of
				 // freedom on the boundary be zero.
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
				   // The first task is trivial:
				   // generate an enumeration of the
				   // degrees of freedom:
  dof_handler.distribute_dofs (fe);

				   // Next task is to construct the
				   // object representing the
				   // constraint that the mean value
				   // of the degrees of freedom on the
				   // boundary shall be zero. For
				   // this, we first want a list of
				   // those nodes which are actually
				   // at the boundary. The
				   // ``DoFTools'' class has a
				   // function that returns an array
				   // of boolean values where ``true''
				   // indicates that the node is at
				   // the boundary. The second
				   // argument denotes a mask
				   // selecting which components of
				   // vector valued finite elements we
				   // want to be considered. Since we
				   // have a scalar finite element
				   // anyway, this mask consists of
				   // only one entry, and its value
				   // must be ``true''.
  std::vector<bool> boundary_dofs (dof_handler.n_dofs(), false);
  DoFTools::extract_boundary_dofs (dof_handler, std::vector<bool>(1,true),
				   boundary_dofs);
  
				   // Let us first pick out the first
				   // boundary node from this list. We
				   // do that by searching for the
				   // first ``true'' value in the
				   // array (note that ``std::find''
				   // returns an iterator to this
				   // element), and computing its
				   // distance to the overall first
				   // element in the array to get its
				   // index:
  const unsigned int first_boundary_dof
    = std::distance (std::find (boundary_dofs.begin(),
				boundary_dofs.end(),
				true),
		     boundary_dofs.begin());
	
  mean_value_constraints.clear ();
  mean_value_constraints.add_line (first_boundary_dof);
  for (unsigned int i=first_boundary_dof+1; i<dof_handler.n_dofs(); ++i)
    if (boundary_dofs[i] == true)
      mean_value_constraints.add_entry (first_boundary_dof,
					i, -1);
  mean_value_constraints.close ();

  CompressedSparsityPattern csp (dof_handler.n_dofs(),
				 dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp);
  mean_value_constraints.condense (csp);

  sparsity_pattern.copy_from (csp);

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};



template <int dim>
void LaplaceProblem<dim>::assemble_and_solve () 
{  
  QGauss2<dim>  cell_quadrature;
  QGauss2<dim-1> face_quadrature;
  MatrixTools::create_laplace_matrix (mapping, dof_handler,
				      cell_quadrature,
				      system_matrix);
  VectorTools::create_right_hand_side (mapping, dof_handler,
				       cell_quadrature,
				       ConstantFunction<dim>(-2),
				       system_rhs);
  
  Vector<double> tmp (system_rhs.size());
  VectorTools::create_boundary_right_hand_side (mapping, dof_handler,
						face_quadrature,
						ConstantFunction<dim>(1),
						tmp);
  system_rhs += tmp;

  mean_value_constraints.condense (system_matrix);
  mean_value_constraints.condense (system_rhs);  

  solve ();
  mean_value_constraints.distribute (solution);
  
  Vector<float> difference_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (mapping, dof_handler,
				     solution,
				     ZeroFunction<dim>(),
				     difference_per_cell,
				     QGauss3<dim>(),
				     H1_seminorm);
  std::cout << "  " << triangulation.n_active_cells() << " cells:  "
	    << "  |u|_1="
	    << difference_per_cell.l2_norm()
	    << ", error="
	    << fabs(difference_per_cell.l2_norm()-sqrt(3.14159265358/2))
	    << std::endl;
};



template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
};



template <int dim>
void LaplaceProblem<dim>::run () 
{
  GridGenerator::hyper_ball (triangulation);
  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);
  
  for (unsigned int cycle=0; cycle<6; ++cycle, triangulation.refine_global(1))
    {
      setup_system ();
      assemble_and_solve ();
    };
};

    

				 // Finally the main function. It's
				 // structure is the same as that used
				 // in several of the previous
				 // examples, so probably needs no
				 // more explanation.
int main () 
{
  try
    {
      deallog.depth_console (0);
      std::cout.precision(5);

				       // This is the main loop, doing
				       // the computations with
				       // mappings of linear through
				       // cubic mappings. Note that
				       // since we need the object of
				       // type ``LaplaceProblem<2>''
				       // only once, we do not even
				       // name it, but create an
				       // unnamed such object and call
				       // the ``run'' function of it,
				       // subsequent to which it is
				       // immediately destroyed again.
      for (unsigned int mapping_degree=1; mapping_degree<=3; ++mapping_degree)
	LaplaceProblem<2>(mapping_degree).run ();
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
