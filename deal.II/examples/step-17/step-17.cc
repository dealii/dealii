/* $Id$ */
/* Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2000, 2004 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


                                 // First the usual assortment of header files
                                 // we have already used in previous example
                                 // programs:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>

                                 // And here come the things that we need
                                 // particularly for this example program and
                                 // that weren't in step-8. First, we are
                                 // going to replace all linear algebra
                                 // components that involve the (global)
                                 // linear system by classes that wrap
                                 // interfaces similar to our own linear
                                 // algebra classes around what PETSc offers
                                 // (PETSc is a library written in C, and
                                 // deal.II comes with wrapper classes that
                                 // provide the PETSc functionality with an
                                 // interface that is similar to the interface
                                 // we already had for our own linear algebra
                                 // classes). In particular, we need vectors
                                 // and matrices that are distributed across
                                 // several processes in MPI programs (and
                                 // simply map to sequential, local vectors
                                 // and matrices if there is only a single
                                 // process, i.e. if you are running on only
                                 // one machine, and without MPI support):
#include <lac/petsc_vector.h>
#include <lac/petsc_parallel_vector.h>
#include <lac/petsc_parallel_sparse_matrix.h>
                                 // Then we also need interfaces for solvers
                                 // and preconditioners that PETSc provides:
#include <lac/petsc_solver.h>
#include <lac/petsc_precondition.h>
                                 // And in addition, we need some algorithms
                                 // for partitioning our meshes so that they
                                 // can be efficiently distributed across an
                                 // MPI network. The partitioning algorithm is
                                 // implemented in the ``GridTools'' class,
                                 // and we need an additional include file for
                                 // a function in ``DoFRenumbering'' that
                                 // allows to sort the indices associated with
                                 // degrees of freedom so that they are
                                 // numbered according to the subdomain they
                                 // are associated with:
#include <grid/grid_tools.h>
#include <dofs/dof_renumbering.h>

                                 // And this is simply C++ again:
#include <fstream>
#include <iostream>


                                 // Now, here comes the declaration of the
                                 // main class and of various other things
                                 // below it. As mentioned in the
                                 // introduction, almost all of this has been
                                 // copied verbatim from step-8, so we only
                                 // comment on the few things that are
                                 // different.
template <int dim>
class ElasticProblem 
{
  public:
    ElasticProblem ();
    ~ElasticProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FESystem<dim>        fe;

    ConstraintMatrix     hanging_node_constraints;

                                     // In step-8, this would have been the
                                     // place where we would have declared the
                                     // member variables for the sparsity
                                     // pattern, the system matrix, right
                                     // hand, and solution vector. We change
                                     // these declarations to use parallel
                                     // PETSc objects instead (note that the
                                     // fact that we use the parallel versions
                                     // is denoted the fact that we use the
                                     // classes from the
                                     // ``PETScWrappers::MPI'' namespace;
                                     // sequential versions of these classes
                                     // are in the ````PETScWrappers''
                                     // namespace, i.e. without the ``MPI''
                                     // part). Note also that we do not use a
                                     // separate sparsity pattern, since PETSc
                                     // manages that as part of its matrix
                                     // data structures.
    PETScWrappers::MPI::SparseMatrix system_matrix;

    PETScWrappers::MPI::Vector       solution;
    PETScWrappers::MPI::Vector       system_rhs;

                                     // The next change is that we have to
                                     // declare a variable that indicates the
                                     // MPI communicator over which we are
                                     // supposed to distribute our
                                     // computations. Note that if this is a
                                     // sequential job without support by MPI,
                                     // then PETSc provides some dummy type
                                     // for ``MPI_Comm'', so we do not have to
                                     // care here whether the job is really a
                                     // parallel one:
    MPI_Comm mpi_communicator;
    
                                     // Then we have two variables that tell
                                     // us where in the parallel world we
                                     // are. The first of the following
                                     // variables, ``n_mpi_processes'' tells
                                     // us how many MPI processes there exist
                                     // in total, while the second one,
                                     // ``this_mpi_process'', indicates which
                                     // is the number of the present process
                                     // within this space of processes. The
                                     // latter variable will have a unique
                                     // value for each process between zero
                                     // and (less than)
                                     // ``n_mpi_processes''. If this program
                                     // is run on a single machine without MPI
                                     // support, then their values are ``1''
                                     // and ``0'', respectively.
    const unsigned int n_mpi_processes;
    const unsigned int this_mpi_process;

                                     // In order to obtain values for the
                                     // above two variables, we need to query
                                     // the MPI subsystem (in case there is no
                                     // MPI running at all, these functions
                                     // automatically query some wrappers that
                                     // PETSc provides and that return default
                                     // values for a single process). We could
                                     // initialize above variables in the
                                     // constructor of this class, but since
                                     // they never change we chose to mark
                                     // them as ``const'', and so they can
                                     // only be initialized if we package all
                                     // the querying functions into auxiliary,
                                     // static functions that return the
                                     // requested values as their return
                                     // value. The argument they take denotes
                                     // the MPI communicator object from which
                                     // they shall query the total number of
                                     // processes, and the rank within this
                                     // communicator:
    static
    unsigned int
    get_n_mpi_processes (const MPI_Comm &mpi_communicator);

    static
    unsigned int
    get_this_mpi_process (const MPI_Comm &mpi_communicator);
};


                                 // The following is again taken from step-8
                                 // without change:
template <int dim>
class RightHandSide :  public Function<dim> 
{
  public:
    RightHandSide ();
    
    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &value_list) const;
};


template <int dim>
RightHandSide<dim>::RightHandSide () :
		Function<dim> (dim)
{}


template <int dim>
inline
void RightHandSide<dim>::vector_value (const Point<dim> &p,
				       Vector<double>   &values) const 
{
  Assert (values.size() == dim, 
	  ExcDimensionMismatch (values.size(), dim));
  Assert (dim >= 2, ExcInternalError());
  
  Point<dim> point_1, point_2;
  point_1(0) = 0.5;
  point_2(0) = -0.5;
  
  if (((p-point_1).square() < 0.2*0.2) ||
      ((p-point_2).square() < 0.2*0.2))
    values(0) = 1;
  else
    values(0) = 0;
  
  if (p.square() < 0.2*0.2)
    values(1) = 1;
  else
    values(1) = 0;    
}



template <int dim>
void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
					    std::vector<Vector<double> >   &value_list) const 
{
  const unsigned int n_points = points.size();

  Assert (value_list.size() == n_points, 
	  ExcDimensionMismatch (value_list.size(), n_points));

  for (unsigned int p=0; p<n_points; ++p)
    RightHandSide<dim>::vector_value (points[p],
				      value_list[p]);
}



                                 // So here first come the two functions that
                                 // query the number of processes associated
                                 // with an MPI communicator object, as well
                                 // as the rank of the present process within
                                 // it. Note again that PETSc provides dummy
                                 // implementations of these functions if no
                                 // MPI support is requested. These dummy
                                 // functions return ``1'' and ``0'' for the
                                 // total number of processes and the rank of
                                 // the present process within the
                                 // communicator, respectively.
template <int dim>
unsigned int
ElasticProblem<dim>::get_n_mpi_processes (const MPI_Comm &mpi_communicator)
{
  int n_jobs;
  MPI_Comm_size (mpi_communicator, &n_jobs);

  return n_jobs;
}



template <int dim>
unsigned int
ElasticProblem<dim>::get_this_mpi_process (const MPI_Comm &mpi_communicator)
{
  int rank;
  MPI_Comm_rank (mpi_communicator, &rank);

  return rank;
}



                                 // The first step in the actual
                                 // implementation of things is the
                                 // constructor of the main class. Apart from
                                 // initializing the same member variables
                                 // that we already had in step-8, we here
                                 // initialize the MPI communicator variable
                                 // we shall use with the global MPI
                                 // communicator linking all processes
                                 // together (in more complex applications,
                                 // one could here use a communicator object
                                 // that only links a subset of all
                                 // processes), and call above helper
                                 // functions to determine the number of
                                 // processes and where the present one fits
                                 // into this picture:
template <int dim>
ElasticProblem<dim>::ElasticProblem ()
                :
		dof_handler (triangulation),
		fe (FE_Q<dim>(1), dim),
                mpi_communicator (MPI_COMM_WORLD),
                n_mpi_processes (get_n_mpi_processes(mpi_communicator)),
                this_mpi_process (get_this_mpi_process(mpi_communicator))
{}



template <int dim>
ElasticProblem<dim>::~ElasticProblem () 
{
  dof_handler.clear ();
}


                                 // The second step is the function in which
                                 // we set up the various variables for the
                                 // global linear system to be solved.
template <int dim>
void ElasticProblem<dim>::setup_system ()
{
                                   // First, we need to generate an
                                   // enumeration for the degrees of freedom
                                   // in our problem. Further below, we will
                                   // show how we assign each cell to one of
                                   // the MPI processes before we even get
                                   // here. What we then need to do is to
                                   // enumerate the degrees of freedom in a
                                   // way so that all degrees of freedom
                                   // associated with cells in subdomain zero
                                   // (which resides on process zero) come
                                   // before all DoFs associated with cells on
                                   // subdomain one, before those on cells on
                                   // process two, and so on. We need this
                                   // since we have to split the global
                                   // vectors for right hand side and
                                   // solution, as well as the matrix into
                                   // contiguous chunks of rows that live on
                                   // each of the processors, and we will want
                                   // to do this in a way that requires
                                   // minimal communication. This is done
                                   // using the following two functions, which
                                   // first generates an initial ordering of
                                   // all degrees of freedom, and then re-sort
                                   // them according to above criterion:
  dof_handler.distribute_dofs (fe);
  DoFRenumbering::subdomain_wise (dof_handler);

                                   // While we're at it, let us also count how
                                   // many degrees of freedom there exist on
                                   // the present process:
  const unsigned int n_local_dofs
    = DoFTools::count_dofs_with_subdomain_association (dof_handler,
                                                       this_mpi_process);  

                                   // Then we initialize the system matrix,
                                   // solution, and right hand side
                                   // vectors. Since they all need to work in
                                   // parallel, we have to pass them an MPI
                                   // communication object, as well as their
                                   // global sizes, and also how many rows out
                                   // of this global size are to be stored
                                   // locally:
  system_matrix.reinit (mpi_communicator,
                        dof_handler.n_dofs(),
                        dof_handler.n_dofs(),
                        n_local_dofs,
                        dof_handler.max_couplings_between_dofs());

  solution.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
  system_rhs.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

                                   // Finally, we need to initialize the
                                   // objects denoting hanging node
                                   // constraints for the present grid. Note
                                   // that since PETSc handles the sparsity
                                   // pattern internally to the matrix, there
                                   // is no need to set up an independent
                                   // sparsity pattern here, and to condense
                                   // it for constraints, as we have done in
                                   // all other example programs.
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  hanging_node_constraints.close ();
}


template <int dim>
void ElasticProblem<dim>::assemble_system () 
{
                                   // xxx move to front
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(dim),
					    boundary_values);

  QGauss2<dim>  quadrature_formula;
  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  std::vector<double>     lambda_values (n_q_points);
  std::vector<double>     mu_values (n_q_points);

  ConstantFunction<dim> lambda(1.), mu(1.);

  RightHandSide<dim>      right_hand_side;
  std::vector<Vector<double> > rhs_values (n_q_points,
					   Vector<double>(dim));


  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();
  for (; cell!=endc; ++cell)
                                     // xxx
    if (cell->subdomain_id() == this_mpi_process)
      {
        cell_matrix.clear ();
        cell_rhs.clear ();

        fe_values.reinit (cell);
      
        lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
        mu.value_list     (fe_values.get_quadrature_points(), mu_values);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int 
              component_i = fe.system_to_component_index(i).first;
	  
            for (unsigned int j=0; j<dofs_per_cell; ++j) 
              {
                const unsigned int 
                  component_j = fe.system_to_component_index(j).first;
	      
                for (unsigned int q_point=0; q_point<n_q_points;
                     ++q_point)
                  {
                    cell_matrix(i,j) 
                      += 
                      (
                        (fe_values.shape_grad(i,q_point)[component_i] *
                         fe_values.shape_grad(j,q_point)[component_j] *
                         lambda_values[q_point])
                        +
                        (fe_values.shape_grad(i,q_point)[component_j] *
                         fe_values.shape_grad(j,q_point)[component_i] *
                         mu_values[q_point])
                        +
                        ((component_i == component_j) ?
                         (fe_values.shape_grad(i,q_point) *
                          fe_values.shape_grad(j,q_point) *
                          mu_values[q_point])  :
                         0)
                        )
                      *
                      fe_values.JxW(q_point);
                  };
              };
          };

        right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
                                           rhs_values);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const unsigned int 
              component_i = fe.system_to_component_index(i).first;
	  
            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_rhs(i) += fe_values.shape_value(i,q_point) *
                             rhs_values[q_point](component_i) *
                             fe_values.JxW(q_point);
          };

        cell->get_dof_indices (local_dof_indices);

                                         //xxx
        MatrixTools::local_apply_boundary_values (boundary_values,
                                                  local_dof_indices,
                                                  cell_matrix,
                                                  cell_rhs,
                                                  false);

                                         // xxx
        hanging_node_constraints
          .distribute_local_to_global (cell_matrix,
                                       local_dof_indices,
                                       system_matrix);

        hanging_node_constraints
          .distribute_local_to_global (cell_rhs,
                                       local_dof_indices,
                                       system_rhs);
      }

                                   //xxx no condense necessary, no apply_b_v
                                   //either


                                   // xxx
  system_matrix.compress ();
  system_rhs.compress ();
}



template <int dim>
void ElasticProblem<dim>::solve () 
{
                                   // xxx
  SolverControl           solver_control (1000, 1e-10);
  PETScWrappers::SolverCG cg (solver_control,
                              mpi_communicator);

  PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);
  
  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  PETScWrappers::Vector localized_solution (solution);
  hanging_node_constraints.distribute (localized_solution);


  std::vector<unsigned int> subdomain_association (dof_handler.n_dofs());
  DoFTools::get_subdomain_association (dof_handler,
                                       subdomain_association);  
  for (unsigned int i=0; i<localized_solution.size(); ++i)
    if (subdomain_association[i] == this_mpi_process)
      solution(i) = static_cast<PetscScalar>(localized_solution(i));
  solution.compress ();  

  if (this_mpi_process == 0)
    std::cout << "   Solver converged in "
              << solver_control.last_step()
              << " iterations." << std::endl;
}



template <int dim>
void ElasticProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

                                   // XXX
  {
    PETScWrappers::Vector localized_solution (solution);
    Vector<float> local_error_per_cell (triangulation.n_active_cells());
    
    typename FunctionMap<dim>::type neumann_boundary;
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss2<dim-1>(),
                                        neumann_boundary,
                                        localized_solution,
                                        local_error_per_cell,
                                        std::vector<bool>(),
                                        0,
                                        multithread_info.n_default_threads,
                                        this_mpi_process);
    
    const unsigned int local_cells
      = (n_mpi_processes == 1 ?
         triangulation.n_active_cells() :
         (this_mpi_process != n_mpi_processes-1 ?
          triangulation.n_active_cells() / n_mpi_processes :
          triangulation.n_active_cells() - triangulation.n_active_cells() / n_mpi_processes * (n_mpi_processes-1)));
    PETScWrappers::MPI::Vector
      global_error_per_cell (mpi_communicator,
                             triangulation.n_active_cells(),
                             local_cells);


    for (unsigned int i=0; i<local_error_per_cell.size(); ++i)
      if (local_error_per_cell(i) != 0)
        global_error_per_cell(i) = local_error_per_cell(i);
    global_error_per_cell.compress ();

    estimated_error_per_cell = global_error_per_cell;
  }

  
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();

                                   // xxx
  GridTools::partition_triangulation (n_mpi_processes, triangulation);
}


template <int dim>
void ElasticProblem<dim>::output_results (const unsigned int cycle) const
{
                                   // xxx
  PETScWrappers::Vector global_solution;
  global_solution = solution;
      
  if (this_mpi_process == 0)
    {
      std::string filename = "solution-";
      filename += ('0' + cycle);
      Assert (cycle < 10, ExcInternalError());
  
      filename += ".gmv";
      std::ofstream output (filename.c_str());

      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);

      std::vector<std::string> solution_names;
      switch (dim)
        {
          case 1:
                solution_names.push_back ("displacement");
                break;
          case 2:
                solution_names.push_back ("x_displacement");	    
                solution_names.push_back ("y_displacement");
                break;
          case 3:
                solution_names.push_back ("x_displacement");	    
                solution_names.push_back ("y_displacement");
                solution_names.push_back ("z_displacement");
                break;
          default:
                Assert (false, ExcInternalError());
        };

                                       // xxx
      std::vector<unsigned int> p (triangulation.n_active_cells());
      GridTools::get_subdomain_association (triangulation, p);
      Vector<double> x(p.begin(), p.end());
      
      data_out.add_data_vector (x, "partitioning");
      data_out.add_data_vector (global_solution, solution_names);
      data_out.build_patches ();
      data_out.write_gmv (output);
    }
}



template <int dim>
void ElasticProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<10; ++cycle)
    {
                                       // xxx
      if (this_mpi_process == 0)
        std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation, -1, 1);
	  triangulation.refine_global (3);

                                           // xxx
          GridTools::partition_triangulation (n_mpi_processes, triangulation);
	}
      else
	refine_grid ();

                                       // xxx
      if (this_mpi_process == 0)
        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl;

      setup_system ();

                                       // xxx
      if (this_mpi_process == 0)
        {
          std::cout << "   Number of degrees of freedom: "
                    << dof_handler.n_dofs()
                    << " (by partition:";
          for (unsigned int partition=0; partition<n_mpi_processes; ++partition)
            std::cout << (partition==0 ? ' ' : '+')
                      << (DoFTools::
                          count_dofs_with_subdomain_association (dof_handler,
                                                                 partition));
          std::cout << ")" << std::endl;
        }
      
      assemble_system ();
      solve ();
      output_results (cycle);
    };
}


int main (int argc, char **argv) 
{
  try
    {
                                       // xxx
      PetscInitialize(&argc,&argv,0,0);
      deallog.depth_console (0);

                                       // xxx localize scope
      {
        ElasticProblem<2> elastic_problem_2d;
        elastic_problem_2d.run ();
      }

                                       // xxx
      PetscFinalize();      
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
