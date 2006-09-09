/*    $Id: step-4.cc,v 1.34 2006/02/06 21:33:10 wolf Exp $       */
/*    Version: $Name:  $                                          */
/*                                                                */
/*    Copyright (C) 2006 by the deal.II authors */
/*    Author: Xing Jin, Wolfgang Bangerth, Texas A&M University, 2006 */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */
                                                          

				 // @sect3{Include files}

				 // The following have all been covered
				 // previously:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/utilities.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>

#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <numerics/data_out.h>
#include <numerics/matrices.h>
#include <numerics/vectors.h>

#include <fstream>
#include <iostream>
#include <sstream>


				 // @sect3{The "forward problem" class template}

				 // The main class is similar to the wave
				 // equation.  The difference is that we add
				 // an absorbing boundary condition. Because
				 // we are only interested in values at
				 // specific locations, we define some
				 // parameters to obtain the coordinates of
				 // those locations.
template <int dim>
class TATForwardProblem
{
  public:
    TATForwardProblem ();
    void run ();
    
  private:
    void setup_system ();
    void solve_p ();
    void solve_v ();
    void output_results () const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix constraints;
    
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;

    Vector<double>       solution_p, solution_v;
    Vector<double>       old_solution_p, old_solution_v;
    Vector<double>       system_rhs_p, system_rhs_v;

    double time, time_step;
    unsigned int timestep_number;
    const double theta;

				     //
    SparseMatrix<double> boundary_matrix;
				     // Number of refinement
    const unsigned int n_refinements;
				     // The acoustic speed in the medium $c_0$
    const double acoustic_speed;

				     // The detector circullarly scan the target region.
				     // The step size of the detector is in angles
    const double step_angle;
				     // The scanning radius
    const double radius;

    const double end_time;
};

 
				 // Declare a class template for the right hand side
				 // of the pressure potential 
template <int dim>
class RightHandSideP : public Function<dim> 
{
  public:
    RightHandSideP () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};

				 // Declare a class template for the right hand side
				 // of the derivative of the pressure potential                  
template <int dim>
class RightHandSideV : public Function<dim> 
{
  public:
    RightHandSideV () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};

				 // Declare a class template for the initial values
				 // of the pressure potential
template <int dim>
class InitialValuesP : public Function<dim> 
{
  public:
    InitialValuesP () : Function<dim>() {};
    
    virtual double value (const Point<dim> &p,
			  const unsigned int  component = 0) const;
};

				 // Declare a class template for the initial values
				 // of the derivative of the pressure potential
template <int dim>
class InitialValuesV : public Function<dim> 
{
  public:
    InitialValuesV () : Function<dim>() {};
    
    virtual double value (const Point<dim> &p, 
			  const unsigned int  component = 0) const;
};

				 // Here is the function to set the right hand side
				 // values to be zero for pressure potential
template <int dim>
double RightHandSideP<dim>::value (const Point<dim> &/*p*/,
				    const unsigned int /*component*/) const 
{
  return 0;
}
				 // Similarly we set the right-hand size of the 
				 // derivative of the pressure potential to be 
				 // zero
template <int dim>
double RightHandSideV<dim>::value (const Point<dim> &/*p*/,
				    const unsigned int /*component*/) const 
{
  return 0;
}


				 // The sources of the thermoacoustic waves 
				 // are small absorbers. We will compare the 
				 // simulation results with the experimental
				 // data.

template <int dim>
double InitialValuesP<dim>::value (const Point<dim> &p,
				    const unsigned int /*component*/) const       
{
                          
  if (std::sqrt(p.square())< 0.025 )
    return 1;
				   // The "distance" function is used to compute
				   // the Euclidian distance between two points.
                                    
  if (p.distance(Point<dim>(-0.135,0))<0.05)
    return 1;
                           
  if (p.distance(Point<dim>(0.17,0))<0.03)
    return 1;

  if (p.distance(Point<dim>(-0.25,0))<0.02)
    return 1;

  if (p.distance(Point<dim>(-0.05,-0.15))<0.015)
    return 1;

  return 0;
}
				 // Initial value for the derivative of
				 // pressure potential is set to zero
template <int dim>   
double InitialValuesV<dim>::value (const Point<dim> &/*p*/,
				    const unsigned int /*component*/) const 
{
 
  return 0;
}


				 // @sect4{Initialize the problem}
				 // Acoustic_speed here is the acoustic speed 
				 // in the medium. Specifically we use acoustic speed
				 // in mineral oil. We use Crank-Nicolson scheme
				 // for our time-dependent problem, therefore theta is
				 // set to be 0.5.  The step size of the detector
				 // is 2.25 degree, which means we need 160 steps 
				 // in order to finish a circular scan. The radius of the
				 // scanning circle is select to be half way between  
				 // the center and the boundary to avoid the reflections 
				 // from the the boundary, so as to miminize the 
				 // interference brought by the inperfect absorbing 
				 // boundary condition. The time step is selected 
				 // to satisfy $k = h/c$                           
template <int dim>
TATForwardProblem<dim>::TATForwardProblem () :
                fe (1),
		dof_handler (triangulation),
		n_refinements (7),
                acoustic_speed (1.437),
                theta (0.5),
	        end_time (0.7),
 		time_step (0.5/std::pow(2.,1.0*n_refinements)/acoustic_speed),  
                step_angle (2.25),  
                radius (0.5)

{}   



				 // @sect4{TATForwardProblem::setup_system}

				 // The following system is pretty much what
				 // we've already done in @ref step_23
				 // "step-23", but with two important
				 // differences. First, we have to create a
				 // circular (or spherical) mesh around the
				 // origin, with a radius of 1. This nothing
				 // new: we've done so before in @ref step_6
				 // "step-6", @ref step_10 "step-10", and @ref
				 // step_11 "step-11", where we also explain
				 // how to attach a boundary object to a
				 // triangulation to be used whenever the
				 // triangulation needs to know where new
				 // boundary points lie when a cell is
				 // refined. Following this, the mesh is
				 // refined <code>n_refinements</code> times
				 // &mdash; this variable was introduced to
				 // make sure the time step size is always
				 // compatible with the cell size, and
				 // therefore satisfies the CFL condition that
				 // was talked about in the introduction of
				 // @ref step_23 "step-23".
				 //
				 // The only other significant change is that
				 // we need to build the boundary mass
				 // matrix. We will comment on this further
				 // down below.
template <int dim>
void TATForwardProblem<dim>::setup_system ()
{
  GridGenerator::hyper_ball (triangulation, Point<dim>(), 1.);
  static const HyperBallBoundary<dim> boundary_description(center);
  triangulation.set_boundary (0,boundary_description);
  triangulation.refine_global (n_refinements);

  std::cout << "Number of active cells: "
	    << triangulation.n_active_cells()
  	    << std::endl;

  dof_handler.distribute_dofs (fe);

  std::cout << "Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl
	    << std::endl;

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
  mass_matrix.reinit (sparsity_pattern);
  laplace_matrix.reinit (sparsity_pattern);

  MatrixCreator::create_mass_matrix (dof_handler, QGauss<dim>(3),
				     mass_matrix);
  MatrixCreator::create_laplace_matrix (dof_handler, QGauss<dim>(3),
					laplace_matrix);

				   // The second difference, as mentioned, to
				   // @ref step_23 "step-23" is that we need
				   // to build the boundary mass matrix that
				   // grew out of the absorbing boundary
				   // conditions.
				   //
				   // A first observation would be that this
				   // matrix is much sparser than the regular
				   // mass matrix, since none of the shape
				   // functions with purely interior support
				   // contributes to this matrix. We could
				   // therefore optimize the storage pattern
				   // to this situation and build up a second
				   // sparsity pattern that only contains the
				   // nonzero entries that we need. There is a
				   // trade-off to make here: first, we would
				   // have to have a second sparsity pattern
				   // object, so that costs memory. Secondly,
				   // the matrix attached to this sparsity
				   // pattern is going to be smaller and
				   // therefore requires less memore; it would
				   // also be faster to perform matrix-vector
				   // multiplications with it. The final
				   // argument, however, is the one that tips
				   // the scale: we are not primarily
				   // interested in performing matrix-vector
				   // with the boundary matrix alone (though
				   // we need to do that for the right hand
				   // side vector once per time step), but
				   // mostly wish to add it up to the other
				   // matrices used in the first of the two
				   // equations since this is the one that is
				   // going to be multiplied with once per
				   // iteration of the CG method,
				   // i.e. significantly more often. It is now
				   // the case that the SparseMatrix::add
				   // class allows to add one matrix to
				   // another, but only if they use the same
				   // sparsity pattern (the reason being that
				   // we can't add nonzero entries to a matrix
				   // after the sparsity pattern has been
				   // created, so we simply require that the
				   // two matrices have the same sparsity
				   // pattern.
				   //
				   // So let's go with that:
  boundary_matrix.reinit (sparsity_pattern);

				   // The second thing to do is to actually
				   // build the matrix. Here, we need to
				   // integrate over faces of cells, so first
				   // we need a quadrature object that works
				   // on <code>dim-1</code> dimensional
				   // objects. Secondly, the FEFaceValues
				   // variant of FEValues that works on faces,
				   // as its name suggest. And finally, the
				   // other variables that are part of the
				   // assembly machinery. All of this we put
				   // between curly braces to limit the scope
				   // of these variables to where we actually
				   // need them.
				   //
				   // The actual act of assembling the matrix
				   // is then fairly straightforward: we loop
				   // over all cells, over all faces of each
				   // of these cells, and then do something
				   // only if that particular face is at the
				   // boundary of the domain. Like this:
  {
    const QGauss<dim-1>  quadrature_formula(3);
    FEFaceValues<dim> fe_values (fe, quadrature_formula, 
				 update_values  |  update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

        
                              
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	if (cell->at_boundary(f))
	  {
	    cell_matrix = 0;

	    fe_values.reinit (cell, f);

	    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  cell_matrix(i,j) += (fe_values.shape_value(i,q_point) *
				       fe_values.shape_value(j,q_point) *
				       fe_values.JxW(q_point));

	    cell->get_dof_indices (local_dof_indices);
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		boundary_matrix.add (local_dof_indices[i],
				     local_dof_indices[j],
				     cell_matrix(i,j));
	  }
  
  }

  system_matrix.copy_from (mass_matrix);
  system_matrix.add (time_step * time_step * theta * theta *
		     acoustic_speed * acoustic_speed,
		     laplace_matrix);
  system_matrix.add (acoustic_speed * theta * time_step, boundary_matrix);
  

  solution_p.reinit (dof_handler.n_dofs());
  old_solution_p.reinit (dof_handler.n_dofs());
  system_rhs_p.reinit (dof_handler.n_dofs());

  solution_v.reinit (dof_handler.n_dofs());
  old_solution_v.reinit (dof_handler.n_dofs());
  system_rhs_v.reinit (dof_handler.n_dofs());

  constraints.close ();
}


				 // @sect4{TATForwardProblem::solve_p and TATForwardProblem::solve_v}

				 // The following two functions, solving the
				 // linear systems for the pressure and the
				 // velocity variable, are taken pretty much
				 // verbatim (with the exception of the change
				 // of name from $u$ to $p$ of the primary
				 // variable) from @ref step_23 "step-23":
template <int dim>
void TATForwardProblem<dim>::solve_p () 
{
  SolverControl           solver_control (1000, 1e-8*system_rhs_p.l2_norm());
  SolverCG<>              cg (solver_control);

  cg.solve (system_matrix, solution_p, system_rhs_p,
	    PreconditionIdentity());

  std::cout << "   p-equation: " << solver_control.last_step()
	    << " CG iterations."
	    << std::endl;
}


template <int dim>
void TATForwardProblem<dim>::solve_v () 
{
  SolverControl           solver_control (1000, 1e-8*system_rhs_v.l2_norm());
  SolverCG<>              cg (solver_control);

  cg.solve (mass_matrix, solution_v, system_rhs_v,
	    PreconditionIdentity());

  std::cout << "   v-equation: " << solver_control.last_step()
	    << " CG iterations."
	    << std::endl;
}



				 // @sect4{TATForwardProblem::output_results}

				 // The same holds here: the function is from
				 // @ref step_23 "step-23".
template <int dim>
void TATForwardProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution_p, "P");
  data_out.add_data_vector (solution_v, "V");

  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-"
	   << Utilities::int_to_string (timestep_number, 3)
	   << ".gnuplot";
  std::ofstream output (filename.str().c_str());
  data_out.write_gnuplot (output);
}


//XXX
                                 // This is the main function
                                 // 
                                      
template <int dim>
void TATForwardProblem<dim>::run () 
{
  setup_system();

  VectorTools::project (dof_handler,constraints, 
			QGauss<dim>(3), InitialValuesP<dim>(),
			old_solution_p);
  VectorTools::project (dof_handler,constraints, 
			QGauss<dim>(3), InitialValuesV<dim>(),
			old_solution_v);


  timestep_number = 1;
  unsigned int n_steps;
  unsigned int n_detectors;
  double scanning_angle;

				   // Number of time steps is defined as the
				   // ratio of the total time to the time step                                 
  n_steps=static_cast<unsigned int>(std::floor(end_time/time_step));     
				   // Number of detector positions is defined          
				   // as the ratio of 360 degrees to the step
				   // angle
  n_detectors=static_cast<unsigned int>(std::ceil(360/step_angle));
				   // Define two vectors to hold the coordinates
				   // of the detectors in the scanning
				   // geometry
  Vector<double> detector_x (n_detectors+1);
  Vector<double> detector_y (n_detectors+1);
				   // Define a vector to hold the value obtained
				   // by the detector
  Vector<double> project_dat (n_steps * n_detectors +1);
				   // Get the coordinates of the detector on the 
				   // different locations of the circle.
				   // Scanning angle is viewing angle at 
				   // current position. The coordinates of
				   // the detectors are calculated from the radius
				   // and scanning angle.
  scanning_angle=0;
  for (unsigned int i=1; i<=n_detectors;  i++){
				     // Scanning clockwisely. We need to change the angles
				     // into radians because std::cos and std:sin accept
				     // values in radian only
    scanning_angle -= step_angle/180 * 3.14159265;   
    detector_x(i) = radius * std::cos(scanning_angle);
    detector_y(i) = radius * std::sin(scanning_angle);
  }

  std::cout<< "Total number of time steps = "<< n_steps <<std::endl;
  std::cout<< "Total number of detectors = "<< n_detectors << std::endl;

				   // Open a file to write the data
				   // obtained by the detectors 
                         
  std::ofstream proj_out;
  proj_out.open("proj.dat");
  

  for (double time = time_step; time<=end_time; time+=time_step, ++timestep_number)
    {
      std::cout << std::endl;                                       
      std::cout<< "time_step " << timestep_number << " @ t=" << time << std::endl;

      Vector<double> tmp1 (solution_p.size());
      Vector<double> tmp2 (solution_v.size());
      Vector<double> F1 (solution_p.size());
      Vector<double> F2 (solution_v.size());
				       // Calculate G1 as defined in the introduction section
                         
      mass_matrix.vmult (tmp1, old_solution_p); 
      mass_matrix.vmult (tmp2, old_solution_v); 
      F1 = tmp1;
      F1.add(time_step * (1-theta), tmp2); 
				       // Calculate G2 as defined in the introduction section
      mass_matrix.vmult (tmp1, old_solution_v);
      laplace_matrix.vmult (tmp2, old_solution_p); 
      F2 = tmp1;
      F2.add(-acoustic_speed*acoustic_speed*time_step*(1-theta), tmp2);
      tmp1=0;
      boundary_matrix.vmult (tmp1,old_solution_p);
      F2.add(acoustic_speed,tmp1);
      
				       // Compute the pressure potential p, the formula
				       // has been presented in the introduction section

      system_rhs_p = F1; 
      system_rhs_p.add(time_step * theta , F2);

      RightHandSideP<dim> rhs_function_p;
      rhs_function_p.set_time (time);

      tmp1=0;
      VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(2),
                                           rhs_function_p, tmp1);
    
      system_rhs_p.add(-theta * theta * time_step * time_step*acoustic_speed*acoustic_speed,tmp1); 
      rhs_function_p.set_time (time-time_step);
      tmp1=0;
      VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(2),
                                           rhs_function_p, tmp1);
      
 
      system_rhs_p.add(-theta * (1-theta) * time_step * time_step*acoustic_speed*acoustic_speed,tmp1); 

      solve_p ();

				       // Compute the derivative potential pressure.
				       // The formula has been presented in the introduction
				       // section. The potential derivative is calculated 
				       // after the potential pressure because the calculation
				       // depends on the current value of the potential 
				       // pressure

      system_rhs_v = F2;
      tmp1 = 0;
      laplace_matrix.vmult (tmp1, solution_p);
      system_rhs_v.add(-time_step * theta*acoustic_speed*acoustic_speed, tmp1);
      tmp1 = 0;
      boundary_matrix.vmult(tmp1, solution_p);
      system_rhs_v.add(-acoustic_speed,tmp1);
      
      RightHandSideV<dim> rhs_function_v;
      rhs_function_v.set_time (time); 

      tmp2 = 0;
      VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(2),
                                           rhs_function_v, tmp2);

      system_rhs_p.add(-theta * time_step*acoustic_speed*acoustic_speed,tmp2); 

      rhs_function_v.set_time (time-time_step);
      tmp2 = 0;
      VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(2),
                                           rhs_function_v, tmp2);
      system_rhs_p.add(-(1-theta)*time_step*acoustic_speed*acoustic_speed,tmp2);
      
      solve_v ();
				       // Compute the energy in the system.By checking
				       // energy change in the system, we can verify
				       // the correctness of the code. 

      double energy = (mass_matrix.matrix_scalar_product(solution_v,solution_v)+
		       acoustic_speed*acoustic_speed*laplace_matrix.matrix_scalar_product(solution_p,solution_p))/2;        
                                                                         
      std::cout << "energy= " << energy << std::endl;

      output_results ();
      
				       //  Evaluate the value at specific locations. 
				       //  For 2-D, it is on a circle. For 1-D, 
				       //  it is a point detector.

      proj_out << time ;

      for (unsigned i=1 ; i<=n_detectors; i++){
        project_dat((timestep_number-1)*n_detectors+i)
	  = VectorTools::point_value (dof_handler,solution_p, 
				      Point<2>(detector_x(i),detector_y(i)));
        proj_out << " "<< project_dat((timestep_number-1)*n_detectors+i)<<" " ;
      }

      proj_out<<std::endl;
          
				       // Update the values for the pressure potential 
				       // and its derivative. 
         
      old_solution_p = solution_p;
      solution_p = 0;
      old_solution_v = solution_v;      
      solution_v = 0;
    }
}



				 // @sect3{The <code>main</code> function}

				 // What remains is the main function of the
				 // program. There is nothing here that hasn't
				 // been shown in several of the previous
				 // programs:
int main () 
{
  try
    {
      deallog.depth_console (0);
      TATForwardProblem<2> forward_problem_solver;
      forward_problem_solver.run ();
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
