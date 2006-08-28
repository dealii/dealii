
/*    $Id: project.cc modified from heat-equation.cc 2006/03/05 $ */
/*    Author: Xing Jin                                            */
/*                                                                */
/*    $Id: step-4.cc,v 1.34 2006/02/06 21:33:10 wolf Exp $        */
/*    Version: $Name:  $                                          */
/*                                                                */
/*    Copyright (c) 1999,2000,2001,2002,2003,2004,2005,2006       */
/*    by the deal.II authors.                                     */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */
                                                          

                           // @sect3{Include files}
                           // Most include files have been covered in 
                           // step-6 and will not be further commented on
#include <grid/tria.h>
#include <dofs/dof_handler.h>
                           // We will need to read the value at a specific 
                           // location. This including file is needed for 
                           // finding a cell that contains a given point
#include <grid/grid_tools.h>
#include <grid/grid_generator.h>
                           // Because the scanning geometry is on a circle,
                           // the boundaries are not straight lines, so
                           // we need some classes to predefine some 
                           // boundary description
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
                         
#include <fe/mapping_q1.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <dofs/dof_constraints.h>

#include <numerics/matrices.h>
#include <numerics/vectors.h>

#include <numerics/data_out.h>
                           // These are for c++
#include <fstream>
#include <iostream>
#include <sstream>
                                                                  
#include <base/logstream.h>
                             
#include <base/point.h>

                           // @sect3{"The forward problem" class template}

                           // The main class is similar to the wave equation.
                           // The difference is that we add an absorbing 
                           // boundary condition. Because we are only interested
                           // in values at specific locations, we define some
                           // parameters to obtain the coordinates of those 
                           // locations.
template <int dim>
class TATForwardProblem
{
  public:
    TATForwardProblem ();
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve_p ();
    void solve_v ();
    void output_results (const unsigned int timestep_number) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix_p;
    SparseMatrix<double> system_matrix_v;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> boundary_matrix;
                           // Number of refinement
    unsigned int n_refinements;
                           // The acoustic speed in the medium $c_0$
    double acoustic_speed;
                           // This parameter is needed for discritizing
                           // time-dependent problem
    double theta; 

                           // The total data collection time
    double total_time;
 
                           // The size of the time step
    double time_step;  
                           // The detector circullarly scan the target region.
                           // The step size of the detector is in angles
    double step_angle;
                           // The scanning radius
    double radius;  
   
    
    Vector<double>       solution_p;
    Vector<double>       old_solution_p;
    Vector<double>       system_rhs_p;

    Vector<double>       solution_v;
    Vector<double>       old_solution_v;
    Vector<double>       system_rhs_v;        
                                          
      
};

 
                            // Declare a class template for the right hand side
                            // of the pressure potential 
template <int dim>
class RightHandSide_p : public Function<dim> 
{
  public:
    RightHandSide_p () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};

                            // Declare a class template for the right hand side
                            // of the derivative of the pressure potential                  
template <int dim>
class RightHandSide_v : public Function<dim> 
{
  public:
    RightHandSide_v () : Function<dim>() {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};

                            // Declare a class template for the initial values
                            // of the pressure potential
template <int dim>
class InitialValues_p : public Function<dim> 
{
  public:
    InitialValues_p () : Function<dim>() {};
    
  virtual double value (const Point<dim> &p,
			  const unsigned int  component = 0) const;
};

                            // Declare a class template for the initial values
                            // of the derivative of the pressure potential
template <int dim>
class InitialValues_v : public Function<dim> 
{
  public:
    InitialValues_v () : Function<dim>() {};
    
    virtual double value (const Point<dim> &p, 
			  const unsigned int  component = 0) const;
};

                             // Here is the function to set the right hand side
                             // values to be zero for pressure potential
template <int dim>
double RightHandSide_p<dim>::value (const Point<dim> &/*p*/,
				  const unsigned int /*component*/) const 
{
  return 0;
}
                              // Similarly we set the right-hand size of the 
                              // derivative of the pressure potential to be 
                              // zero
template <int dim>
double RightHandSide_v<dim>::value (const Point<dim> &/*p*/,
				  const unsigned int /*component*/) const 
{
  return 0;
}


                           // The sources of the thermoacoustic waves 
                           // are small absorbers. We will compare the 
                           // simulation results with the experimental
                           // data.

template <int dim>
double InitialValues_p<dim>::value (const Point<dim> &p,
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
double InitialValues_v<dim>::value (const Point<dim> &/*p*/,
				   const unsigned int /*component*/) const 
{
 
 return 0;
}

                            // Evaluate point values at arbitrary locations 
                            // In real situation, we collect data by placing
                            // a detector in the medium. By scanning the detector,
                            // we obtain a series projections of the target
                            // from different viewing angles. By using a 
                            // circular radon transform, we can reconstruct
                            // the energy distribution in the target area from 
                            // the measurements obtained by the detectors.
                                 
template <int dim> 
double point_value (const DoFHandler<dim> &dof,
	     const Vector<double>  &fe_function,
	     const Point<dim>      &point)
{                              
                            // Define a map that maps the unit cell to a 
                            // a general grid cell with straight lines in 
                            // dim dimensions
  static const MappingQ1<dim> mapping;
  const FiniteElement<dim>& fe = dof.get_fe();

  Assert(fe.n_components() == 1,
	 ExcMessage ("Finite element is not scalar as is necessary for this function"));

                             // First find the cell in which this point
                             // is, initialize a quadrature rule with
                             // it, and then a FEValues object
                             // The algorithm will first look for the 
                             // surrounding cell on a coarse grid, and
                             // then recersively checking its sibling
                             // cells.
  const typename DoFHandler<dim>::active_cell_iterator cell = GridTools::find_active_cell_around_point (dof, point);
              
  const Point<dim> unit_point = mapping.transform_real_to_unit_cell(cell, point);
  Assert (GeometryInfo<dim>::is_inside_unit_cell (unit_point),
          ExcInternalError());

  const Quadrature<dim> quadrature (unit_point);
  FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
  fe_values.reinit(cell);

                             // Then use this to get at the values of
                             // the given fe_function at this point
  std::vector<double> u_value(1);
  fe_values.get_function_values(fe_function, u_value);

  return u_value[0];
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
	        total_time (0.7),
 		time_step (0.5/std::pow(2.,1.0*n_refinements)/acoustic_speed),  
                step_angle (2.25),  
                radius (0.5)

{}   

                              // This is similar to step-6 except that
                              // the mesh generated is a hyper_ball. We select
                              // hyper_ball instead of hyper_cube because of
                              // our data collection geometry is on a circular in 
                              // 2-D, and on a sphere in 3-D. 
template <int dim>
void TATForwardProblem<dim>::make_grid_and_dofs ()
{
                              // In two dimensional domain. The center of the
                              // circle shall be the point (0,0) and the radius
                              // is 1
  const Point<2> center (0,0);
  GridGenerator::hyper_ball (triangulation, center, 1);
                              // Accordingly, we use hyper ball boundary
                              // instead of hyper cube.
  static const HyperBallBoundary<dim> boundary_description(center);
  triangulation.set_boundary (0,boundary_description);
                              //  The mesh is refined n_refinements times
  triangulation.refine_global (n_refinements);

  std::cout << "   Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

  dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;


  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
                              // We will do the following for both
                              // the pressure potential and its derivative
  system_matrix_p.reinit (sparsity_pattern);
  system_matrix_v.reinit (sparsity_pattern);

  mass_matrix.reinit (sparsity_pattern);
  laplace_matrix.reinit (sparsity_pattern);
  boundary_matrix.reinit (sparsity_pattern);

  solution_p.reinit (dof_handler.n_dofs());
  old_solution_p.reinit (dof_handler.n_dofs());
  system_rhs_p.reinit (dof_handler.n_dofs());

  solution_v.reinit (dof_handler.n_dofs());
  old_solution_v.reinit (dof_handler.n_dofs());
  system_rhs_v.reinit (dof_handler.n_dofs());
 
}


                               // @sect3{ Assemble system}
                               // Because we used absorbing boundary condition in the
                               // simulation, a new boundary matrix is introduced.
                               // We need to assemble boundary matrix. The detailed
                               // description for assembling matrix is discussed in
                               // step-3. 
template <int dim>
void TATForwardProblem<dim>::assemble_system () 
{  
  MatrixCreator::create_mass_matrix (dof_handler, QGauss<dim>(3),
				     mass_matrix);
  MatrixCreator::create_mass_matrix (dof_handler, QTrapez<dim>(),
				     mass_matrix);
  mass_matrix /= 2;

  MatrixCreator::create_laplace_matrix (dof_handler, QGauss<dim>(3),
					laplace_matrix);
  MatrixCreator::create_laplace_matrix (dof_handler, QTrapez<dim>(),
					laplace_matrix);
  laplace_matrix /= 2;

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

                                 // The system matrix of pressure potential 
                                 // is shown in the introduction 
  system_matrix_p = 0;
  system_matrix_p.copy_from (mass_matrix);
  system_matrix_p.add (time_step*time_step*theta*theta*acoustic_speed*acoustic_speed, laplace_matrix);
  system_matrix_p.add (acoustic_speed*theta*time_step, boundary_matrix);
                                 // The system matrix of the derivative
                                 // of the pressure potential is same as 
                                 // the mass matrix
  system_matrix_v = 0;
  system_matrix_v.copy_from (mass_matrix);

}


                                 // We will solve two equations. 
                                 // We first solve for pressure potential
                                 // at a time step, then solve the derivative
                                 // of the pressure potential.
template <int dim>
void TATForwardProblem<dim>::solve_p () 
{
  SolverControl           solver_control (1000, 1e-10);
  SolverCG<>              cg (solver_control);
  cg.solve (system_matrix_p, solution_p, system_rhs_p,
	    PreconditionIdentity());

  std::cout << "   " << solver_control.last_step()
	    << " CG iterations needed to obtain convergence."
	    << std::endl;
}
                       
                                 // To solve the derivative of the pressure potential
template <int dim>
void TATForwardProblem<dim>::solve_v () 
{
  SolverControl           solver_control (1000, 1e-10);
  SolverCG<>              cg (solver_control);
 
  cg.solve (system_matrix_v, solution_v, system_rhs_v,
	    PreconditionIdentity());

  std::cout << "   " << solver_control.last_step()
	    << " CG iterations needed to obtain convergence."
	    << std::endl;
}

                                 // We output the solution for pressure potential
                                 // at each time step in "vtk" format.
template <int dim>
void TATForwardProblem<dim>::output_results (const unsigned int timestep_number) const
{
  
  DataOut<dim> data_out; 
  
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution_p, "P");
  data_out.add_data_vector (solution_v, "V");

  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-"
	   << timestep_number<<".vtk";
  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);

  
}

                                 // This is the main function
                                 // 
                                      
template <int dim>
void TATForwardProblem<dim>::run () 
{
  std::cout << "Solving problem in " << dim << " space dimensions." << std::endl;
  
  make_grid_and_dofs();
  assemble_system ();

  ConstraintMatrix constraints;
  constraints.close();

  VectorTools::project (dof_handler,constraints, 
                            QGauss<dim>(3), InitialValues_p<dim>(),
                            old_solution_p);
  VectorTools::project (dof_handler,constraints, 
                            QGauss<dim>(3), InitialValues_v<dim>(),
                            old_solution_v);


  unsigned int timestep_number = 1;
  unsigned int n_steps;
  unsigned int n_detectors;
  double scanning_angle;

                                // Number of time steps is defined as the
                                // ratio of the total time to the time step                                 
  n_steps=static_cast<unsigned int>(std::floor(total_time/time_step));     
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
  

  for (double time = time_step; time<=total_time; time+=time_step, ++timestep_number)
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

      RightHandSide_p<dim> rhs_function_p;
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
      
      RightHandSide_v<dim> rhs_function_v;
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

      //  output_results (timestep_number);
      
                               //  Evaluate the value at specific locations. 
                               //  For 2-D, it is on a circle. For 1-D, 
                               //  it is a point detector.

      proj_out << time ;

      for (unsigned i=1 ; i<=n_detectors; i++){
        project_dat((timestep_number-1)*n_detectors+i)=point_value (dof_handler,solution_p, 
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
      proj_out.close();
                             
    
}
                              // @sect3{The "main" function}
                              // The main function calls the above functions
                              // in the order of their appearances.

int main ()  
{
  deallog.depth_console (0);
  {
    TATForwardProblem<2> TAT_forward_2d;
    TAT_forward_2d.run ();
  }
  
  return 0;
}
