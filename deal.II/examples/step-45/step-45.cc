/* $Id$ */
/* Author: Markus Buerg, University of Karlsruhe, 2010 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2010 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


                                 // @sect3{Include files}

				 // The include files are already known. The
				 // one critical for the current program is
				 // the one that contains the ConstraintMatrix
				 // in the <code>lac/</code> directory:
#include <base/function.h>
#include <base/quadrature_lib.h>

#include <lac/constraint_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/solver_control.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>

#include <grid/grid_generator.h>
#include <grid/tria.h>

#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <numerics/data_out.h>
#include <numerics/vectors.h>

#include <fstream>


using namespace dealii;

                                 // @sect3{The <code>LaplaceProblem</code> class}

				 // The class <code>LaplaceProblem</code> is
				 // the main class of this problem. As
				 // mentioned in the introduction, it is
				 // fashioned after the corresponding class in
				 // step-3. Correspondingly, the documentation
				 // from that tutorial program applies here as
				 // well. The only new member variable is the
				 // <code>constraints</code> variables that
				 // will hold the constraints from the
				 // periodic boundary condition. We will
				 // initialize it in the
				 // <code>setup_system()</code> function.
class LaplaceProblem
{
  public:
    LaplaceProblem ();
    void run ();

  private:
    Triangulation<2> triangulation;

    FE_Q<2> fe;
    DoFHandler<2> dof_handler;

    ConstraintMatrix constraints;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double> system_rhs;
    Vector<double> solution;

    void assemble_system ();
    void output_results ();
    void setup_system ();
    void solve ();   
};


                                 // @sect3{The <code>RightHandSide</code> class}

				 // The following implements the right hand
				 // side function discussed in the
				 // introduction. Its implementation is
				 // obvious given what has been shown in
				 // step-4 before:
class RightHandSide: public Function<2>
{
  public:
    RightHandSide ();

    virtual double value (const Point<2>& p,
			  const unsigned int component = 0) const;
};


RightHandSide::RightHandSide ()
		:
		Function<2> ()
{}


double
RightHandSide::value (const Point<2>&p,
		      const unsigned int) const
{
  return (numbers::PI * numbers::PI *
	  std::sin (numbers::PI * p (0)) *
	  std::sin (numbers::PI * p (1)));
}

                                 // @sect3{Implementation of the <code>LaplaceProblem</code> class}

				 // The first part of implementing the main
				 // class is the constructor. It is unchanged
				 // from step-3 and step-4:
LaplaceProblem::LaplaceProblem ()
		:
		fe (1),
		dof_handler (triangulation)
{}


void LaplaceProblem::setup_system ()
{
  GridGenerator::hyper_cube (triangulation);
				   // We change the boundary indicator on the
				   // parts of the boundary, where we have
				   // Dirichlet boundary conditions, to one
				   // such that we can distinguish between the
				   // parts of the boundary, where periodic
				   // and where Dirichlet boundary conditions
				   // hold.
  Triangulation<2>::active_cell_iterator cell = triangulation.begin_active ();
   
  cell->face (2)->set_boundary_indicator (1);
  cell->face (3)->set_boundary_indicator (1);   
  triangulation.refine_global (5);
				   // Here the degrees of freedom are
				   // distributed.
  dof_handler.distribute_dofs (fe);
  std::cout << "Number of active cells: "
	    << triangulation.n_active_cells ()
	    << std::endl
	    << "Degrees of freedom: " << dof_handler.n_dofs ()
	    << std::endl;
				   // Now it is the time for the constraint
				   // matrix. The first constraints we put in
				   // are the periodic boundary
				   // conditions. For this let us consider the
				   // constraints we have to take care of in
				   // more detail first: We want to identify
				   // all degrees of freedom located on the
				   // right part of the boundary with the ones
				   // located on the left part. Thus, first we
				   // select a degree of freedom on the right
				   // part of the boundary. Then we look for
				   // the corresponding one the left-hand side
				   // and identify these two. Since we are
				   // using finite elements of order 1 here,
				   // finding the corresponding degree of
				   // freedom on the other side is quite easy:
				   // All degrees of freedom are located on
				   // vertices and thus we simply can take the
				   // degree of freedom on the other side,
				   // which is located on the vertex with the
				   // same y-component.  Here starts the
				   // implementation: First we declare a
				   // vector, which stores the global index of
				   // the boundary degree of freedom together
				   // with the y-component of the vertex on
				   // which it is located.
  constraints.clear ();
   
  std::vector<std::pair<unsigned int, double> > dof_locations;
   
  dof_locations.reserve (dof_handler.n_boundary_dofs ());

				   // Then we loop over all active cells and
				   // check, whether the cell is located at
				   // the boundary. If this is the case, we
				   // check, if it is located on the right
				   // part of the boundary, hence, if face 1
				   // is at the boundary.
  for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active (); cell != dof_handler.end (); ++cell)
    if (cell->at_boundary () && cell->face (1)->at_boundary ()) {
				       // Since each degree of freedom of the
				       // face is located on one vertex, we
				       // can access them directly and store
				       // the global index of the degree of
				       // freedom togehter with the
				       // y-component of the vertex on which
				       // it is located.
      dof_locations.push_back (std::pair<unsigned int, double> (cell->vertex_dof_index (1, 0), cell->vertex (1) (1)));
      dof_locations.push_back (std::pair<unsigned int, double> (cell->vertex_dof_index (3, 0), cell->vertex (3) (1)));
    }
				   // Now we have to find the corresponding
				   // degrees of freedom on the left part of
				   // the boundary. Therefore we loop over all
				   // cells again and choose the ones, where
				   // face 0 is at the boundary.
  for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active (); cell != dof_handler.end (); ++cell)
    if (cell->at_boundary () && cell->face (0)->at_boundary ()) {
				       // For every degree of freedom on this
				       // face we add a new line to the
				       // constraint matrix. Then we identify
				       // it with the corresponding degree of
				       // freedom on the right part of the
				       // boundary.
      constraints.add_line (cell->vertex_dof_index (0, 0));
         
      for (unsigned int i = 0; i < dof_locations.size (); ++i)
	if (dof_locations[i].second == cell->vertex (0) (1)) {
	  constraints.add_entry (cell->vertex_dof_index (0, 0), dof_locations[i].first, 1.0);
	  break;
	}
         
      constraints.add_line (cell->vertex_dof_index (2, 0));
         
      for (unsigned int i = 0; i < dof_locations.size (); ++i)
	if (dof_locations[i].second == cell->vertex (2) (1)) {
	  constraints.add_entry (cell->vertex_dof_index (2, 0), dof_locations[i].first, 1.0);
	  break;
	}
    }
				   // Finally we have to set the homogeneous
				   // Dirichlet boundary conditions on the
				   // upper and lower parts of the boundary
				   // and close the
				   // <code>ConstraintMatrix</code> object.
  VectorTools::interpolate_boundary_values (dof_handler, 1, ZeroFunction<2> (), constraints);
  constraints.close ();
				   // Then we create the sparsity pattern and
				   // the system matrix and initialize the
				   // solution and right-hand side vectors.
  const unsigned int n_dofs = dof_handler.n_dofs ();
   
  sparsity_pattern.reinit (n_dofs, n_dofs, dof_handler.max_couplings_between_dofs ());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern, constraints, false);
  sparsity_pattern.compress ();
  system_matrix.reinit (sparsity_pattern);
  system_rhs.reinit (n_dofs);
  solution.reinit (n_dofs);
}


                                 // @sect4{LaplaceProblem::assemble_system}

				 // Assembling the system matrix and the
				 // right-hand side vector is done as in other
				 // tutorials before.
void LaplaceProblem::assemble_system ()
{
  QGauss<2>  quadrature_formula(2);
  FEValues<2> fe_values (fe, quadrature_formula, 
			   update_values   | update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const RightHandSide right_hand_side;

  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
				      endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
				   fe_values.shape_grad (j, q_point) *
				   fe_values.JxW (q_point));

	    cell_rhs(i) += (fe_values.shape_value (i, q_point) *
			    right_hand_side.value (fe_values.quadrature_point (q_point)) *
			    fe_values.JxW (q_point));
	  }
      
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (cell_matrix, cell_rhs,
					      local_dof_indices,
					      system_matrix, system_rhs);
    }
}


                                 // @sect4{LaplaceProblem::solve}

				 // To solve the linear system of equations
				 // $Au=b$ we use the CG solver with an
				 // SSOR-preconditioner. This is, again,
				 // copied almost verbatim from step-4, with
				 // the exception of the preconditioner. As in
				 // step-6, we need to make sure that
				 // constrained degrees of freedom get their
				 // correct values after solving by calling
				 // the ConstraintMatrix::distribute function:
void LaplaceProblem::solve ()
{
  SolverControl solver_control (dof_handler.n_dofs (), 1e-12);
  PreconditionSSOR<SparseMatrix<double> > precondition;
   
  precondition.initialize (system_matrix);
   
  SolverCG<> cg (solver_control);
   
  cg.solve (system_matrix, solution, system_rhs, precondition);
  constraints.distribute (solution);
}


                                 // @sect4{LaplaceProblem::output_results}

				 // This is another function copied from
				 // previous tutorial programs. It generates
				 // graphical output in VTK format:
void LaplaceProblem::output_results ()
{
  DataOut<2> data_out;
   
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "u");
  data_out.build_patches ();
   
  std::ofstream output ("solution.vtk");
   
  data_out.write_vtk (output);
}



                                 // @sect4{LaplaceProblem::run}

				 // And another function copied from previous
				 // programs:
void LaplaceProblem::run () {
  setup_system ();
  assemble_system ();
  solve ();
  output_results ();
}

                                 // @sect3{The <code>main</code> function}

				 // And at the end we have the main function
				 // as usual, this time copied from step-6:
int main ()
{
  try
    {
      deallog.depth_console (0);

      LaplaceProblem laplace_problem;
      laplace_problem.run ();
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
