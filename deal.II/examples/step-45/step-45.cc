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

				 // The include files are already known.
#include <base/function.h>
#include <base/quadrature_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <grid/grid_generator.h>
#include <grid/tria.h>
#include <lac/constraint_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/solver_control.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <numerics/data_out.h>
#include <numerics/vectors.h>

using namespace dealii;
				 // The RightHandSide class is a function
				 // object representing the right-hand side of
				 // the problem.
class RightHandSide: public Function<2> {
  public:
    RightHandSide ();
    virtual double value (const Point<2>& p, const unsigned int component = 0) const;
};

				 // This function returns the value of the
				 // right-hand side at a given point
				 // <code>p</code>.
double RightHandSide::value (const Point<2>&p, const unsigned int) const {
  return numbers::PI * numbers::PI * std::sin (numbers::PI * p (0)) * std::sin (numbers::PI * p (1));
}

				 // Here comes the constructor of the class
				 // <code>RightHandSide</code>.
RightHandSide::RightHandSide (): Function<2> () {
}

				 // The class <code>LaplaceProblem</code> is
				 // the main class, which solves the problem.
class LaplaceProblem {
  private:
    ConstraintMatrix constraints;
    const RightHandSide right_hand_side;
    Triangulation<2> triangulation;
    DoFHandler<2> dof_handler;
    FE_Q<2> fe;
    SparseMatrix<double> A;
    SparsityPattern sparsity_pattern;
    Vector<double> b;
    Vector<double> u;
    void assemble_system ();
    void output_results ();
    void setup_system ();
    void solve ();
   
  public:
    LaplaceProblem ();
    void run ();
};

				 // The constructor of the class, where the
				 // <code>DoFHandler</code> and the finite
				 // element object are initialized.
LaplaceProblem::LaplaceProblem (): dof_handler (triangulation), fe (1) {
}

				 //Assembling the system matrix and the
				 //right-hand side vector is done as in other
				 //tutorials before.
void LaplaceProblem::assemble_system () {
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  QGauss<2> quadrature (2);
  const unsigned int n_quadrature_points = quadrature.size ();
  double JxW;
  FEValues<2> fe_values (fe, quadrature, update_gradients | update_JxW_values | update_quadrature_points | update_values);
  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<Point<2> > quadrature_points;
  std::vector<unsigned int> cell_dof_indices (dofs_per_cell);
  Vector<double> cell_rhs (dofs_per_cell);
   
  for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active (); cell != dof_handler.end (); ++cell) {
    cell_rhs = 0;
    fe_values.reinit (cell);
    quadrature_points = fe_values.get_quadrature_points ();
      
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int q_point = 0; q_point < n_quadrature_points; ++q_point) {
	JxW = fe_values.JxW (q_point);
	cell_rhs (i) += JxW * fe_values.shape_value (i, q_point) * right_hand_side.value (quadrature_points[q_point]);
            
	for (unsigned int j = 0; j < dofs_per_cell; ++j)
	  cell_matrix (i, j) += JxW * fe_values.shape_grad (i, q_point) * fe_values.shape_grad (j, q_point);
      }
      
    cell->get_dof_indices (cell_dof_indices);
    constraints.distribute_local_to_global (cell_matrix, cell_rhs, cell_dof_indices, A, b);
  }
}

void LaplaceProblem::setup_system () {
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
  std::cout << "Number of active cells: " << triangulation.n_active_cells () << std::endl << "Degrees of freedom: " << dof_handler.n_dofs () << std::endl;
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
   
  unsigned int dofs_per_face = fe.dofs_per_face;
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
  A.reinit (sparsity_pattern);
  b.reinit (n_dofs);
  u.reinit (n_dofs);
}
				 // To solve the linear system of equations
				 // $Au=b$ we use the CG solver with an
				 // SSOR-preconditioner.
void LaplaceProblem::solve () {
  SolverControl solver_control (dof_handler.n_dofs (), 1e-15);
  PreconditionSSOR<SparseMatrix<double> > precondition;
   
  precondition.initialize (A);
   
  SolverCG<> cg (solver_control);
   
  cg.solve (A, u, b, precondition);
  constraints.distribute (u);
}

void LaplaceProblem::output_results () {
				   // As graphical output we create vtk-file
				   // of the computed solution.
  DataOut<2> data_out;
   
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (u, "u");
  data_out.build_patches ();
   
  std::ofstream output ("solution.vtk");
   
  data_out.write_vtk (output);
}
				 // This function manages the solving process
				 // of the problem.
void LaplaceProblem::run () {
  setup_system ();
  assemble_system ();
  solve ();
  output_results ();
}
				 // And at the end we have the main function
				 // as usual.
int main () {
  LaplaceProblem laplace_problem;
   
  laplace_problem.run ();
}
