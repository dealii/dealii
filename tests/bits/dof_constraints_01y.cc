//----------------------------  dof_constraints_01y.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_constraints_01y.cc  ---------------------------


// check DoFConstraints::distribute_local_to_global for matrices and
// vectors in the presence of both constraints and boundary
// values. this proved to be a really tricky can of worms. this test
// checks things on a rather small domain, the sister test
// dof_constraints_01y does it on complicated meshes

#include "../tests.h"
#include <base/function_lib.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <fe/fe_q.h>
#include <fstream>
#include <iostream>


template <int dim>
void test ()
{
  deallog << dim << "D" << std::endl;
  
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);

				   // set boundary conditions on part of the boundary
  for (unsigned int i=0; i<dim; ++i)
    triangulation.begin_active()->face(i)->set_boundary_indicator (1);

                                   // refine the mesh in a random way so as to
                                   // generate as many hanging node
                                   // constraints as possible
  triangulation.refine_global (4-dim);
  for (unsigned int i=0; i<11-2*dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active();
      for (unsigned int index=0; cell != triangulation.end(); ++cell, ++index)
        if (index % (3*dim) == 0)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement ();
    }
  
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;
  
                                   // set up a DoFHandler and compute hanging
                                   // node constraints
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  deallog << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  constraints.close ();
  deallog << "Number of constraints: " << constraints.n_constraints() << std::endl;

				   // compute some boundary values
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler, 0,
					    ConstantFunction<dim> (1.),
					    boundary_values);
  
                                   // then set up a sparsity pattern and two
                                   // matrices on top of it
  SparsityPattern sparsity (dof_handler.n_dofs(),
                            dof_handler.n_dofs(),
                            dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity);
  constraints.condense (sparsity);
  SparseMatrix<double> A(sparsity), B(sparsity);
  Vector<double> a(dof_handler.n_dofs()), b(dof_handler.n_dofs());

                                   // then fill the two matrices by setting up
                                   // bogus matrix entries and (1) writing
                                   // them into the matrix and condensing away
                                   // hanging node constraints later on, or
                                   // (2) distributing them right away
  std::vector<unsigned int> local_dofs (fe.dofs_per_cell);
  FullMatrix<double> local_matrix (fe.dofs_per_cell, fe.dofs_per_cell);
  Vector<double> local_vector (fe.dofs_per_cell);
  
  for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      cell->get_dof_indices (local_dofs);
      local_matrix = 0;
      local_vector = 0;
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
          local_matrix(i,j) = 1;

      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	local_vector(i) = i;
      
                                       // copy local to global by ourselves
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
          A.add (local_dofs[i], local_dofs[j], local_matrix(i,j));
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	a(local_dofs[i]) += local_vector(i);
      
                                       // or let other functions do that
      MatrixTools::local_apply_boundary_values (boundary_values,
						local_dofs,
						local_matrix,
						local_vector,
						true);
      constraints.distribute_local_to_global (local_matrix, local_dofs,
					      boundary_values, B);
      constraints.distribute_local_to_global (local_vector, local_dofs,
					      boundary_values, b);
    }

                                   // now condense away constraints from A
  constraints.condense (A);
  Vector<double> x (dof_handler.n_dofs());
  MatrixTools::apply_boundary_values (boundary_values, A, x, a);

                                   // we haven't yet set the diagonal entries
                                   // for constrained nodes. we can do so at
                                   // will, since these values don't matter
                                   // anyway
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    if (constraints.is_constrained(i))
      B.set(i,i,A(i,i));

				   // the diagonal elements of
				   // boundary nodes are computed in
				   // different ways, so make sure
				   // that the boundary values will be
				   // correct
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    if (boundary_values.find(i) != boundary_values.end())
      {
	B.diag_element(i) = A.diag_element(i);
	b(i) = a(i);
      }
	
				   // output constraints and boundary values
  deallog << "CONSTRAINTS:" << std::endl;
  constraints.print (deallog.get_file_stream());

  deallog << "BOUNDARY VALUES:" << std::endl;
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    if (boundary_values.find(i) != boundary_values.end())
      deallog << i << ' ' << boundary_values[i]
	      << std::endl;
  
				   // output both matrices. at the
				   // time of writing this test, the
				   // matrices were different and we
				   // needed to figure out what was
				   // going on. since the matrices
				   // aren't big, we just leave the
				   // output here:
  deallog << "MATRIX A:" << std::endl;
  A.print (deallog.get_file_stream());

  deallog << "MATRIX B:" << std::endl;
  B.print (deallog.get_file_stream());
    
                                   // now comes the check: we subtract B from
                                   // A, and make sure that the result is zero
  A.add_scaled (-1., B);

  deallog << "MATRIX DIFFERENCE:" << std::endl;
  A.print (deallog.get_file_stream());

  deallog << "RHS VECTORS:" << std::endl;
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    deallog << a(i) << ' ' << b(i)
	    << (boundary_values.find(i) != boundary_values.end() ? " b" : "  ")
	    << (constraints.is_constrained(i) ? " c" : "")
	    << std::endl;

  deallog << "|A|=" << A.frobenius_norm() << std::endl;
  deallog << "|B|=" << B.frobenius_norm() << std::endl;
  Assert (A.frobenius_norm() < 1e-12*B.frobenius_norm(),
          ExcInternalError());
}



int main ()
{
  std::ofstream logfile("dof_constraints_01y.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      test<2> ();
      test<3> ();
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
}
