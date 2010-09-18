//----------------------------  create_mass_matrix_03.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  create_mass_matrix_03.cc  ---------------------------


// The MatrixCreator::create_mass_matrix function overload that also assembles
// a right hand side vector had a bug in that the statement that assembled the
// rhs vector was nested in the wrong loop. this was fixed by Moritz' commit
// 14428
//
// the function internally has four branches, with different code used
// for the cases with/without coefficient and scalar/vector-valued
// finite element. we test these four cases through the _01, _02, _03,
// and _04 tests. the version without creating a right hand side vector is tested in the
// _0[1234]a tests, and versions without computing a right hand side
// vectors with and without coefficient in the _0[1234][bc] tests



#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function_lib.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/mapping_q.h>
#include <numerics/matrices.h>

#include <fstream>



template <int dim>
void
check ()
{
  Triangulation<dim> tr;  
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

				   // create a system element composed
				   // of one Q1 and one Q2 element
  FE_Q<dim> element(2);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

				   // use a more complicated mapping
				   // of the domain and a quadrature
				   // formula suited to the elements
				   // we have here
  MappingQ<dim> mapping (3);
  QGauss<dim> quadrature(6);

				   // create sparsity pattern. note
				   // that different components should
				   // not couple, so use pattern
  SparsityPattern sparsity (dof.n_dofs(), dof.n_dofs());
  DoFTools::make_sparsity_pattern (dof, sparsity);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close ();
  constraints.condense (sparsity);
  sparsity.compress ();
  
  SparseMatrix<double> matrix;
  matrix.reinit (sparsity);

  Functions::ExpFunction<dim> coefficient;
  
  MatrixTools::
    create_mass_matrix (mapping, dof,
			quadrature, matrix,
			&coefficient);

				   // since we only generate
				   // output with two digits after
				   // the dot, and since matrix
				   // entries are usually in the
				   // range of 1 or below,
				   // multiply matrix by 100 to
				   // make test more sensitive
  deallog << "Matrix: " << std::endl;
  for (unsigned int i=0; i<matrix.n_nonzero_elements(); ++i)
    deallog << matrix.global_entry(i) * 100
	    << std::endl;
}



int main ()
{
  std::ofstream logfile ("create_mass_matrix_03b/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
