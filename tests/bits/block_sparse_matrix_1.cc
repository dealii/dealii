//----------------------------  block_sparse_matrix_1.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2002, 2003, 2004 by the deal.II authors and Brian Carnes
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_sparse_matrix_1.cc  ---------------------------

// test by Brian: check some of the scaling operations on matrices

#include "../tests.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <base/quadrature_lib.h>
#include <base/function.h>

#include <fe/fe_q_hierarchical.h>
#include <fe/fe_values.h>

#include <dofs/dof_tools.h>

#include <lac/full_matrix.h>
#include <lac/pointer_matrix.h>
#include <lac/block_matrix_array.h>
#include <lac/block_sparse_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/sparsity_pattern.h>
#include <lac/block_vector.h>
#include <lac/vector.h>
#include <lac/solver_cg.h>
#include <lac/sparse_ilu.h>
#include <lac/precondition.h>
#include <lac/solver_control.h>
#include <lac/vector_memory.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>



int main()
{
  std::ofstream logfile("block_sparse_matrix_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  BlockSparseMatrix<double> B;

  Triangulation<2> tria;  
  GridGenerator::hyper_cube (tria,0,1);
  tria.refine_global (1);

  FE_Q_Hierarchical<2> fe (1);
  DoFHandler<2> dof_handler (tria);
  dof_handler.distribute_dofs (fe);
    
  BlockSparsityPattern sparsity_pattern;
  sparsity_pattern.reinit (2,2);
  sparsity_pattern.collect_sizes ();
  
  sparsity_pattern.block(0,0).reinit (dof_handler.n_dofs(),
				      dof_handler.n_dofs(),
				      dof_handler.max_couplings_between_dofs());
  sparsity_pattern.block(0,1).reinit (dof_handler.n_dofs(),
				      1,
				      1);
  sparsity_pattern.block(1,0).reinit (1,
				      dof_handler.n_dofs(),
				      dof_handler.n_dofs());
  sparsity_pattern.block(1,1).reinit (1,
				      1,
				      1);
  sparsity_pattern.collect_sizes ();
  
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern.block(0,0));
  
  for (unsigned int j=0; j<dof_handler.n_dofs (); ++j)
  {
    sparsity_pattern.block(0,1).add (j,0);
    sparsity_pattern.block(1,0).add (0,j);
  }
  sparsity_pattern.block(1,1).add (0,0);
  sparsity_pattern.compress ();
  
  B.reinit (sparsity_pattern);  
  
  QGauss<2> qr (2); 
  FEValues<2> fe_values (dof_handler.get_fe(), 
			 qr, UpdateFlags(update_values    |
					 update_gradients |
					 update_q_points  |
					 update_JxW_values));

  MatrixTools::create_laplace_matrix (dof_handler, qr, B.block(0,0));
  B.block(1,1).add (0,0,1.);

  B.print_formatted (deallog.get_file_stream (),3,false);

  B.block(0,0) *= 10.;
  B.print_formatted (deallog.get_file_stream (),3,false);
  
  B *= 10.;
  B.print_formatted (deallog.get_file_stream (),3,false);
  
  B /= 10.;
  B.print_formatted (deallog.get_file_stream (),3,false);
  
  B.block(0,0) /= 10.;
  B.print_formatted (deallog.get_file_stream (),3,false);

  B = 0;
  
  return 0;
}
