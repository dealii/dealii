//----------------------------  block_sparse_matrix_iterator_02.cc  ---------------------------
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
//----------------------------  block_sparse_matrix_iterator_02.cc  ---------------------------


// I believed that this test would trigger a bug. Alas, it doesn't, but it
// doesn't hurt to test some anyway

#include "../tests.h"
#include <lac/block_sparsity_pattern.h>
#include <lac/block_sparse_matrix.h>
#include <lac/block_vector.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_q.h>
#include <fstream>
#include <iostream>


void test ()
{
  const int dim = 2;
  
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);

                                   // refine once, then refine first cell to
                                   // create hanging nodes
  triangulation.refine_global (1);
  triangulation.execute_coarsening_and_refinement ();
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;
  
                                   // set up a DoFHandler and compute hanging
                                   // node constraints for a Q2 element
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  deallog << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

                                   // then set up a sparsity pattern and a
                                   // matrix on top of it
  std::vector<unsigned int> block_sizes(2);
  block_sizes[0] = dof_handler.n_dofs()/3;
  block_sizes[1] = dof_handler.n_dofs() - block_sizes[0];

  BlockSparsityPattern sparsity(2,2);
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      sparsity.block(i,j).reinit (block_sizes[i], block_sizes[j],
                                  dof_handler.max_couplings_between_dofs());
  sparsity.collect_sizes();
  
  DoFTools::make_sparsity_pattern (dof_handler, sparsity);
  sparsity.compress ();
  BlockSparseMatrix<double> A(sparsity);

  const BlockSparseMatrix<double>::const_iterator
    begin = A.begin(),
    end = A.end();

  deallog << begin->row() << ' ' << begin->index() << ' '
          << begin->column() << ' ' << begin->block_row() << ' '
          << begin->block_column()
          << std::endl;
  deallog << end->row() << ' ' << end->index() << ' '
          << end->column() << ' ' << end->block_row() << ' '
          << end->block_column()
          << std::endl;

                                   // this matrix certainly has entries
  Assert (begin != end, ExcInternalError());
  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("block_sparse_matrix_iterator_02.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      test ();
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
