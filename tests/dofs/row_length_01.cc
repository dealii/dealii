//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include "../lib/dof_tools_frame.h"
#include "../lib/test_grids.h"

#include <lac/sparsity_pattern.h>
#include <numeric>

template <int dim>
void
check_this (const DoFHandler<dim> &dof)
{
  const unsigned int n = dof.n_dofs();
  deallog << "dofs: " << n
	  << "\tcell: " << dof.get_fe().dofs_per_cell
	  << std::endl;
  
  std::vector<unsigned int> row_length(n);
  DoFTools::compute_row_length_vector(dof, row_length);
  SparsityPattern sparsity(n, n, row_length);
  DoFTools::make_sparsity_pattern(dof, sparsity);
  sparsity.compress();
  unsigned int sum = std::accumulate(row_length.begin(), row_length.end(), 0U);
  deallog << std::endl << "Entries estimated/actual " << sum
	  <<  '/' << sparsity.n_nonzero_elements() << std::endl;
//  sparsity.clear();
  output_vector(row_length);
  
  DoFTools::compute_row_length_vector(dof, row_length, DoFTools::always);
  output_vector(row_length);
}


template <int dim>
void check()
{
  Triangulation<dim> tr;
  TestGrids::hypercube(tr);
  deallog << "cube" << dim << std::endl;
  check_grid(tr);
  tr.refine_global(1);
  deallog << "refined cube" << dim << std::endl;
  check_grid(tr);  
}


int
main()
{
  try
    {
      std::ofstream logfile("row_length_01/output");
      logfile.precision (2);
      deallog.attach(logfile);
//      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);
      check<2>();
      return 0;
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

