//----------------------------  trilinos_sparse_matrix_iterator_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_sparse_matrix_iterator_01.cc  ---------------------------


// test TrilinosWrappers::MatrixBase::const_iterator

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <fstream>
#include <iostream>


void test ()
{
  TrilinosWrappers::SparseMatrix m(5U,5U,5U);
  m.set (0,0,1);
  m.set (1,1,2);
  m.set (1,2,3);
  m.compress ();
  TrilinosWrappers::SparseMatrix::const_iterator i = m.begin();
  deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  ++i;
  deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;
  i++;
  deallog << i->row() << ' ' << i->column() << ' ' << i->value() << std::endl;

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("sparse_matrix_iterator_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
        test ();
      }
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
