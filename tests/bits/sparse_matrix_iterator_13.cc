//----------------------------  sparse_matrix_iterator_13.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix_iterator_13.cc  ---------------------------

// test SparseMatrix::iterator::operator-

#include "../tests.h"
#include <deal.II/lac/sparse_matrix.h>
#include <fstream>
#include <iomanip>


void test ()
{
  SparsityPattern sp (5,5,3);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if (((i+2*j+1) % 3 == 0)
          ||
          (i==j))
        sp.add (i,j);
  sp.compress ();

  SparseMatrix<double> m(sp);

  for (unsigned int row=0; row<sp.n_rows(); ++row)
    Assert (m.begin(row)-m.begin(row) == 0,
	    ExcInternalError());

  for (unsigned int row=0; row<sp.n_rows(); ++row)
    Assert (m.end(row)-m.begin(row) == (int)sp.row_length(row),
	    ExcInternalError());
  for (unsigned int row=0; row<sp.n_rows(); ++row)
    Assert (m.begin(row)-m.end(row) == -(int)sp.row_length(row),
	    ExcInternalError());

  {
    unsigned int counter = 0;
    for (unsigned int row=0; row<sp.n_rows(); ++row)
      {
	Assert (m.begin(row)-m.begin(0) == (int)counter,
		ExcInternalError());
	Assert (m.begin(0)-m.begin(row) == -(int)counter,
		ExcInternalError());
	counter += sp.row_length(row);
      }
  }

  Assert (m.begin() - m.begin(0) == 0, ExcInternalError());
  Assert (m.begin(0) - m.begin() == 0, ExcInternalError());
  Assert (m.end(sp.n_rows()-1) - m.end() == 0, ExcInternalError());
  Assert (m.end() - m.end(sp.n_rows()-1) == 0, ExcInternalError());
  Assert (m.end() - m.begin() == (int)sp.n_nonzero_elements(),
	  ExcInternalError());
  Assert (m.begin() - m.end() == -(int)sp.n_nonzero_elements(),
	  ExcInternalError());

  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("sparse_matrix_iterator_13/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
