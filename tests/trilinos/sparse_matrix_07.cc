//-----------------  trilinos_sparse_matrix_01.cc  -------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------  trilinos_sparse_matrix_01.cc  -------------------------


// test TrilinosWrappers::SparseMatrix::reinit with a dealii::SparseMatrix
// with a separate sparsity pattern that is partly subset and partly superset

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <fstream>
#include <iostream>


int main (int argc,char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  std::ofstream logfile("sparse_matrix_07/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  SparsityPattern sparsity (5,5,5);
  sparsity.add (1,2);
  sparsity.add (2,3);
  sparsity.add (3,4);
  sparsity.add (4,3);
  sparsity.compress();
  SparseMatrix<double> matrix(sparsity);
  {
    double value = 1;
    for (SparseMatrix<double>::iterator p=matrix.begin();
	 p != matrix.end(); ++p, ++value)
      p->value() = value;
  }
  deallog << "Original:" << std::endl;
  matrix.print_formatted (deallog.get_file_stream());

  // create a separate sparsity pattern to use. note that this sparsity
  // pattern stores the elements explicitly added here but also the diagonal
  // entries (that's what SparsityPattern does for square matrices)
  SparsityPattern xsparsity (5,5,5);
  xsparsity.add (1,2);
  xsparsity.add (2,3);
  xsparsity.add (2,4);
  xsparsity.add (2,1);
  xsparsity.compress();


  // now copy everything into a Trilinos matrix
  Epetra_Map map(5,5,0,Utilities::Trilinos::comm_world());
  TrilinosWrappers::SparseMatrix tmatrix;
  tmatrix.reinit (map, map, matrix, 0, true, &xsparsity);

  deallog << "Copy structure only:" << std::endl;
  tmatrix.print (deallog.get_file_stream());
}
