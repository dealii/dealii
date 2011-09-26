//-----------------  trilinos_sparse_matrix_01.cc  -------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------  trilinos_sparse_matrix_01.cc  -------------------------


// TrilinosWrappers::SparseMatrix::el() and operator() had problems
// looking up indices in rectangular matrices because they
// accidentally used the row instead of the column map. verify that
// this is now fixed
//
// this testcase is reduced from one contributed by Habib Talavatifard

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <fstream>
#include <iostream>


int main (int argc,char **argv)
{
  std::ofstream logfile("sparse_matrix_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  IndexSet row_partitioning (3);
  IndexSet col_partitioning (4);

  row_partitioning.add_range(0, 3);
  col_partitioning.add_range(0, 4);

				   // Add element (2,3) to the matrix
  TrilinosWrappers::SparsityPattern sp (row_partitioning, col_partitioning);
  sp.add (2,3);
  sp.compress();

  TrilinosWrappers::SparseMatrix A (sp);
  A.add (2, 3, 2.0);
  A.compress();

				   // verify that entry (2,3) is
				   // indeed what we expect. verify
				   // that both methods of accessing
				   // the entry work
  Assert (A.el(2, 3) == 2, ExcInternalError());
  Assert (A(2, 3) == 2, ExcInternalError());

  deallog << "OK" << std::endl;
}
