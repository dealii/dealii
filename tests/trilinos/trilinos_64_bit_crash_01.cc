//----------------------------  trilinos_trilinos_64_bit_crash_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005, 2008, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_trilinos_64_bit_crash_01.cc  ---------------------------

// test the BiCGStab solver using the Trilinos matrix and vector classes



#include "../tests.h"
#include <fstream>
#include <iomanip>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>


typedef dealii::types::global_dof_index size_type;

#ifndef DEAL_II_USE_LARGE_INDEX_TYPE
// define a helper function that queries the size of an Epetra_Map object
// by calling either the 32- or 64-bit function necessary, and returns the
// result in the correct data type so that we can use it in calling other
// Epetra member functions that are overloaded by index type
int n_global_elements (const Epetra_BlockMap &map)
{
  return map.NumGlobalElements();
}

int min_my_gid(const Epetra_BlockMap &map)
{
  return map.MinMyGID();
}

int max_my_gid(const Epetra_BlockMap &map)
{
  return map.MaxMyGID();
}

int n_global_cols(const Epetra_CrsGraph &graph)
{
  return graph.NumGlobalCols();
}

int global_column_index(const Epetra_CrsMatrix &matrix, int i)
{
  return matrix.GCID(i);
}

int global_row_index(const Epetra_CrsMatrix &matrix, int i)
{
  return matrix.GRID(i);
}
#else
// define a helper function that queries the size of an Epetra_Map object
// by calling either the 32- or 64-bit function necessary, and returns the
// result in the correct data type so that we can use it in calling other
// Epetra member functions that are overloaded by index type
long long int n_global_elements (const Epetra_BlockMap &map)
{
  return map.NumGlobalElements64();
}

long long int min_my_gid(const Epetra_BlockMap &map)
{
  return map.MinMyGID64();
}

long long int max_my_gid(const Epetra_BlockMap &map)
{
  return map.MaxMyGID64();
}

long long int n_global_cols(const Epetra_CrsGraph &graph)
{
  return graph.NumGlobalCols64();
}

long long int global_column_index(const Epetra_CrsMatrix &matrix, int i)
{
  return matrix.GCID64(i);
}

long long int global_row_index(const Epetra_CrsMatrix &matrix, int i)
{
  return matrix.GRID64(i);
}
#endif


template <typename Sparsity>
void copy_row (const Sparsity        &csp,
	       const size_type        row,
	       std::vector<TrilinosWrappers::types::int_type> &row_indices)
{
  typename Sparsity::row_iterator col_num = csp.row_begin (row);
  for (size_type col=0; col_num != csp.row_end (row); ++col_num, ++col)
    row_indices[col] = *col_num;
}


int main(int argc, char **argv)
{
  std::ofstream logfile("trilinos_64_bit_crash_01/output");
  deallog.attach(logfile);
  deallog << std::setprecision(4);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  CompressedSimpleSparsityPattern sparsity_pattern (1,1);
  sparsity_pattern.add(0,0);
  sparsity_pattern.compress();

  const Epetra_Map input_row_map (static_cast<TrilinosWrappers::types::int_type>(sparsity_pattern.n_rows()),
				  0,
				  Utilities::Trilinos::comm_self());
  const Epetra_Map input_col_map (static_cast<TrilinosWrappers::types::int_type>(sparsity_pattern.n_cols()),
				  0,
				  Utilities::Trilinos::comm_self());


  std_cxx1x::shared_ptr<Epetra_FECrsMatrix> matrix;

  std_cxx1x::shared_ptr<Epetra_Map> column_space_map(new Epetra_Map (input_col_map));

  const size_type first_row = min_my_gid(input_row_map),
		   last_row = max_my_gid(input_row_map)+1;
  std::vector<int> n_entries_per_row(last_row-first_row);

  for (size_type row=first_row; row<last_row; ++row)
    n_entries_per_row[row-first_row] = sparsity_pattern.row_length(row);

  std_cxx1x::shared_ptr<Epetra_CrsGraph> graph;
  graph.reset (new Epetra_CrsGraph (Copy, input_row_map, input_col_map,
				    &n_entries_per_row[0], true));

  TrilinosWrappers::types::int_type row_indices[1] = {0};

  graph->Epetra_CrsGraph::InsertGlobalIndices (size_type(0), 1, &row_indices[0]);

  graph->FillComplete(input_col_map, input_row_map);
  graph->OptimizeStorage();

  deallog << "OK" << std::endl;
}
