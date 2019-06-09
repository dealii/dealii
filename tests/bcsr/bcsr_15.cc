// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// check FEEvaluation typetraits with BCSR

#include <deal.II/matrix_free/type_traits.h>

#include <deal.II/lac/block_csr_matrix.h>

#include <fstream>
#include <iostream>

using namespace dealii;


template <typename Number=double>
void test()
{
  deallog << "has_local_element:" << std::endl
          << "BlockCSRMatrixIterators::RowsAccessor<Number, true> = "
          << dealii::internal::has_local_element<
               BlockCSRMatrixIterators::RowsAccessor<Number, true>>::value
          << std::endl
          << "BlockCSRMatrixIterators::RowsAccessor<Number, false> = "
          << dealii::internal::has_local_element<
               BlockCSRMatrixIterators::RowsAccessor<Number, false>>::value
          << std::endl;

  // now check has_partitioners_are_compatible:
  deallog << "has_partitioners_are_compatible:" << std::endl
          << "BlockCSRMatrixIterators::RowsAccessor<Number, true> = "
          << dealii::internal::has_partitioners_are_compatible<
               BlockCSRMatrixIterators::RowsAccessor<Number, true>>::value
          << std::endl
          << "BlockCSRMatrixIterators::RowsAccessor<Number, false> = "
          << dealii::internal::has_partitioners_are_compatible<
               BlockCSRMatrixIterators::RowsAccessor<Number, false>>::value
          << std::endl;

  // check has_begin:
  deallog
    << "has_begin:" << std::endl
    << "BlockCSRMatrixIterators::RowsAccessor<Number, true> = "
    << dealii::internal::has_begin<
          BlockCSRMatrixIterators::RowsAccessor<Number, true>>::value
    << std::endl
    << "BlockCSRMatrixIterators::RowsAccessor<Number, false> = "
    << dealii::internal::has_begin<
          BlockCSRMatrixIterators::RowsAccessor<Number, false>>::value
    << std::endl;

  // check is_vectorizable:
  deallog
    << "is_vectorizable:" << std::endl
    << "BlockCSRMatrixIterators::RowsAccessor<Number, true> = "
    << dealii::internal::is_vectorizable<
          BlockCSRMatrixIterators::RowsAccessor<Number, true>, Number>::value
    << std::endl
    << "BlockCSRMatrixIterators::RowsAccessor<Number, false> = "
    << dealii::internal::is_vectorizable<
          BlockCSRMatrixIterators::RowsAccessor<Number, false>, Number>::value
    << std::endl;


  deallog
    << "has_add_local_element:" << std::endl
    << "BlockCSRMatrixIterators::RowsAccessor<Number, true> = "
    << dealii::internal::has_add_local_element<
          BlockCSRMatrixIterators::RowsAccessor<Number, true>>::value
    << std::endl
    << "BlockCSRMatrixIterators::RowsAccessor<Number, false> = "
    << dealii::internal::has_add_local_element<
          BlockCSRMatrixIterators::RowsAccessor<Number, false>>::value
    << std::endl;

  deallog
    << "has_set_local_element:" << std::endl
    << "BlockCSRMatrixIterators::RowsAccessor<Number, true> = "
    << dealii::internal::has_set_local_element<
          BlockCSRMatrixIterators::RowsAccessor<Number, true>>::value
    << std::endl
    << "BlockCSRMatrixIterators::RowsAccessor<Number, false> = "
    << dealii::internal::has_set_local_element<
          BlockCSRMatrixIterators::RowsAccessor<Number, false>>::value
    << std::endl;

  deallog << "OK" << std::endl;
}

int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  test();
}
