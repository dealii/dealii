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

// check MatrixFree typetraits with BCSR

#include <deal.II/lac/block_csr_matrix.h>
#include <deal.II/lac/block_vector_base.h>

#include <deal.II/matrix_free/type_traits.h>

#include <fstream>
#include <iostream>

using namespace dealii;


template <typename Number = double>
void
test()
{
  deallog
    << "has_update_ghost_values_start:" << std::endl
    << "BlockCSRMatrix<Number> = "
    << dealii::internal::has_update_ghost_values_start<
         BlockCSRMatrix<Number>>::value
    << std::endl
    << "has_compress_start:" << std::endl
    << "BlockCSRMatrix<Number> = "
    << dealii::internal::has_compress_start<BlockCSRMatrix<Number>>::value
    << std::endl
    << "has_exchange_on_subset:" << std::endl
    << "BlockCSRMatrix<Number> = "
    << dealii::internal::has_exchange_on_subset<BlockCSRMatrix<Number>>::value
    << std::endl
    << "has_communication_block_size:" << std::endl
    << "BlockCSRMatrix<Number> = "
    << dealii::internal::has_communication_block_size<
         BlockCSRMatrix<Number>>::value
    << std::endl
    << "is_serial_or_dummy:" << std::endl
    << "BlockCSRMatrix<Number> = "
    << dealii::internal::is_serial_or_dummy<BlockCSRMatrix<Number>>::value
    << std::endl
    << "IsBlockVector:" << std::endl
    << "BlockCSRMatrix<Number> = "
    << dealii::IsBlockVector<BlockCSRMatrix<Number>>::value << std::endl;

  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  std::ofstream                            logfile("output");
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  test();
}
