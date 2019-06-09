// check MatrixFree typetraits with BCSR

#include <deal.II/lac/block_vector_base.h>
#include <deal.II/matrix_free/type_traits.h>

#include <RFAStDFT/block_csr_matrix.h>

#include <fstream>
#include <iostream>

using namespace dealii;
using namespace RealFAStDFT;

template <typename Number=double>
void test()
{
  deallog << "has_update_ghost_values_start:" << std::endl
          << "BlockCSRMatrix<Number> = "
          << dealii::internal::has_update_ghost_values_start<
               BlockCSRMatrix<Number>>::value
          << std::endl
          << "has_compress_start:" << std::endl
          << "BlockCSRMatrix<Number> = "
          << dealii::internal::has_compress_start<
               BlockCSRMatrix<Number>>::value
          << std::endl
          << "has_exchange_on_subset:" << std::endl
          << "BlockCSRMatrix<Number> = "
          << dealii::internal::has_exchange_on_subset<
               BlockCSRMatrix<Number>>::value
          << std::endl
          << "has_communication_block_size:" << std::endl
          << "BlockCSRMatrix<Number> = "
          << dealii::internal::has_communication_block_size<
               BlockCSRMatrix<Number>>::value
          << std::endl
          << "is_serial_or_dummy:" << std::endl
          << "BlockCSRMatrix<Number> = "
          << dealii::internal::is_serial_or_dummy<
               BlockCSRMatrix<Number>>::value
          << std::endl
          << "IsBlockVector:" << std::endl
          << "BlockCSRMatrix<Number> = "
          << dealii::IsBlockVector<
               BlockCSRMatrix<Number>>::value
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
