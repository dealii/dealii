// check BlockCSRMatrix::copy_to() as well as write accessors to elements in blocks

#include <deal.II/base/logstream.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <RFAStDFT/block_csr_matrix.h>

#include <fstream>
#include <iostream>
#include <numeric>

using namespace RealFAStDFT;
using namespace dealii;

void test ()
{
  // number of blocks:
  const std::vector<unsigned int> row_blocks = {{3,2,1}};
  const std::vector<unsigned int> col_blocks = {{2,2,1}};
  const unsigned int M = row_blocks.size();
  const unsigned int N = col_blocks.size();

  std::vector<types::global_dof_index> row_offset;
  std::vector<types::global_dof_index> col_offset;

  auto setup_offset = [](const std::vector<unsigned int> &blocks,
                         std::vector<types::global_dof_index> &offset) {
      offset.resize(blocks.size()+1, 0);
      std::partial_sum(blocks.begin(), blocks.end(), ++offset.begin());
    };

  setup_offset(row_blocks, row_offset);
  setup_offset(col_blocks, col_offset);


  DynamicSparsityPattern dsp(M,N);
  dsp.add(0,0);
  dsp.add(0,2);
  dsp.add(1,1);
  dsp.add(2,0);
  dsp.add(2,2);

  std::shared_ptr<BlockIndices> rb =
    std::make_shared<BlockIndices>(row_blocks);
  std::shared_ptr<BlockIndices> cb =
    std::make_shared<BlockIndices>(col_blocks);

  auto bcsr_part =
    std::make_shared<dealii::Utilities::MPI::Partitioner>(rb->size());

  // setup matrices
  BlockCSRMatrix<double> A;
  A.reinit(dsp, rb, cb, bcsr_part);

  // setup
  {
    A(row_offset[0]+0,col_offset[0]+0) = 11;
    A(row_offset[0]+1,col_offset[0]+0) = 21;
    A(row_offset[0]+2,col_offset[0]+0) = 31;

    A(row_offset[0]+0,col_offset[0]+1) = 12;
    A(row_offset[0]+1,col_offset[0]+1) = 22;
    A(row_offset[0]+2,col_offset[0]+1) = 32;

    A(row_offset[0]+0,col_offset[2]+0) = 15;
    A(row_offset[0]+1,col_offset[2]+0) = 25;
    A(row_offset[0]+2,col_offset[2]+0) = 35;

    A(row_offset[1]+0,col_offset[1]+0) = 43;
    A(row_offset[1]+1,col_offset[1]+0) = 54;

    A(row_offset[1]+0,col_offset[1]+1) = 44;
    A(row_offset[1]+1,col_offset[1]+1) = 55;

    A(row_offset[2]+0,col_offset[0]+0) = 61;
    A(row_offset[2]+0,col_offset[0]+1) = 62;

    A(row_offset[2]+0,col_offset[2]+0) = 65;
  }

  LAPACKFullMatrix<double> full;
  A.copy_to(full);

  full.print_formatted(deallog.get_file_stream(),0,false, 3);

  deallog << "Ok" << std::endl;
}


int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile,/*do not print job id*/false);
  dealii::deallog.depth_console(0);

  test ();
}
