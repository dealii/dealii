
#include "ParallelHelperNoMPIStub.hpp"

namespace MESQUITE_NS {

int get_parallel_rank()
{
  int rank=0;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int get_parallel_size()
{
  int nprocs=0;
  //MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  return nprocs;
}

double reduce_parallel_max(double value) {return value;}

void parallel_barrier() {}


}
