#include <fstream>

#include "../tests.h"
#include <deal.II/base/mpi.h>
#include <deal.II/lac/psctoolkit.h>

int main(int argc, char **argv){
    int iam, nproc;
    int iam_mpi, nproc_mpi;
   
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    psb_c_ctxt *cctxt = PSCToolkit::Communicator::InitFromMPI(MPI_COMM_WORLD);
    
    nproc_mpi = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    iam_mpi = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  
    PSCToolkit::Communicator::Info(cctxt, &iam, &nproc);
       
    // We don't handle mpi ourselves, so we output stuff manually
    std::string   ofname = "output_" + std::to_string(iam);
    std::ofstream output(ofname);
    output << "Process PSCToolkit " << iam << " of " << nproc << " is running." << std::endl << "MPI Process " << iam_mpi << " of " << nproc_mpi << " is running." << std::endl;
    output.close();    

    // Add MPI barrier
    PSCToolkit::Communicator::Barrier(cctxt);

    // Concatenate output_i to the file "output"
    if (iam == 0)
    {
        std::ofstream final_output("output");
        for (int i = 0; i < nproc; i++)
        {
            std::string   ofname = "output_" + std::to_string(i);
            std::ifstream input(ofname);
            final_output << input.rdbuf();
            input.close();
            std::remove(ofname.c_str());
        }
    }

    return 0;
}