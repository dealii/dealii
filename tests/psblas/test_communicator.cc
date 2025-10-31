#include <fstream>

#include "../tests.h"
#include <deal.II/lac/psctoolkit.h>

int main(int argc, char **argv){



    int iam, nproc;
    psb_c_ctxt *cctxt = PSCToolkit::Communicator::Init();
    
    PSCToolkit::Communicator::Info(cctxt, &iam, &nproc);
    
   
    // We don't handle mpi ourselves, so we output stuff manually
    std::string   ofname = "output_" + std::to_string(iam);
    std::ofstream output(ofname);
    output << "Process " << iam << " of " << nproc << " is running." << std::endl;
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

    PSCToolkit::Communicator::Exit(cctxt);

    return 0;
}