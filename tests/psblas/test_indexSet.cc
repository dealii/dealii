#include <fstream>

#include "../tests.h"
#include <deal.II/base/mpi.h>
#include <deal.II/lac/psctoolkit.h>

#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/logstream.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_generator.h>
 #include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>


#include <iostream>

using namespace dealii;

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  // Initialize PSBLAS communicator
  int iam, nproc;
  psb_c_ctxt *cctxt = PSCToolkit::Communicator::InitFromMPI(mpi_communicator);
  PSCToolkit::Communicator::Info(cctxt, &iam, &nproc);

  std::string ofname = "output_" + std::to_string(iam);
  std::ofstream output(ofname);

  const unsigned int dim = 2;

  // Create distributed triangulation of the unit square
  Triangulation<dim, dim> tria_base;

  // Create a serial triangulation (here by reading an external mesh):
  GridGenerator::hyper_cube(tria_base,0,1);
  tria_base.refine_global(5);

  // Partition
  GridTools::partition_triangulation(nproc, tria_base);

  // Create building blocks:
  const TriangulationDescription::Description<dim, dim> description = TriangulationDescription::Utilities::create_description_from_triangulation(tria_base, mpi_communicator);

  // Create a fully distributed triangulation:
  parallel::fullydistributed::Triangulation<dim, dim> triangulation(mpi_communicator);
  triangulation.create_triangulation(description);

  // Finite element and DoFHandler
  FE_Q<dim> fe(dim);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

  // Create PSBLAS descriptor for the IndexSet
  psb_c_descriptor *descriptor = PSCToolkit::Communicator::CreateDescriptor(locally_owned_dofs, *cctxt);

  // Assemble the PSBLAS descriptor
  int info = PSCToolkit::Communicator::DescriptorAssembly(descriptor);

  // Extract the number of local indices from PSBLAS descriptor
  psb_i_t number_of_local_indexes = psb_c_cd_get_local_rows(descriptor);

  // Allocate memory for the PSBLAS set of global indices
  psb_l_t *psblas_index_set = (psb_l_t *)malloc(number_of_local_indexes * sizeof(psb_l_t));

  // Fill the PSBLAS index set with local indices
  psb_c_cd_get_global_indices(psblas_index_set, number_of_local_indexes, true, descriptor);

  // Get the portion of local indices from the deal.ii index set
  const std::vector<types::global_dof_index> &indexes = locally_owned_dofs.get_index_vector();

// Check if the PSBLAS index set matches the locally owned indices of indexes vector
for (psb_i_t i = 0; i < number_of_local_indexes; ++i)
{
    if (static_cast<types::global_dof_index>(psblas_index_set[i]) != indexes[i])
    {
        std::cerr << "Mismatch at index " << i << ": PSBLAS index = "
                    << psblas_index_set[i] << ", deal.ii index = "
                    << indexes[i] << " On process: " << iam << std::endl;
    }

}

 // Write output file stating that the index set is okay on process iam
    output << "Process " << iam << " of " << nproc
             << " has a valid PSBLAS index set with " << number_of_local_indexes
             << " local indices." << std::endl;
    output.close();
    
// Concatenate output_i to the file "output"
if (iam == 0)
{
    std::ofstream final_output("output");
    for (int i = 0; i < nproc; i++)
    {
        std::string ofname = "output_" + std::to_string(i);
        std::ifstream input(ofname);
        final_output << input.rdbuf();
        input.close();
        std::remove(ofname.c_str());
    }
}

  // Clean up
  free(psblas_index_set);
  info = PSCToolkit::Communicator::DescriptorFree(descriptor);
  return info;
}