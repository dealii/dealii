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

#include <deal.II/meshworker/mesh_loop.h>


#include <iostream>

using namespace dealii;

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  // Initialize PSBLAS communicator
  int iam, nproc, info;
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
  IndexSet locally_relevant_dofs; 
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // Create PSBLAS descriptor for the IndexSet
  psb_c_descriptor *descriptor = PSCToolkit::Communicator::CreateDescriptor(locally_owned_dofs, *cctxt);

  // Create a PSBLAS Sparse matrix on the descriptor and communicator
  psb_c_dspmat *psblas_sparse_matrix = PSCToolkit::Matrix::CreateSparseMatrix(descriptor);

  // Create a PSBLAS vector for the right-hand side 
  psb_c_dvector *psblas_rhs_vector = PSCToolkit::PSBVector::CreateVector(descriptor);

  // Assemble system
  QGauss<dim> quadrature_formula(fe.degree + 1);
  FEValues<dim> fe_values(fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators()){
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        cell_matrix = 0.;
        cell_rhs    = 0.;

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double rhs_value =
              (fe_values.quadrature_point(q_point)[1] >
                    0.5 +
                      0.25 * std::sin(4.0 * numbers::PI *
                                      fe_values.quadrature_point(q_point)[0]) ?
                  1. :
                  -1.);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
                                        fe_values.shape_grad(j, q_point) *
                                        fe_values.JxW(q_point);

                cell_rhs(i) += rhs_value *                         
                                fe_values.shape_value(i, q_point) * 
                                fe_values.JxW(q_point);
              }
          }

        cell->get_dof_indices(local_dof_indices);
        // Distribute local indices and values to the PSBLAS sparse matrix and vector

        info = PSCToolkit::Matrix::distribute_local_to_global(local_dof_indices,
                                                        cell_matrix,
                                                        cell_rhs,
                                                        psblas_sparse_matrix,
                                                        psblas_rhs_vector,
                                                        descriptor);
        if (info != 0)
          deallog << "Error distributing local to global: " << info << std::endl;
      }
    }


// Assemble the Descriptor, the PSBLAS sparse matrix and the PSBLAS vector
info = PSCToolkit::Communicator::DescriptorAssembly(descriptor);
if (info != 0) {
  deallog << "Error assembling PSBLAS descriptor: " << info << std::endl;
}
info = PSCToolkit::Matrix::AssembleSparseMatrix(psblas_sparse_matrix, descriptor);
if (info != 0) {
  deallog << "Error assembling PSBLAS sparse matrix: " << info << std::endl;
}
info = PSCToolkit::PSBVector::AssembleVector(psblas_rhs_vector, descriptor);
if (info != 0) {
  deallog << "Error assembling PSBLAS vector: " << info << std::endl;
}

output << "Process " << iam << " of  " << nproc
       << " I have assembled a matrix of  " << psb_c_dmat_get_nrows(psblas_sparse_matrix) << " x "
       << psb_c_dmat_get_ncols(psblas_sparse_matrix) << " size with "
       << psb_c_dnnz(psblas_sparse_matrix,descriptor) << " non-zero entries."
       << " The locally owned dofs are: " << locally_owned_dofs.n_elements()
       << " and the locally relevant dofs are: " << locally_relevant_dofs.n_elements() << std::endl;
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
  info = PSCToolkit::Matrix::FreeSparseMatrix(psblas_sparse_matrix, descriptor);
  if (info != 0) {
    deallog << "Error freeing PSBLAS sparse matrix: " << info << std::endl;
  }
  info = PSCToolkit::PSBVector::FreeVector(psblas_rhs_vector, descriptor);
  if (info != 0) {
    deallog << "Error freeing PSBLAS vector: " << info << std::endl;
  }
  // Free the PSBLAS descriptor
  info = PSCToolkit::Communicator::DescriptorFree(descriptor);
  if (info != 0) {
    deallog << "Error freeing PSBLAS descriptor: " << info << std::endl;
  }
  return info;
}