#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/numerics/vector_tools.h>


#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/lac/parpack_solver.h>

#include <fstream>
#include <iostream>

// test Parpack on Step-36 with PETSc algebra

const unsigned int dim = 2;//run in 2d to save time

const double eps = 1e-10;

template <typename DoFHandlerType>
std::vector<dealii::IndexSet>
locally_owned_dofs_per_subdomain (const DoFHandlerType  &dof_handler)
{
  std::vector< dealii::types::subdomain_id > subdomain_association (dof_handler.n_dofs ());
  dealii::DoFTools::get_subdomain_association (dof_handler, subdomain_association);

  const unsigned int n_subdomains = 1 + (*max_element (subdomain_association.begin (),
                                                       subdomain_association.end ()   ));

  std::vector<dealii::IndexSet> index_sets (n_subdomains,dealii::IndexSet(dof_handler.n_dofs()));

  // loop over subdomain_association and populate IndexSet when a
  // change in subdomain ID is found
  dealii::types::global_dof_index i_min          = 0;
  dealii::types::global_dof_index this_subdomain = subdomain_association[0];

  for (dealii::types::global_dof_index index = 1;
       index < subdomain_association.size (); ++index)
    {
      //found index different from the current one
      if (subdomain_association[index] != this_subdomain)
        {
          index_sets[this_subdomain].add_range (i_min, index);
          i_min = index;
          this_subdomain = subdomain_association[index];
        }
    }

  // the very last element is of different index
  if (i_min == subdomain_association.size () - 1)
    {
      index_sets[this_subdomain].add_index (i_min);
    }

  // otherwise there are at least two different indices
  else
    {
      index_sets[this_subdomain].add_range (
        i_min, subdomain_association.size ());
    }

  for (unsigned int i = 0; i < n_subdomains; i++)
    index_sets[i].compress ();

  return index_sets;
} //locally_owned_dofs_per_subdomain

class PETScInverse
{
public:
  PETScInverse(const dealii::PETScWrappers::MatrixBase &A, dealii::SolverControl &cn,const MPI_Comm &mpi_communicator = PETSC_COMM_SELF):
    solver(cn,mpi_communicator),
    matrix(A),
    preconditioner(matrix)
  {

  }

  void vmult  ( dealii::PETScWrappers::MPI::Vector         &dst,
                const dealii::PETScWrappers::MPI::Vector   &src) const
  {
    ;
    solver.solve(matrix, dst, src,preconditioner);
  }


private:
  mutable dealii::PETScWrappers::SolverCG solver;
  const dealii::PETScWrappers::MatrixBase &matrix;
  PETScWrappers::PreconditionBlockJacobi preconditioner;

};

void test ()
{
  const unsigned int global_mesh_refinement_steps = 5;
  const unsigned int number_of_eigenvalues        = 5;

  MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  const unsigned int n_mpi_processes = dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);
  const unsigned int this_mpi_process = dealii::Utilities::MPI::this_mpi_process(mpi_communicator);


  dealii::Triangulation<dim> triangulation;
  dealii::DoFHandler<dim> dof_handler(triangulation);
  dealii::FE_Q<dim> fe(1);
  dealii::ConstraintMatrix constraints;
  dealii::IndexSet locally_owned_dofs;
  dealii::IndexSet locally_relevant_dofs;

  std::vector<dealii::PETScWrappers::MPI::Vector> eigenfunctions;
  std::vector<PetscScalar>                        eigenvalues;
  dealii::PETScWrappers::MPI::SparseMatrix        stiffness_matrix, mass_matrix;

  dealii::GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (global_mesh_refinement_steps);

  // we do not use metis but rather partition by hand below.
  //dealii::GridTools::partition_triangulation (n_mpi_processes, triangulation);
  {
    const double x0 = -1.0;
    const double x1 =  1.0;
    const double dL = (x1-x0) / n_mpi_processes;

    dealii::Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
    for (; cell!=endc; ++cell)
      {
        const dealii::Point<dim> &center = cell->center();
        const double x = center[0];

        const unsigned int id = std::floor ( (x-x0)/dL);
        cell->set_subdomain_id (id);
      }
  }

  dof_handler.distribute_dofs (fe);
  dealii::DoFRenumbering::subdomain_wise (dof_handler);
  std::vector<dealii::IndexSet> locally_owned_dofs_per_processor
    = locally_owned_dofs_per_subdomain (dof_handler);
  locally_owned_dofs = locally_owned_dofs_per_processor[this_mpi_process];
  locally_relevant_dofs.clear();
  dealii::DoFTools::extract_locally_relevant_dofs (dof_handler,
                                                   locally_relevant_dofs);

  constraints.clear();
  constraints.reinit (locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  dealii::VectorTools::interpolate_boundary_values (dof_handler,
                                                    0,
                                                    dealii::ZeroFunction<dim> (),
                                                    constraints);
  constraints.close ();

  dealii::CompressedSimpleSparsityPattern csp (locally_relevant_dofs);
  // Fill in ignoring all cells that are not locally owned
  dealii::DoFTools::make_sparsity_pattern (dof_handler, csp,
                                           constraints,
                                           /* keep constrained dofs */ true);
  std::vector<dealii::types::global_dof_index> n_locally_owned_dofs(n_mpi_processes);
  for (unsigned int i = 0; i < n_mpi_processes; i++)
    n_locally_owned_dofs[i] = locally_owned_dofs_per_processor[i].n_elements();

  dealii::SparsityTools::distribute_sparsity_pattern
  (csp,
   n_locally_owned_dofs,
   mpi_communicator,
   locally_relevant_dofs);

  // Initialise the stiffness and mass matrices
  stiffness_matrix.reinit (locally_owned_dofs,
                           locally_owned_dofs,
                           csp,
                           mpi_communicator);

  mass_matrix.reinit (locally_owned_dofs,
                      locally_owned_dofs,
                      csp,
                      mpi_communicator);

  eigenfunctions.resize (5);
  for (unsigned int i=0; i<eigenfunctions.size (); ++i)
    eigenfunctions[i].reinit (locally_owned_dofs, mpi_communicator);//without ghost dofs

  eigenvalues.resize (eigenfunctions.size ());


  // ready for assembly
  stiffness_matrix = 0;
  mass_matrix = 0;

  dealii::QGauss<dim>   quadrature_formula(2);
  dealii::FEValues<dim> fe_values (fe, quadrature_formula,
                                   dealii::update_values |
                                   dealii::update_gradients |
                                   dealii::update_quadrature_points |
                                   dealii::update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  dealii::FullMatrix<double> cell_stiffness_matrix (dofs_per_cell, dofs_per_cell);
  dealii::FullMatrix<double> cell_mass_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);

  typename dealii::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active (),
  endc = dof_handler.end ();
  for (; cell!=endc; ++cell)
    if (cell->subdomain_id() == this_mpi_process)
      {
        fe_values.reinit (cell);
        cell_stiffness_matrix = 0;
        cell_mass_matrix      = 0;

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                cell_stiffness_matrix (i, j)
                += (fe_values.shape_grad (i, q_point) *
                    fe_values.shape_grad (j, q_point)
                   ) * fe_values.JxW (q_point);

                cell_mass_matrix (i, j)
                += (fe_values.shape_value (i, q_point) *
                    fe_values.shape_value (j, q_point)
                   ) * fe_values.JxW (q_point);
              }

        cell->get_dof_indices (local_dof_indices);

        constraints
        .distribute_local_to_global (cell_stiffness_matrix,
                                     local_dof_indices,
                                     stiffness_matrix);
        constraints
        .distribute_local_to_global (cell_mass_matrix,
                                     local_dof_indices,
                                     mass_matrix);
      }

  stiffness_matrix.compress (dealii::VectorOperation::add);
  mass_matrix.compress (dealii::VectorOperation::add);

  // test Arpack
  {
    std::vector<std::complex<double> > lambda(eigenfunctions.size());

    for (unsigned int i=0; i < eigenvalues.size(); i++)
      eigenfunctions[i] = PetscScalar();

    dealii::SolverControl solver_control (dof_handler.n_dofs(), 1e-9,/*log_history*/false,/*log_results*/false);
    PETScInverse inverse(stiffness_matrix,solver_control,mpi_communicator);
    const unsigned int num_arnoldi_vectors = 2*eigenvalues.size() + 2;

    dealii::PArpackSolver<dealii::PETScWrappers::MPI::Vector>::AdditionalData
    additional_data(num_arnoldi_vectors,
                    dealii::PArpackSolver<dealii::PETScWrappers::MPI::Vector>::largest_magnitude,
                    true);

    dealii::PArpackSolver<dealii::PETScWrappers::MPI::Vector> eigensolver (solver_control,
        mpi_communicator,
        additional_data);
    eigensolver.reinit(locally_owned_dofs);
    eigenfunctions[0] = 1.;
    eigensolver.set_initial_vector(eigenfunctions[0]);
    eigensolver.solve (stiffness_matrix,
                       mass_matrix,
                       inverse,
                       lambda,
                       eigenfunctions,
                       eigenvalues.size());

    for (unsigned int i = 0; i < lambda.size(); i++)
      eigenvalues[i] = lambda[i].real();

    for (unsigned int i=0; i < eigenvalues.size(); i++)
      dealii::deallog << eigenvalues[i] << std::endl;

    // make sure that we have eigenvectors and they are mass-orthonormal:
    // a) (A*x_i-\lambda*B*x_i).L2() == 0
    // b) x_j*B*x_i=\delta_{ij}
    {
      const double precision = 1e-7;
      PETScWrappers::MPI::Vector Ax(eigenfunctions[0]), Bx(eigenfunctions[0]);
      for (unsigned int i=0; i < eigenfunctions.size(); ++i)
        {
          mass_matrix.vmult(Bx,eigenfunctions[i]);

          for (unsigned int j=0; j < eigenfunctions.size(); j++)
            Assert( std::abs( eigenfunctions[j] * Bx - (i==j))< precision,
                    ExcMessage("Eigenvectors " +
                               Utilities::int_to_string(i) +
                               " and " +
                               Utilities::int_to_string(j) +
                               " are not orthonormal!"));

          stiffness_matrix.vmult(Ax,eigenfunctions[i]);
          Ax.add(-1.0*std::real(lambda[i]),Bx);
          Assert (Ax.l2_norm() < precision,
                  ExcMessage(Utilities::to_string(Ax.l2_norm())));
        }
    }
  }


  dof_handler.clear ();
  dealii::deallog << "Ok"<<std::endl;
}


int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile,/*do not print job id*/false);
  dealii::deallog.threshold_double(eps);

  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        test ();
      }

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
