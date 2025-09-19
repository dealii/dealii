#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <complex>
#include <iostream>

#include "../tests.h"


template <int dim, int degree_p, typename VectorType>
class MatrixFreeTest
{
public:
  using CellIterator = typename DoFHandler<dim>::active_cell_iterator;
  using Number       = double;

  MatrixFreeTest(const MatrixFree<dim, double> &data_in)
    : data(data_in){};

  void
  local_apply(const MatrixFree<dim, double>               &data,
              VectorType                                  &dst,
              const VectorType                            &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, degree_p + 1, degree_p + 2, dim, double> velocity(data,
                                                                        0);
    FEEvaluation<dim, degree_p, degree_p + 2, 1, double> pressure(data, 1);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        velocity.reinit(cell);
        velocity.read_dof_values(src.block(0));
        pressure.reinit(cell);
        pressure.read_dof_values(src.block(1));
        local_integrate(velocity, pressure);

        velocity.distribute_local_to_global(dst.block(0));
        pressure.distribute_local_to_global(dst.block(1));
      }
  }
  void
  local_integrate(
    FEEvaluation<dim, degree_p + 1, degree_p + 2, dim, double> &velocity,
    FEEvaluation<dim, degree_p, degree_p + 2, 1, double>       &pressure) const
  {
    using vector_t = VectorizedArray<double>;
    pressure.evaluate(EvaluationFlags::values);
    velocity.evaluate(EvaluationFlags::gradients);

    for (unsigned int q = 0; q < velocity.n_q_points; ++q)
      {
        SymmetricTensor<2, dim, vector_t> sym_grad_u =
          velocity.get_symmetric_gradient(q);
        vector_t pres = pressure.get_value(q);
        vector_t div  = -trace(sym_grad_u);
        pressure.submit_value(div, q);

        // subtract p * I
        for (unsigned int d = 0; d < dim; ++d)
          sym_grad_u[d][d] -= pres;

        velocity.submit_symmetric_gradient(sym_grad_u, q);
      }
    velocity.integrate(EvaluationFlags::gradients);
    pressure.integrate(EvaluationFlags::values);
  }


  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0;
    data.cell_loop(&MatrixFreeTest<dim, degree_p, VectorType>::local_apply,
                   this,
                   dst,
                   src);
  };

private:
  const MatrixFree<dim, double> &data;
};



template <int dim, int fe_degree>
void
test(const FESystem<dim> &fe_u,
     const FESystem<dim> &fe_p,
     const MPI_Comm       comm_global)
{
  SphericalManifold<dim>                    manifold;
  parallel::distributed::Triangulation<dim> triangulation(comm_global);
  GridGenerator::hyper_shell(triangulation, Point<dim>(), 0.5, 1., 96, true);
  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, manifold);
  triangulation.begin_active()->set_refine_flag();
  triangulation.last()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  triangulation.refine_global(3 - dim);
  triangulation.last()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  MappingQ<dim>   mapping(3);
  DoFHandler<dim> dof_handler_u(triangulation);
  DoFHandler<dim> dof_handler_p(triangulation);
  dof_handler_u.distribute_dofs(fe_u);
  dof_handler_p.distribute_dofs(fe_p);

  MatrixFree<dim, double> mf_data;

  AffineConstraints<double> constraints_u, constraints_p;
  std::vector<const dealii::AffineConstraints<double> *> constraints{
    &constraints_u, &constraints_p};

  IndexSet locally_relevant_dofs_u;
  IndexSet locally_relevant_dofs_p;
  DoFTools::extract_locally_relevant_dofs(dof_handler_u,
                                          locally_relevant_dofs_u);
  DoFTools::extract_locally_relevant_dofs(dof_handler_p,
                                          locally_relevant_dofs_p);
  constraints_u.reinit(dof_handler_u.locally_owned_dofs(),
                       locally_relevant_dofs_u);
  constraints_p.reinit(dof_handler_p.locally_owned_dofs(),
                       locally_relevant_dofs_p);
  const std::set<types::boundary_id> no_normal_flux_boundaries = {0, 1};
  DoFTools::make_hanging_node_constraints(dof_handler_u, constraints_u);
  DoFTools::make_hanging_node_constraints(dof_handler_p, constraints_p);
  if (fe_u.dofs_per_vertex > 0)
    VectorTools::compute_no_normal_flux_constraints(
      dof_handler_u, 0, no_normal_flux_boundaries, constraints_u, mapping);
  constraints_u.close();
  constraints_p.close();

  // setup matrix-free structure
  {
    QGauss<1>                            quad(fe_degree + 2);
    std::vector<const DoFHandler<dim> *> dof_handlers(
      {&dof_handler_u, &dof_handler_p});

    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim, double>::AdditionalData::TasksParallelScheme::none;
    additional_data.mapping_update_flags = update_values | update_gradients;
    // no parallelism
    mf_data.reinit(mapping, dof_handlers, constraints, quad, additional_data);
  }


  LinearAlgebra::distributed::BlockVector<double> solution(2);
  LinearAlgebra::distributed::BlockVector<double> system_rhs(2);
  LinearAlgebra::distributed::BlockVector<double> mf_solution(2);
  for (unsigned int i = 0; i < 2; ++i)
    {
      mf_data.initialize_dof_vector(solution.block(i), i);
      mf_data.initialize_dof_vector(system_rhs.block(i), i);
      mf_data.initialize_dof_vector(mf_solution.block(i), i);
    }

  dealii::TrilinosWrappers::BlockSparsityPattern sparsity_pattern;
  std::vector<dealii::IndexSet>                  row_partitioning{
    dof_handler_u.locally_owned_dofs(), dof_handler_p.locally_owned_dofs()};
  std::vector<dealii::IndexSet> column_partitioning{
    dof_handler_u.locally_owned_dofs(), dof_handler_p.locally_owned_dofs()};
  std::vector<dealii::IndexSet> writeable_rows{locally_relevant_dofs_u,
                                               locally_relevant_dofs_p};
  const auto                    subdomain =
    dealii::Utilities::MPI::this_mpi_process(dof_handler_u.get_communicator());
  sparsity_pattern.reinit(row_partitioning,
                          column_partitioning,
                          writeable_rows,
                          dof_handler_u.get_communicator());

  DoFTools::make_block_sparsity_pattern_block(dof_handler_u,
                                              dof_handler_u,
                                              sparsity_pattern.block(0, 0),
                                              constraints_u,
                                              constraints_u,
                                              true,
                                              subdomain);
  DoFTools::make_block_sparsity_pattern_block(dof_handler_u,
                                              dof_handler_p,
                                              sparsity_pattern.block(0, 1),
                                              constraints_u,
                                              constraints_p,
                                              true,
                                              subdomain);
  DoFTools::make_block_sparsity_pattern_block(dof_handler_p,
                                              dof_handler_u,
                                              sparsity_pattern.block(1, 0),
                                              constraints_p,
                                              constraints_u,
                                              true,
                                              subdomain);
  sparsity_pattern.compress();
  dealii::TrilinosWrappers::BlockSparseMatrix system_matrix;
  system_matrix.reinit(sparsity_pattern);

  // fill system_rhs with random numbers
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < system_rhs.block(i).size(); ++j)
      if (system_rhs.block(i).in_local_range(j) &&
          constraints[i]->is_constrained(j) == false)
        {
          const double val =
            -1 + 2. * (double)Testing::rand() / double(RAND_MAX);
          system_rhs.block(i)(j) = val;
        }


  MatrixFreeTest<dim,
                 fe_degree,
                 LinearAlgebra::distributed::BlockVector<double>>
    mf(mf_data);
  MatrixFreeTools::compute_matrix(
    mf_data,
    constraints,
    system_matrix,
    &MatrixFreeTest<
      dim,
      fe_degree,
      LinearAlgebra::distributed::BlockVector<double>>::local_integrate,
    &mf,
    {0, 1},
    {0, 0});
  system_matrix.vmult(solution, system_rhs);
  mf.vmult(mf_solution, system_rhs);

  // Verification
  mf_solution -= solution;
  const double error    = mf_solution.linfty_norm();
  const double relative = solution.linfty_norm();
  deallog << "Verification " << fe_u.get_name() << "-" << fe_p.get_name()
          << ": " << error / relative << std::endl
          << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;
  MPI_Comm                         comm = MPI_COMM_WORLD;
  {
    deallog << std::endl << "Test with doubles" << std::endl << std::endl;
    deallog.push("2d");
    test<2, 1>(FESystem<2>(FE_Q<2>(2), 2), FESystem<2>(FE_Q<2>(1), 1), comm);
    test<2, 2>(FESystem<2>(FE_Q<2>(3), 2), FESystem<2>(FE_Q<2>(2), 1), comm);
    test<2, 3>(FESystem<2>(FE_Q<2>(4), 2), FESystem<2>(FE_Q<2>(3), 1), comm);
    test<2, 1>(FESystem<2>(FE_DGQ<2>(2), 2),
               FESystem<2>(FE_DGQ<2>(1), 1),
               comm);
    test<2, 0>(FESystem<2>(FE_DGQ<2>(1), 2),
               FESystem<2>(FE_DGQ<2>(0), 1),
               comm);
    test<2, 3>(FESystem<2>(FE_DGQ<2>(4), 2),
               FESystem<2>(FE_DGQ<2>(3), 1),
               comm);
    deallog.pop();
    deallog.push("3d");
    test<3, 1>(FESystem<3>(FE_Q<3>(2), 3), FESystem<3>(FE_Q<3>(1), 1), comm);
    test<3, 1>(FESystem<3>(FE_DGQ<3>(2), 3),
               FESystem<3>(FE_DGQ<3>(1), 1),
               comm);
    deallog.pop();
  }
}
