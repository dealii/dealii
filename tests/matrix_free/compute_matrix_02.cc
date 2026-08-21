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
  using FECellIntegrator =
    FEEvaluation<dim, degree_p, degree_p + 1, dim, double>;
  using FEFaceIntegrator =
    FEFaceEvaluation<dim, degree_p, degree_p + 1, dim, double>;

  MatrixFreeTest(const MatrixFree<dim, double> &data_in)
    : data(data_in){};



  void
  do_cell_integral_range(
    const MatrixFree<dim, Number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const
  {
    FECellIntegrator electric(matrix_free, 0);
    FECellIntegrator magnetic(matrix_free, 1);
    for (unsigned int cell = range.first; cell < range.second; ++cell)
      {
        electric.reinit(cell);
        electric.read_dof_values(src.block(0));
        magnetic.reinit(cell);
        magnetic.read_dof_values(src.block(1));

        do_cell_integral_local(electric, magnetic);

        electric.distribute_local_to_global(dst.block(0));
        magnetic.distribute_local_to_global(dst.block(1));
      }
  }


  void
  do_cell_integral_local(FECellIntegrator &electric,
                         FECellIntegrator &magnetic) const
  {
    electric.evaluate(EvaluationFlags::gradients);
    magnetic.evaluate(EvaluationFlags::gradients);

    for (unsigned int q = 0; q < electric.n_q_points; ++q)
      {
        auto curl_h = magnetic.get_curl(q);
        auto curl_e = electric.get_curl(q);
        electric.submit_value(-curl_h, q);
        magnetic.submit_value(curl_e, q);
      }

    electric.integrate(EvaluationFlags::values);
    magnetic.integrate(EvaluationFlags::values);
  }

  void
  do_face_integral_range(
    const MatrixFree<dim, Number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const
  {
    FEFaceIntegrator electric_in = FEFaceIntegrator(matrix_free, true, 0, 0);
    FEFaceIntegrator electric_ex = FEFaceIntegrator(matrix_free, false, 0, 0);
    FEFaceIntegrator magnetic_in = FEFaceIntegrator(matrix_free, true, 1, 0);
    FEFaceIntegrator magnetic_ex = FEFaceIntegrator(matrix_free, false, 1, 0);
    std::pair<FEFaceIntegrator &, FEFaceIntegrator &> electric(electric_in,
                                                               electric_ex);
    std::pair<FEFaceIntegrator &, FEFaceIntegrator &> magnetic(magnetic_in,
                                                               magnetic_ex);

    for (unsigned int face = range.first; face < range.second; ++face)
      {
        electric_in.reinit(face);
        electric_ex.reinit(face);
        electric_in.read_dof_values(src.block(0));
        electric_ex.read_dof_values(src.block(0));
        magnetic_in.reinit(face);
        magnetic_ex.reinit(face);
        magnetic_in.read_dof_values(src.block(1));
        magnetic_ex.read_dof_values(src.block(1));

        do_face_integral_local(electric, magnetic);

        electric_in.distribute_local_to_global(dst.block(0));
        electric_ex.distribute_local_to_global(dst.block(0));
        magnetic_in.distribute_local_to_global(dst.block(1));
        magnetic_ex.distribute_local_to_global(dst.block(1));
      }
  }
  void
  do_face_integral_local(
    const std::pair<FEFaceIntegrator &, FEFaceIntegrator &> &electric,
    const std::pair<FEFaceIntegrator &, FEFaceIntegrator &> &magnetic) const
  {
    FEFaceIntegrator &electric_in = electric.first;
    FEFaceIntegrator &electric_ex = electric.second;
    FEFaceIntegrator &magnetic_in = magnetic.first;
    FEFaceIntegrator &magnetic_ex = magnetic.second;
    electric_in.evaluate(EvaluationFlags::values);
    magnetic_in.evaluate(EvaluationFlags::values);
    electric_ex.evaluate(EvaluationFlags::values);
    magnetic_ex.evaluate(EvaluationFlags::values);

    for (unsigned int q = 0; q < electric_in.n_q_points; ++q)
      {
        auto       normal = electric_in.normal_vector(q);
        const auto jump_exn =
          cross_product_3d(electric_in.get_value(q) - electric_ex.get_value(q),
                           normal);
        const auto jump_hxn =
          cross_product_3d(magnetic_in.get_value(q) - magnetic_ex.get_value(q),
                           normal);
        ;

        magnetic_in.submit_value(-jump_exn * Number(0.5), q);
        magnetic_ex.submit_value(-jump_exn * Number(0.5), q);

        electric_in.submit_value(jump_hxn * Number(0.5), q);
        electric_ex.submit_value(jump_hxn * Number(0.5), q);
      }

    electric_in.integrate(EvaluationFlags::values);
    magnetic_in.integrate(EvaluationFlags::values);
    electric_ex.integrate(EvaluationFlags::values);
    magnetic_ex.integrate(EvaluationFlags::values);
  }

  void
  do_boundary_integral_range(
    const MatrixFree<dim, Number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const
  {
    FEFaceIntegrator electric(matrix_free, true, 0, 0);
    FEFaceIntegrator magnetic(matrix_free, true, 1, 0);
    for (unsigned int face = range.first; face < range.second; ++face)
      {
        electric.reinit(face);
        electric.read_dof_values(src.block(0));
        magnetic.reinit(face);
        magnetic.read_dof_values(src.block(1));
        do_boundary_integral_local(electric, magnetic);
        electric.distribute_local_to_global(dst.block(0));
      }
  }

  void
  do_boundary_integral_local(FEFaceIntegrator &electric,
                             FEFaceIntegrator &magnetic) const
  {
    electric.evaluate(EvaluationFlags::values);
    magnetic.evaluate(EvaluationFlags::values);

    for (unsigned int q = 0; q < electric.n_q_points; ++q)
      {
        auto       normal = electric.normal_vector(q);
        const auto hxn    = cross_product_3d(magnetic.get_value(q), normal);
        electric.submit_value(hxn, q);
      }

    electric.integrate(EvaluationFlags::values);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    dst = 0;
    data.loop(
      &MatrixFreeTest<dim, degree_p, VectorType>::do_cell_integral_range,
      &MatrixFreeTest<dim, degree_p, VectorType>::do_face_integral_range,
      &MatrixFreeTest<dim, degree_p, VectorType>::do_boundary_integral_range,
      this,
      dst,
      src);
  };

private:
  const MatrixFree<dim, double> &data;
};



template <int dim, int fe_degree>
void
test(const FESystem<dim> &fe, const MPI_Comm comm_global)
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
  DoFHandler<dim> dof_handler_e(triangulation);
  DoFHandler<dim> dof_handler_h(triangulation);
  dof_handler_e.distribute_dofs(fe);
  dof_handler_h.distribute_dofs(fe);

  MatrixFree<dim, double> mf_data;

  AffineConstraints<double> constraints_e, constraints_h;
  std::vector<const dealii::AffineConstraints<double> *> constraints{
    &constraints_e, &constraints_h};

  IndexSet locally_relevant_dofs_e;
  IndexSet locally_relevant_dofs_h;
  DoFTools::extract_locally_relevant_dofs(dof_handler_e,
                                          locally_relevant_dofs_e);
  DoFTools::extract_locally_relevant_dofs(dof_handler_h,
                                          locally_relevant_dofs_h);
  constraints_e.reinit(dof_handler_e.locally_owned_dofs(),
                       locally_relevant_dofs_e);
  constraints_h.reinit(dof_handler_h.locally_owned_dofs(),
                       locally_relevant_dofs_h);
  DoFTools::make_hanging_node_constraints(dof_handler_e, constraints_e);
  DoFTools::make_hanging_node_constraints(dof_handler_h, constraints_h);
  constraints_e.close();
  constraints_h.close();

  // setup matrix-free structure
  QGauss<1>                            quad(fe_degree + 1);
  std::vector<Quadrature<1>>           quads{quad, quad};
  std::vector<const DoFHandler<dim> *> dof_handlers(
    {&dof_handler_e, &dof_handler_h});

  // no parallelism
  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::TasksParallelScheme::none;
  additional_data.mapping_update_flags = update_values | update_gradients;
  additional_data.mapping_update_flags_inner_faces =
    update_values | update_normal_vectors | update_quadrature_points;
  additional_data.mapping_update_flags_boundary_faces =
    update_values | update_normal_vectors | update_quadrature_points;
  mf_data.reinit(mapping, dof_handlers, constraints, quads, additional_data);


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
    dof_handler_e.locally_owned_dofs(), dof_handler_h.locally_owned_dofs()};
  std::vector<dealii::IndexSet> column_partitioning{
    dof_handler_e.locally_owned_dofs(), dof_handler_h.locally_owned_dofs()};
  std::vector<dealii::IndexSet> writeable_rows{locally_relevant_dofs_e,
                                               locally_relevant_dofs_h};
  const auto                    subdomain =
    dealii::Utilities::MPI::this_mpi_process(dof_handler_e.get_communicator());
  sparsity_pattern.reinit(row_partitioning,
                          column_partitioning,
                          writeable_rows,
                          dof_handler_e.get_communicator());
  dealii::DoFTools::make_flux_sparsity_pattern(*dof_handlers[0],
                                               sparsity_pattern.block(0, 1),
                                               constraints_h,
                                               true,
                                               subdomain);
  dealii::DoFTools::make_flux_sparsity_pattern(*dof_handlers[1],
                                               sparsity_pattern.block(1, 0),
                                               constraints_e,
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
      LinearAlgebra::distributed::BlockVector<double>>::do_cell_integral_local,
    &MatrixFreeTest<
      dim,
      fe_degree,
      LinearAlgebra::distributed::BlockVector<double>>::do_face_integral_local,
    &MatrixFreeTest<dim,
                    fe_degree,
                    LinearAlgebra::distributed::BlockVector<double>>::
      do_boundary_integral_local,
    &mf,
    {0, 1},
    {0, 0});
  system_matrix.vmult(solution, system_rhs);
  mf.vmult(mf_solution, system_rhs);

  // Verification
  mf_solution -= solution;
  const double error    = mf_solution.linfty_norm();
  const double relative = solution.linfty_norm();
  deallog << "Verification " << fe.get_name() << "-" << fe.get_name() << ": "
          << error / relative << std::endl
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
    deallog.push("3d");
    test<3, 1>(FESystem<3>(FE_DGQ<3>(1), 3), comm);
    test<3, 2>(FESystem<3>(FE_DGQ<3>(2), 3), comm);
    test<3, 3>(FESystem<3>(FE_DGQ<3>(3), 3), comm);
    deallog.pop();
  }
}
