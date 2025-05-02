
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/tensor_product_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_patch_evaluation.h>
#include <deal.II/matrix_free/patch_storage.h>
#include <deal.II/matrix_free/patch_distributors.h>


#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>



#include <fstream>
#include <iostream>



namespace Operators
{
  using namespace dealii;

  template <int dim, int fe_degree, typename number>
  class LaplaceOperator
    : public MatrixFreeOperators::
        Base<dim, LinearAlgebra::distributed::Vector<number>>
  {
  public:
    using value_type = number;
    using VectorType = LinearAlgebra::distributed::Vector<number>;


    LaplaceOperator();

    void clear() override;

    virtual void compute_diagonal() override;

  private:
    virtual void apply_add(VectorType       &dst,
                           const VectorType &src) const override;

    void
    local_apply(const MatrixFree<dim, number>               &data,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const;

    void local_compute_diagonal(
      const MatrixFree<dim, number>               &data,
      VectorType                                  &dst,
      const unsigned int                          &dummy,
      const std::pair<unsigned int, unsigned int> &cell_range) const;
  };

  template <int dim, int fe_degree, typename number>
  LaplaceOperator<dim, fe_degree, number>::LaplaceOperator()
    : MatrixFreeOperators::Base<dim, VectorType>()
  {}



  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::clear()
  {
    MatrixFreeOperators::Base<dim, VectorType>::clear();
  }


  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::local_apply(
    const MatrixFree<dim, number>               &data,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);
        phi.evaluate(EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate(EvaluationFlags::gradients);
        phi.distribute_local_to_global(dst);
      }
  }


  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::apply_add(
    VectorType       &dst,
    const VectorType &src) const
  {
    this->data->cell_loop(&LaplaceOperator::local_apply, this, dst, src);
  }
  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::compute_diagonal()
  {
    this->inverse_diagonal_entries.reset(new DiagonalMatrix<VectorType>());
    VectorType &inverse_diagonal = this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);
    unsigned int dummy = 0;
    this->data->cell_loop(&LaplaceOperator::local_compute_diagonal,
                          this,
                          inverse_diagonal,
                          dummy);

    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i) =
          1. / inverse_diagonal.local_element(i);
      }
  }

  template <int dim, int fe_degree, typename number>
  void LaplaceOperator<dim, fe_degree, number>::local_compute_diagonal(
    const MatrixFree<dim, number> &data,
    VectorType                    &dst,
    const unsigned int &,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data);

    AlignedVector<VectorizedArray<number>> diagonal(phi.dofs_per_cell);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
              phi.submit_dof_value(VectorizedArray<number>(), j);
            phi.submit_dof_value(make_vectorized_array<number>(1.), i);

            phi.evaluate(EvaluationFlags::gradients);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate(EvaluationFlags::gradients);
            diagonal[i] = phi.get_dof_value(i);
          }
        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi.submit_dof_value(diagonal[i], i);
        phi.distribute_local_to_global(dst);
      }
  }



  template <int dim, int fe_degree, typename number, typename VectorizedValues>
  class PatchTensorInverse
  {
  public:
    //  using value_type = number;


    const constexpr static unsigned int n_rows_1d = 2 * fe_degree - 1;

    void clear()
    {}

    void initialize(std::shared_ptr<const MatrixFree<dim, number>> mf);

    template <typename VectorType>
    void vmult(VectorType &dst, const VectorType &src) const;

    void Tvmult(LinearAlgebra::distributed::Vector<number>       &dst,
                const LinearAlgebra::distributed::Vector<number> &src,
                AlignedVector<VectorizedValues> &tmp_memory) const
    {
      vmult(dst, src, tmp_memory);
    }

  private:
    std::shared_ptr<const MatrixFree<dim, number>> data;

    // replace number with VectorizedArray<Number> if you want
    // vectorization over patches.
    TensorProductMatrixSymmetricSum<dim, VectorizedValues, n_rows_1d>
      patch_matrix;

    std::array<Table<2, VectorizedValues>, dim> patch_mass;
    std::array<Table<2, VectorizedValues>, dim> patch_laplace;
  };

  template <int dim, int fe_degree, typename number, typename VectorizedValues>
  void PatchTensorInverse<dim, fe_degree, number, VectorizedValues>::initialize(
    std::shared_ptr<const MatrixFree<dim, number>> mf)
  {
    data = mf;

    auto fe_1d = std::make_unique<FE_DGQ<1>>(fe_degree);

    const unsigned int                          N = fe_degree + 1;
    FullMatrix<number>                          laplace_unscaled(N, N);
    std::array<Table<2, VectorizedValues>, dim> mass_matrices;
    std::array<Table<2, VectorizedValues>, dim> laplace_matrices;
    for (unsigned int d = 0; d < dim; ++d)
      {
        mass_matrices[d].reinit(N, N);
        laplace_matrices[d].reinit(N, N);
      }

    QGauss<1> quadrature(N);
    // assemble 1d matrices on reference cell
    for (unsigned int i = 0; i < N; ++i)
      for (unsigned int j = 0; j < N; ++j)
        {
          double sum_mass = 0, sum_laplace = 0;
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              sum_mass += (fe_1d->shape_value(i, quadrature.point(q)) *
                           fe_1d->shape_value(j, quadrature.point(q))) *
                          quadrature.weight(q);
              sum_laplace += (fe_1d->shape_grad(i, quadrature.point(q))[0] *
                              fe_1d->shape_grad(j, quadrature.point(q))[0]) *
                             quadrature.weight(q);
            }
          for (unsigned int d = 0; d < dim; ++d)
            mass_matrices[d](i, j) = sum_mass;

          laplace_unscaled(i, j) = sum_laplace;
        }

    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*data);

    phi.reinit(0);

    auto inverse_jacobian = phi.inverse_jacobian(0);

    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        if (d != e)
          AssertThrow(inverse_jacobian[d][e][0] == 0., ExcNotImplemented());

    // number -> VectorizedArray<Number> for vectorization
    // Fixme ? not sure if number ...
    number jacobian_determinant = inverse_jacobian[0][0][0];
    for (unsigned int e = 1; e < dim; ++e)
      jacobian_determinant *= inverse_jacobian[e][e][0];
    jacobian_determinant = 1. / jacobian_determinant;


    for (unsigned int d = 0; d < dim; ++d)
      {
        // FIXME later:
        const number scaling_factor = inverse_jacobian[d][d][0] *
                                      inverse_jacobian[d][d][0] *
                                      jacobian_determinant;

        for (unsigned int i = 0; i < N; ++i)
          for (unsigned int j = 0; j < N; ++j)
            laplace_matrices[d](i, j) =
              scaling_factor * VectorizedValues(laplace_unscaled(i, j));
      }

    const unsigned int NN = n_rows_1d;

    for (unsigned int d = 0; d < dim; ++d)
      {
        patch_mass[d].reinit(NN, NN);
        patch_laplace[d].reinit(NN, NN);
      }

    // assemble tensor product on one patch
    for (unsigned int p = 0; p < 2; ++p)
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int i = 0; i < fe_degree; ++i)
          for (unsigned int j = 0; j < fe_degree; ++j)
            {
              patch_mass[d](i + p * (fe_degree - 1), j + p * (fe_degree - 1)) +=
                mass_matrices[d](i + 1 - p, j + 1 - p);
              patch_laplace[d](i + p * (fe_degree - 1),
                               j + p * (fe_degree - 1)) +=
                laplace_matrices[d](i + 1 - p, j + 1 - p);
            }

    patch_matrix.reinit(patch_mass, patch_laplace);

#ifdef DEBUG_INVERSE
    for (unsigned int i = 0; i < NN; ++i)
      {
        for (unsigned int j = 0; j < NN; ++j)
          std::cout << patch_mass[0](i, j) << " ";
      }
    std::cout << std::endl;
#endif
  }

  template <int dim, int fe_degree, typename number, typename VectorizedValues>
  template <typename VectorType>
  void PatchTensorInverse<dim, fe_degree, number, VectorizedValues>::vmult(
    VectorType       &dst,
    const VectorType &src) const
  {
    const unsigned int n_dofs_interior =
      Utilities::pow<unsigned int>(n_rows_1d, dim);
    patch_matrix.apply_inverse(
      ArrayView<VectorizedValues>(dst.data(), n_dofs_interior),
      ArrayView<const VectorizedValues>(src.data(), n_dofs_interior));
  }



  template <int dim, int fe_degree, typename number>
  class LaplacePatchSmoother : public ::PatchOperator::SmootherBase<dim, number>
  {
  public:
    using BaseType         = PatchOperator::SmootherBase<dim, number>;
    using VectorType       = LinearAlgebra::distributed::Vector<number>;
    using MatrixFreeType   = MatrixFree<dim, number>;
    using PatchStorageType = PatchOperator::PatchStorage<MatrixFreeType>;
    using InverseType      = PatchTensorInverse<dim, fe_degree, number, number>;

    using AdditionalData = std::shared_ptr<PatchStorageType>;

    void local_apply(const PatchStorageType                      &patch_storage,
                     VectorType                                  &dst,
                     const VectorType                            &src,
                     const typename PatchStorageType::PatchRange &patch_range,
                     const bool &do_forward) const override;

    void initialize(std::shared_ptr<PatchStorageType> &patch_storage) override;

    void clear() override;

  private:
    InverseType inverse;
  };

  template <int dim, int fe_degree, typename number>
  void LaplacePatchSmoother<dim, fe_degree, number>::initialize(
    std::shared_ptr<PatchStorageType> &patch_storage)
  {
    BaseType::initialize(patch_storage);
    if (patch_storage->n_patches() != 0)
      inverse.initialize(patch_storage->get_matrix_free());
  }


  template <int dim, int fe_degree, typename number>
  void LaplacePatchSmoother<dim, fe_degree, number>::clear()
  {
    BaseType::clear();
    inverse.clear();
  }

  template <int dim, int fe_degree, typename number>
  void LaplacePatchSmoother<dim, fe_degree, number>::local_apply(
    const PatchStorageType                      &patch_storage,
    VectorType                                  &dst,
    const VectorType                            &src,
    const typename PatchStorageType::PatchRange &patch_range,
    const bool                                  &do_forward) const
  {
    using FEEval = FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number>;
    using PatchEval =
      PatchOperator::FEPatchEvaluation<FEEval,
                                       DistributorLookup<dim, fe_degree>>;


    PatchEval patch_eval(patch_storage,
                         FEEval(*patch_storage.get_matrix_free()));

    const auto begin = do_forward ? patch_range.first : patch_range.second;
    const auto end   = do_forward ? patch_range.second : patch_range.first;
    const auto step  = do_forward ? 1 : -1;


    for (auto patch_index = begin; patch_index != end; patch_index += step)
      {
        patch_eval.reinit(patch_index);

        AlignedVector<number> local_rhs(patch_eval.n_patch_dofs());
        AlignedVector<number> local_residual(patch_eval.n_patch_dofs());
        AlignedVector<number> local_correction(patch_eval.n_patch_dofs());

        patch_eval.read_dof_values(src);

        patch_eval.gather_local_to_patch(ArrayView<double>(local_rhs), false);

        patch_eval.read_dof_values(dst);
        for (auto &phi : patch_eval.fe_evaluations)
          {
            phi.evaluate(EvaluationFlags::gradients);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_gradient(phi.get_gradient(q), q);
            phi.integrate(EvaluationFlags::gradients);
          }

        patch_eval.gather_local_to_patch(ArrayView<double>(local_residual),
                                         true);

        // residual = A*solution - rhs
        for (unsigned int i = 0; i < local_residual.size(); ++i)
          local_residual[i] -= local_rhs[i];


        inverse.vmult(local_correction, local_residual);

        // multiply local_correction by -1
        for (unsigned int i = 0; i < local_correction.size(); ++i)
          local_correction[i] *= -1;
        patch_eval.distribute_patch_to_local(
          ArrayView<const number>(local_correction), false);


        patch_eval.distribute_local_to_global(dst);
      }
  }


} // namespace Operators


namespace LaplaceSolver
{
  using namespace dealii;


  template <int dim, int fe_degree>
  class LaplaceSolver
  {
  public:
    using Number              = double;
    using LevelNumber         = double;
    using VectorType          = LinearAlgebra::distributed::Vector<Number>;
    using LevelVectorType     = LinearAlgebra::distributed::Vector<LevelNumber>;
    using LevelMatrixFreeType = MatrixFree<dim, LevelNumber>;
    LaplaceSolver();
    void initialize(const unsigned int n_refinements = 4);

    void print_timings() const;

    void solve();

    void output_results(const unsigned int cycle) const;

  private:
    void setup_system();
    void assemble_rhs();

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim> triangulation;
#else
    Triangulation<dim> triangulation;
#endif

    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;

    MappingQ1<dim> mapping;

    AffineConstraints<Number> constraints;
    using SystemMatrixType = Operators::LaplaceOperator<dim, fe_degree, Number>;
    SystemMatrixType system_matrix;

    MGConstrainedDoFs mg_constrained_dofs;
    using LevelMatrixType =
      Operators::LaplaceOperator<dim, fe_degree, LevelNumber>;
    MGLevelObject<LevelMatrixType> mg_matrices;

    VectorType solution;
    VectorType system_rhs;

    double             setup_time;
    ConditionalOStream pcout;
    ConditionalOStream time_details;


    mutable std::vector<std::vector<
      std::pair<double, std::chrono::time_point<std::chrono::system_clock>>>>
      all_mg_timers;
  };



  template <int dim, int fe_degree>
  LaplaceSolver<dim, fe_degree>::LaplaceSolver()
    :
#ifdef DEAL_II_WITH_P4EST
    triangulation(
      MPI_COMM_WORLD,
      Triangulation<dim>::limit_level_difference_at_vertices,
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
    ,
#else
    triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    ,
#endif
    fe(fe_degree)
    , dof_handler(triangulation)
    , setup_time(0.)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    , time_details(std::cout,
                   false &&
                     Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {}


  template <int dim, int fe_degree>
  void LaplaceSolver<dim, fe_degree>::setup_system()
  {
    Timer time;
    setup_time = 0;

    system_matrix.clear();
    mg_matrices.clear_elements();

    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();

    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs()
          << std::endl;

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
    constraints.close();
    setup_time += time.wall_time();
    time_details << "Distribute DoFs & B.C.     (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
    time.restart();

    {
      typename MatrixFree<dim, Number>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, Number>::AdditionalData::none;
      additional_data.mapping_update_flags =
        (update_gradients | update_JxW_values | update_quadrature_points);
      std::shared_ptr<MatrixFree<dim, Number>> system_mf_storage(
        new MatrixFree<dim, Number>());
      system_mf_storage->reinit(mapping,
                                dof_handler,
                                constraints,
                                QGauss<1>(fe.degree + 1),
                                additional_data);


      system_mf_storage->reinit(mapping,
                                dof_handler,
                                constraints,
                                QGauss<1>(fe.degree + 1),
                                additional_data);

      system_matrix.initialize(system_mf_storage);
    }

    system_matrix.initialize_dof_vector(solution);
    system_matrix.initialize_dof_vector(system_rhs);

    setup_time += time.wall_time();
    time_details << "Setup matrix-free system   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
    time.restart();

    const unsigned int nlevels = triangulation.n_global_levels();
    mg_matrices.resize(0, nlevels - 1);

    std::set<types::boundary_id> dirichlet_boundary;
    dirichlet_boundary.insert(0);
    mg_constrained_dofs.initialize(dof_handler);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                       dirichlet_boundary);

    for (unsigned int level = 0; level < nlevels; ++level)
      {
        IndexSet relevant_dofs;
        DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                      level,
                                                      relevant_dofs);
        AffineConstraints<Number> level_constraints;
        level_constraints.reinit(relevant_dofs);
        level_constraints.add_lines(
          mg_constrained_dofs.get_boundary_indices(level));
        level_constraints.close();

        typename MatrixFree<dim, LevelNumber>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim, LevelNumber>::AdditionalData::none;
        additional_data.mapping_update_flags =
          (update_gradients | update_JxW_values | update_quadrature_points);
        additional_data.mg_level          = level;
        additional_data.store_ghost_cells = true;
        std::shared_ptr<MatrixFree<dim, LevelNumber>> mg_mf_storage_level(
          new MatrixFree<dim, LevelNumber>());
        mg_mf_storage_level->reinit(mapping,
                                    dof_handler,
                                    level_constraints,
                                    QGauss<1>(fe.degree + 1),
                                    additional_data);

        mg_matrices[level].initialize(mg_mf_storage_level,
                                      mg_constrained_dofs,
                                      level);

        mg_matrices[level].compute_diagonal();
      }
    setup_time += time.wall_time();
    time_details << "Setup matrix-free levels   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
  }



  template <int dim, int fe_degree>
  void LaplaceSolver<dim, fe_degree>::assemble_rhs()
  {
    Timer time;

    system_rhs = 0;
    FEEvaluation<dim, fe_degree> phi(*system_matrix.get_matrix_free());
    for (unsigned int cell = 0;
         cell < system_matrix.get_matrix_free()->n_cell_batches();
         ++cell)
      {
        phi.reinit(cell);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_value(make_vectorized_array<Number>(1.0), q);
        phi.integrate(EvaluationFlags::values);
        phi.distribute_local_to_global(system_rhs);
      }
    system_rhs.compress(VectorOperation::add);

    setup_time += time.wall_time();
    time_details << "Assemble right hand side   (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s" << std::endl;
  }

  template <int dim, int fe_degree>
  void LaplaceSolver<dim, fe_degree>::solve()
  {
    using StorageType = PatchOperator::PatchStorage<MatrixFree<dim, double>>;



    const unsigned max_lvl = triangulation.n_global_levels() - 1;
    MGTransferMatrixFree<dim, LevelNumber> mg_transfer(mg_constrained_dofs);
    mg_transfer.build(dof_handler);

    using SmootherType =
      Operators::LaplacePatchSmoother<dim, fe_degree, double>;
    mg::SmootherRelaxation<SmootherType, LevelVectorType> mg_smoother;
    mg_smoother.resize(0, max_lvl);

    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      {
        std::shared_ptr<StorageType> storage =
          std::make_shared<StorageType>(mg_matrices[level].get_matrix_free());
        // FIXME: Assert not initialized
        storage->initialize();

        mg_smoother[level].initialize(storage);
      }

    mg_smoother[max_lvl].get_storage()->output_patches("patches");
    mg_smoother[max_lvl].get_storage()->output_centerpoints("centerpoints");

    mg_smoother.set_steps(1);

    MGCoarseGridApplySmoother<LevelVectorType> mg_coarse;
    mg_coarse.initialize(mg_smoother);

    mg::Matrix<LevelVectorType> mg_matrix(mg_matrices);

    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
      mg_interface_matrices;
    mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
    for (unsigned int level = 0; level < triangulation.n_global_levels();
         ++level)
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    mg::Matrix<LevelVectorType> mg_interface(mg_interface_matrices);

    Multigrid<LevelVectorType> mg(
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
    mg.set_edge_matrices(mg_interface, mg_interface);

    PreconditionMG<dim, LevelVectorType, MGTransferMatrixFree<dim, LevelNumber>>
      preconditioner(dof_handler, mg, mg_transfer);

    SolverControl        solver_control(15, 1e-8 * system_rhs.l2_norm());
    SolverCG<VectorType> cg(solver_control);

    Timer time;
    time.reset();
    constraints.set_zero(solution);
    // mg_smoother[max_lvl].step(solution, system_rhs);
    cg.solve(system_matrix, solution, system_rhs, preconditioner);
    constraints.distribute(solution);

    pcout << "Time solve (" << solver_control.last_step() << " iterations)"
          << (solver_control.last_step() < 10 ? "  " : " ") << "(CPU/wall) "
          << time.cpu_time() << "s/" << time.wall_time() << "s\n";
  }



  template <int dim, int fe_degree>
  void
  LaplaceSolver<dim, fe_degree>::output_results(const unsigned int cycle) const
  {
    Timer time;
    if (triangulation.n_global_active_cells() > 1000000)
      return;

    DataOut<dim> data_out;

    solution.update_ghost_values();
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");
    data_out.build_patches(mapping, fe_degree);

    DataOutBase::VtkFlags flags;
    // flags.compression_level = DataOutBase::VtkFlags::best_speed;
    data_out.set_flags(flags);
    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, MPI_COMM_WORLD, 3);

    time_details << "Time write output          (CPU/wall) " << time.cpu_time()
                 << "s/" << time.wall_time() << "s\n";
  }


  template <int dim, int fe_degree>
  void
  LaplaceSolver<dim, fe_degree>::initialize(const unsigned int n_refinements)
  {
    {
      const unsigned int n_vect_doubles = VectorizedArray<double>::size();
      const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

      pcout << "Vectorization over " << n_vect_doubles
            << " doubles = " << n_vect_bits << " bits ("
            << Utilities::System::get_current_vectorization_level() << ")"
            << std::endl;

      // pcout << "Running with  "
      // << n_refinements
      // << "  refinements and task_chunk_size = " << task_chunk_size
      // << std::endl;
    }


    {
      GridGenerator::subdivided_hyper_cube(triangulation, 2);
    }
    triangulation.refine_global(n_refinements);
    setup_system();
    assemble_rhs();
    pcout << std::endl;
  }


} // namespace LaplaceSolver

int main(int argc, char *argv[])
{
  using namespace dealii;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  const unsigned int dim       = 2;
  const unsigned int fe_degree = 1;

  LaplaceSolver::LaplaceSolver<dim, fe_degree> solver;
  solver.initialize(2);
  solver.solve();
  solver.output_results(0);

  return 0;
}
