// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * Tests p-multigrid with a block system of two Laplacians, implemented through
 * a matrix-free operator. Considers a hexahedral cubic mesh.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

namespace
{
  using namespace dealii;

  // Physical dimension.
  static constexpr int dim = 3;

  // Number of blocks.
  static constexpr unsigned int n_blocks = 2;

  // Alias for block vectors.
  using BlockVector = LinearAlgebra::distributed::BlockVector<double>;

  // Alias for non-block vectors.
  using Vector = BlockVector::BlockType;

  // Alias for FEEvaluation.
  using FEEval = FEEvaluation<dim, -1, 1, 1, double>;

  // Block-wise Laplace operator.
  class Operator : public MatrixFreeOperators::Base<dim, BlockVector>
  {
  public:
    // Constructor.
    Operator() = default;

    // Compute the diagonal.
    virtual void
    compute_diagonal() override
    {
      inverse_diagonal_entries.reset(new DiagonalMatrix<BlockVector>());
      BlockVector &inverse_diagonal =
        this->inverse_diagonal_entries->get_vector();

      inverse_diagonal.reinit(n_blocks);
      for (unsigned int b = 0; b < n_blocks; ++b)
        data->initialize_dof_vector(inverse_diagonal.block(b), b);
      inverse_diagonal.collect_sizes();

      std::vector<unsigned int> dofs_per_cell(n_blocks);
      std::vector<AlignedVector<VectorizedArray<double>>> local_diagonal(
        n_blocks);
      for (unsigned int b = 0; b < n_blocks; ++b)
        {
          dofs_per_cell[b] = data->get_dofs_per_cell(b);
          local_diagonal[b].resize(dofs_per_cell[b]);
        }

      FEEval phi_0(*data, 0);
      FEEval phi_1(*data, 1);

      auto cell_operation =
        [&](const MatrixFree<dim, double> & /*data*/,
            BlockVector &dst,
            const unsigned int & /*dummy*/,
            const std::pair<unsigned int, unsigned int> &cell_range) {
          for (unsigned int c = cell_range.first; c < cell_range.second; ++c)
            {
              phi_0.reinit(c);

              for (unsigned int i = 0; i < dofs_per_cell[0]; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell[0]; ++j)
                    phi_0.submit_dof_value(i == j ? 1.0 : 0.0, j);

                  phi_0.evaluate(EvaluationFlags::gradients);
                  for (unsigned int q = 0; q < phi_0.n_q_points; ++q)
                    phi_0.submit_gradient(phi_0.get_gradient(q), q);
                  phi_0.integrate(EvaluationFlags::gradients);

                  local_diagonal[0][i] = phi_0.get_dof_value(i);
                }

              for (unsigned int i = 0; i < dofs_per_cell[0]; ++i)
                phi_0.submit_dof_value(local_diagonal[0][i], i);
              phi_0.distribute_local_to_global(dst.block(0));

              phi_1.reinit(c);

              for (unsigned int i = 0; i < dofs_per_cell[1]; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell[1]; ++j)
                    phi_1.submit_dof_value(i == j ? 1.0 : 0.0, j);

                  phi_1.evaluate(EvaluationFlags::gradients);
                  for (unsigned int q = 0; q < phi_1.n_q_points; ++q)
                    phi_1.submit_gradient(phi_1.get_gradient(q), q);
                  phi_1.integrate(EvaluationFlags::gradients);

                  local_diagonal[1][i] = phi_1.get_dof_value(i);
                }

              for (unsigned int i = 0; i < dofs_per_cell[1]; ++i)
                phi_1.submit_dof_value(local_diagonal[1][i], i);
              phi_1.distribute_local_to_global(dst.block(1));
            }
        };

      unsigned int dummy = 0;
      data->cell_loop<BlockVector, unsigned int>(cell_operation,
                                                 inverse_diagonal,
                                                 dummy);

      set_constrained_entries_to_one(inverse_diagonal);

      for (const auto &idx : inverse_diagonal.locally_owned_elements())
        inverse_diagonal[idx] = 1.0 / inverse_diagonal[idx];
      inverse_diagonal.compress(VectorOperation::insert);
    }

  protected:
    // Compute the matrix-vector product.
    virtual void
    apply_add(BlockVector &dst, const BlockVector &src) const override
    {
      FEEval phi_0(*data, 0);
      FEEval phi_1(*data, 1);

      auto local_apply =
        [this,
         &phi_0,
         &phi_1](const MatrixFree<dim, double> & /*data*/,
                 BlockVector                                 &dst,
                 const BlockVector                           &src,
                 const std::pair<unsigned int, unsigned int> &cell_range) {
          for (unsigned int cell = cell_range.first; cell < cell_range.second;
               ++cell)
            {
              {
                phi_0.reinit(cell);
                phi_0.read_dof_values(src.block(0));
                phi_0.evaluate(EvaluationFlags::gradients);

                for (unsigned int q = 0; q < phi_0.n_q_points; ++q)
                  phi_0.submit_gradient(phi_0.get_gradient(q), q);

                phi_0.integrate(EvaluationFlags::gradients);
                phi_0.distribute_local_to_global(dst.block(0));
              }

              {
                phi_1.reinit(cell);
                phi_1.read_dof_values(src.block(1));
                phi_1.evaluate(EvaluationFlags::gradients);

                for (unsigned int q = 0; q < phi_1.n_q_points; ++q)
                  phi_1.submit_gradient(phi_1.get_gradient(q), q);

                phi_1.integrate(EvaluationFlags::gradients);
                phi_1.distribute_local_to_global(dst.block(1));
              }
            }
        };

      data->cell_loop<BlockVector, BlockVector>(local_apply, dst, src);
    }
  };

  // Alias for the multigrid smoother type.
  using SmootherType =
    PreconditionChebyshev<Operator, BlockVector, DiagonalMatrix<BlockVector>>;

  // Custom MGTransfer class for block vectors.
  //
  // This class has very similar functionality to that of
  // MGTransferBlockMatrixFree, except that it initializes the underlying
  // transfers for the single blocks using MGTwoLevelTransfer (whereas
  // MGTransferBlockMatrixFree only allows to build based on the constraints
  // and DoF handlers, i.e. does not expose
  // MGTransferMatrixFree::initialize_two_level_transfers with a possibly
  // different transfer for each block).
  class MGTransferBlockCustom
    : public MGTransferBlockMatrixFreeBase<dim,
                                           double,
                                           MGTransferMatrixFree<dim, double>>
  {
  public:
    // This can be used for block vectors.
    static const bool supports_dof_handler_vector = true;

    // Alias for the multigrid transfer type on each block.
    using TransferType = MGTransferMatrixFree<dim, double>;

    // Constructor.
    MGTransferBlockCustom(const unsigned int &n_blocks_)
      : MGTransferBlockMatrixFreeBase(/* same_for_all = */ false)
      , n_blocks(n_blocks_)
      , block_transfers(n_blocks)
      , two_level_transfers(n_blocks)
    {}

    // Reinitialize based on the provided operators.
    //
    // The operator type is assumed to expose a method get_matrix_free(), which
    // is the case as long as it is derived from MatrixFreeOperators::Base.
    template <class OperatorType>
    void
    reinit(const MGLevelObject<std::unique_ptr<OperatorType>> &operators)
    {
      const unsigned int min_level = operators.min_level();
      const unsigned int max_level = operators.max_level();

      for (unsigned int b = 0; b < n_blocks; ++b)
        {
          two_level_transfers[b].resize(min_level, max_level);

          for (unsigned int level = min_level + 1; level <= max_level; ++level)
            {
              two_level_transfers[b][level].reinit(
                *operators[level]->get_matrix_free(),
                b,
                *operators[level - 1]->get_matrix_free(),
                b);
            }

          block_transfers[b].initialize_two_level_transfers(
            two_level_transfers[b]);

          // Here we assume that each block only has one component. If that were
          // not the case, we'd need to call
          // block_transfers[b].set_component_to_block_map(...).

          block_transfers[b].build();
        }
    }

  protected:
    // Get the transfer operator for a given block.
    virtual const TransferType &
    get_matrix_free_transfer(const unsigned int b) const override
    {
      return block_transfers[b];
    }

    // Number of blocks.
    const unsigned int n_blocks;

    // Vector of transfer operators, one for each block.
    std::vector<TransferType> block_transfers;

    // Vectors of two level transfers. The outer vector runs over the
    // blocks, the inner MGLevelObject runs over the levels.
    std::vector<MGLevelObject<MGTwoLevelTransfer<dim, Vector>>>
      two_level_transfers;
  };

  // p-multigrid test class.
  class PMGTest
  {
  public:
    // Constructor.
    PMGTest(const unsigned int &n_refs_, const unsigned int &fe_degree_)
      : n_refs(n_refs_)
      , fe_degree(fe_degree_)
      , mg_min_level(0)
      , mg_max_level(fe_degree - 1)
      , mg_n_levels(mg_max_level - mg_min_level + 1)
    {}

    // Run the example.
    void
    run()
    {
      const std::string separator_section(50, '=');
      const std::string separator_subsection(50, '-');

      // Create the mesh.
      {
        deallog << separator_section << std::endl;
        deallog << "Generating mesh..." << std::flush;

        GridGenerator::hyper_cube(triangulation,
                                  0,
                                  1,
                                  /* colorize = */ true);
        triangulation.refine_global(n_refs);

        deallog << " done." << std::endl;
        deallog << "Number of active cells: " << triangulation.n_active_cells()
                << std::endl;
      }

      // Setup p-multigrid.
      {
        deallog << separator_section << std::endl;
        deallog << "Initializing p-multigrid hierarchy:" << std::endl;

        const unsigned int min_level = 0;
        const unsigned int max_level = fe_degree - 1;
        const unsigned int n_levels  = max_level - min_level + 1;

        deallog << "  Number of levels: " << n_levels << std::endl;

        mg_dof_handlers.resize(n_levels);
        mg_constraints.resize(n_levels);
        mg_operators.resize(min_level, max_level);

        // Boundary functions.
        std::map<types::boundary_id, const Function<dim> *> boundary_functions;
        for (unsigned int i = 0; i < 2 * dim; ++i)
          boundary_functions[i] = &zero_function;

        for (unsigned int l = mg_min_level; l <= mg_max_level; ++l)
          {
            deallog << separator_subsection << std::endl;
            deallog << "  Level " << l << ":" << std::endl;

            const unsigned int degree = l + 1;
            deallog << "    Degree: " << degree << std::endl;

            // Construct the DoF handlers and constraints, one for each
            // block.
            mg_dof_handlers[l].resize(n_blocks);
            mg_constraints[l].resize(n_blocks);

            for (unsigned int b = 0; b < n_blocks; ++b)
              {
                const FE_Q<dim> fe(degree);

                mg_dof_handlers[l][b] = std::make_unique<DoFHandler<dim>>();
                mg_dof_handlers[l][b]->reinit(triangulation);
                mg_dof_handlers[l][b]->distribute_dofs(fe);

                deallog << "    Block " << b << ": "
                        << mg_dof_handlers[l][b]->n_dofs() << " DoFs"
                        << std::endl;

                mg_constraints[l][b].reinit(
                  mg_dof_handlers[l][b]->locally_owned_dofs(),
                  DoFTools::extract_locally_relevant_dofs(
                    *mg_dof_handlers[l][b]));
                VectorTools::interpolate_boundary_values(mapping,
                                                         *mg_dof_handlers[l][b],
                                                         boundary_functions,
                                                         mg_constraints[l][b]);
                mg_constraints[l][b].close();
              }

            // Setup the MatrixFree object for this level.
            typename MatrixFree<dim, double>::AdditionalData additional_data;
            additional_data.tasks_parallel_scheme =
              MatrixFree<dim, double>::AdditionalData::none;
            additional_data.mapping_update_flags =
              update_gradients | update_JxW_values | update_values |
              update_quadrature_points;

            std::vector<const DoFHandler<dim> *>           dof_handler_ptrs;
            std::vector<const AffineConstraints<double> *> constraint_ptrs;
            for (unsigned int i = 0; i < n_blocks; ++i)
              {
                dof_handler_ptrs.push_back(mg_dof_handlers[l][i].get());
                constraint_ptrs.push_back(&mg_constraints[l][i]);
              }

            const QGauss<dim> quadrature(degree + 1);
            auto mf_level = std::make_shared<MatrixFree<dim, double>>();
            mf_level->reinit(mapping,
                             dof_handler_ptrs,
                             constraint_ptrs,
                             quadrature,
                             additional_data);

            // Construct the operator for this level.
            mg_operators[l] = std::make_unique<Operator>();
            mg_operators[l]->initialize(mf_level);
          }

        mg_matrix.initialize(mg_operators);

        // Setup the transfer operators.
        {
          deallog << separator_subsection << std::endl;
          deallog << "  Initializing transfer operators" << std::endl;

          mg_transfer = std::make_unique<MGTransferBlockCustom>(n_blocks);
          mg_transfer->reinit(mg_operators);
        }

        // Setup the smoothers.
        {
          deallog << "  Initializing smoothers" << std::endl;

          MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
            mg_min_level, mg_max_level);

          for (unsigned int l = mg_min_level; l <= mg_max_level; ++l)
            {
              if (l > 0)
                {
                  smoother_data[l].smoothing_range     = 15.;
                  smoother_data[l].degree              = 5;
                  smoother_data[l].eig_cg_n_iterations = 10;
                }
              else
                {
                  smoother_data[l].smoothing_range = 1e-3;
                  smoother_data[l].degree = numbers::invalid_unsigned_int;
                  smoother_data[l].eig_cg_n_iterations = mg_operators[0]->m();
                }

              mg_operators[l]->compute_diagonal();
              smoother_data[l].preconditioner =
                mg_operators[l]->get_matrix_diagonal_inverse();
            }

          mg_smoother.initialize(mg_operators, smoother_data);
        }

        // Setup the coarse solver.
        {
          deallog << "  Initializing coarse solver" << std::endl;

          mg_coarse.initialize(mg_smoother);
        }

        // Build the multigrid preconditioner object.
        {
          deallog << "  Initializing multigrid preconditioner" << std::endl;

          mg = std::make_unique<Multigrid<BlockVector>>(
            mg_matrix, mg_coarse, *mg_transfer, mg_smoother, mg_smoother);

          std::vector<const DoFHandler<dim> *> dof_handler_ptrs(n_blocks);
          for (unsigned int i = 0; i < n_blocks; ++i)
            dof_handler_ptrs[i] = mg_dof_handlers[mg_max_level][i].get();

          mg_preconditioner = std::make_unique<
            PreconditionMG<dim, BlockVector, MGTransferBlockCustom>>(
            dof_handler_ptrs, *mg, *mg_transfer);
        }
      }

      // Apply the multigrid preconditioner to a vector.
      {
        deallog << separator_section << std::endl;
        deallog << "Applying multigrid preconditioner" << std::endl;

        BlockVector x(n_blocks);
        for (unsigned int b = 0; b < n_blocks; ++b)
          mg_operators[mg_max_level]->get_matrix_free()->initialize_dof_vector(
            x.block(b), b);
        x.collect_sizes();
        x(1) = 1.0;

        BlockVector y(n_blocks);
        y.reinit(x);

        mg_preconditioner->vmult(y, x);

        deallog << "  ||x||_2 = " << x.l2_norm() << std::endl;
        deallog << "  ||y||_2 = " << y.l2_norm() << std::endl;
      }
    }

  protected:
    // Number of mesh refinements.
    const unsigned int n_refs;

    // Finite element degree.
    const unsigned int fe_degree;

    // Triangulation.
    Triangulation<dim> triangulation;

    // Minimum multigrid level.
    const unsigned int mg_min_level;

    // Maximum multigrid level.
    const unsigned int mg_max_level;

    // Number of multigrid levels.
    const unsigned int mg_n_levels;

    // Mapping.
    MappingQ1<dim> mapping;

    // DoF handlers for the multigrid levels. The outer vector runs over the
    // levels, the inner vector runs over the blocks.
    std::vector<std::vector<std::unique_ptr<DoFHandler<dim>>>> mg_dof_handlers;

    // Constraints for the multigrid levels. The outer vector runs over the
    // levels, the inner vector runs over the blocks.
    std::vector<std::vector<AffineConstraints<double>>> mg_constraints;

    // Operators on the multigrid levels.
    MGLevelObject<std::unique_ptr<Operator>> mg_operators;

    // Transfer between multigrid levels.
    std::unique_ptr<MGTransferBlockCustom> mg_transfer;

    // Smoothers for the multigrid preconditioner.
    mg::SmootherRelaxation<SmootherType, BlockVector> mg_smoother;

    // Multigrid matrices.
    mg::Matrix<BlockVector> mg_matrix;

    // Coarse solver.
    MGCoarseGridApplySmoother<BlockVector> mg_coarse;

    // Multigrid.
    std::unique_ptr<Multigrid<BlockVector>> mg;

    // Multigrid preconditioner.
    std::unique_ptr<PreconditionMG<dim, BlockVector, MGTransferBlockCustom>>
      mg_preconditioner;

    // Zero function, used to set Dirichlet BCs.
    Functions::ZeroFunction<dim> zero_function;
  };
} // namespace

int
main(int /* argc */, char * /* argv */[])
{
  initlog();

  PMGTest example(/* n_refs = */ 2, /* fe_degree = */ 2);
  example.run();

  return 0;
}
