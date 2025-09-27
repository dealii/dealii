// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Similar to matrix_vector_hp.cc but in parallel.

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "matrix_vector_mf.h"



template <int dim, typename Number>
class MatrixFreeTestHP
{
public:
  MatrixFreeTestHP(const MatrixFree<dim, Number> &data_in)
    : data(data_in){};

  void
  local_apply(const MatrixFree<dim, Number>                    &data,
              LinearAlgebra::distributed::Vector<Number>       &dst,
              const LinearAlgebra::distributed::Vector<Number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    helmholtz_operator_no_template<dim>(data, dst, src, cell_range);
  }

  void
  vmult(LinearAlgebra::distributed::Vector<Number>       &dst,
        const LinearAlgebra::distributed::Vector<Number> &src) const
  {
    dst = 0;
    data.cell_loop(&MatrixFreeTestHP<dim, Number>::local_apply, this, dst, src);
  };

private:
  const MatrixFree<dim, Number> &data;
};



template <int dim, int fe_degree>
void
test()
{
  if (fe_degree > 1)
    return;

  using number = double;
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_ball(tria);
  tria.refine_global(1);

  // refine a few cells
  for (unsigned int i = 0; i < 11 - 3 * dim; ++i)
    {
      unsigned int counter = 0;
      for (const auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned() && ((counter++) % (7 - i) == 0))
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  const unsigned int max_degree = 9 - 2 * dim;

  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim>  quadrature_collection;
  hp::QCollection<1>    quadrature_collection_mf;

  for (unsigned int deg = 1; deg <= max_degree; ++deg)
    {
      fe_collection.push_back(FE_Q<dim>(QGaussLobatto<1>(deg + 1)));
      quadrature_collection.push_back(QGauss<dim>(deg + 1));
      quadrature_collection_mf.push_back(QGauss<1>(deg + 1));
    }

  DoFHandler<dim> dof(tria);
  // set the active FE index in a random order
  {
    for (const auto &cell : dof.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        const unsigned int fe_index = Testing::rand() % max_degree;
        cell->set_active_fe_index(fe_index);
      }
  }

  // setup DoFs
  dof.distribute_dofs(fe_collection);
  AffineConstraints<double> constraints;
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  TrilinosWrappers::SparsityPattern dsp(dof.locally_owned_dofs(),
                                        MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern(dof, dsp, constraints, false);
  dsp.compress();

  TrilinosWrappers::SparseMatrix system_matrix;
  system_matrix.reinit(dsp);

  // std::cout << "Number of cells: " <<
  // dof.get_triangulation().n_active_cells() << std::endl; std::cout << "Number
  // of degrees of freedom: " << dof.n_dofs() << std::endl; std::cout << "Number
  // of constraints: " << constraints.n_constraints() << std::endl;

  // set up MatrixFree
  MatrixFree<dim, number>                          mf_data;
  typename MatrixFree<dim, number>::AdditionalData data;
  data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
  mf_data.reinit(
    MappingQ1<dim>{}, dof, constraints, quadrature_collection_mf, data);
  MatrixFreeTestHP<dim, number> mf(mf_data);

  // assemble sparse matrix with (\nabla v,
  // \nabla u) + (v, 10 * u)
  {
    hp::FEValues<dim>                    hp_fe_values(fe_collection,
                                   quadrature_collection,
                                   update_values | update_gradients |
                                     update_JxW_values);
    FullMatrix<double>                   cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;

    for (const auto &cell : dof.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_matrix = 0;
        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
             ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) += ((fe_values.shape_grad(i, q_point) *
                                         fe_values.shape_grad(j, q_point) +
                                       10. * fe_values.shape_value(i, q_point) *
                                         fe_values.shape_value(j, q_point)) *
                                      fe_values.JxW(q_point));
            }
        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(cell_matrix,
                                               local_dof_indices,
                                               system_matrix);
      }
  }
  system_matrix.compress(VectorOperation::values::add);

  // fill a right hand side vector with random
  // numbers in unconstrained degrees of freedom
  LinearAlgebra::distributed::Vector<double> src;
  mf_data.initialize_dof_vector(src);


  LinearAlgebra::distributed::Vector<double> result_spmv(src), result_mf(src);

  for (unsigned int i = 0; i < dof.n_locally_owned_dofs(); ++i)
    {
      if (constraints.is_constrained(
            src.get_partitioner()->local_to_global(i)) == false)
        src.local_element(i) = random_value<double>();
    }

  // now perform matrix-vector product and check
  // its correctness
  system_matrix.vmult(result_spmv, src);
  mf.vmult(result_mf, src);

  result_mf -= result_spmv;
  const double diff_norm = result_mf.linfty_norm();
  deallog << "Norm of difference: " << diff_norm << std::endl << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  mpi_initlog();

  {
    deallog.push("2d");
    test<2, 1>();
    test<2, 2>();
    deallog.pop();
    deallog.push("3d");
    test<3, 1>();
    test<3, 2>();
    deallog.pop();
  }
}
