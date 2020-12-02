// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Same as matrix_free_01 but testing mixed meshes (and also pure simplex and
// hypercube mesh as special case of mixed meshs).

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/simplex/fe_lib.h>
#include <deal.II/simplex/grid_generator.h>
#include <deal.II/simplex/quadrature_lib.h>

#include "../tests.h"

using namespace dealii;

namespace dealii
{
  namespace GridGenerator
  {
    template <int dim, int spacedim>
    void
    subdivided_hyper_rectangle_with_simplices_mix(
      Triangulation<dim, spacedim> &   tria,
      const std::vector<unsigned int> &repetitions,
      const Point<dim> &               p1,
      const Point<dim> &               p2,
      const bool                       colorize = false)
    {
      AssertDimension(dim, spacedim);

      AssertThrow(colorize == false, ExcNotImplemented());

      std::vector<Point<spacedim>> vertices;
      std::vector<CellData<dim>>   cells;

      if (dim == 2)
        {
          // determine cell sizes
          const Point<dim> dx((p2[0] - p1[0]) / repetitions[0],
                              (p2[1] - p1[1]) / repetitions[1]);

          // create vertices
          for (unsigned int j = 0; j <= repetitions[1]; ++j)
            for (unsigned int i = 0; i <= repetitions[0]; ++i)
              vertices.push_back(
                Point<spacedim>(p1[0] + dx[0] * i, p1[1] + dx[1] * j));

          // create cells
          for (unsigned int j = 0; j < repetitions[1]; ++j)
            for (unsigned int i = 0; i < repetitions[0]; ++i)
              {
                // create reference QUAD cell
                std::array<unsigned int, 4> quad{{
                  (j + 0) * (repetitions[0] + 1) + i + 0, //
                  (j + 0) * (repetitions[0] + 1) + i + 1, //
                  (j + 1) * (repetitions[0] + 1) + i + 0, //
                  (j + 1) * (repetitions[0] + 1) + i + 1  //
                }};                                       //

                if (j < repetitions[1] / 2 && i < repetitions[0] / 2)
                  {
                    CellData<dim> quad_;
                    quad_.vertices = {quad[0], quad[1], quad[2], quad[3]};
                    cells.push_back(quad_);

                    continue;
                  }

                // TRI cell 0
                {
                  CellData<dim> tri;
                  tri.vertices = {quad[0], quad[1], quad[2]};
                  cells.push_back(tri);
                }

                // TRI cell 1
                {
                  CellData<dim> tri;
                  tri.vertices = {quad[3], quad[2], quad[1]};
                  cells.push_back(tri);
                }
              }
        }
      else
        {
          AssertThrow(colorize == false, ExcNotImplemented());
        }

      // actually create triangulation
      tria.create_triangulation(vertices, cells, SubCellData());
    }


    template <int dim, int spacedim>
    void
    subdivided_hyper_cube_with_simplices_mix(Triangulation<dim, spacedim> &tria,
                                             const unsigned int repetitions,
                                             const double       p1 = 0.0,
                                             const double       p2 = 1.0,
                                             const bool colorize   = false)
    {
      if (dim == 2)
        {
          subdivided_hyper_rectangle_with_simplices_mix(
            tria, {{repetitions, repetitions}}, {p1, p1}, {p2, p2}, colorize);
        }
      else if (dim == 3)
        {
          subdivided_hyper_rectangle_with_simplices_mix(
            tria,
            {{repetitions, repetitions, repetitions}},
            {p1, p1, p1},
            {p2, p2, p2},
            colorize);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented())
        }
    }
  } // namespace GridGenerator
} // namespace dealii

template <int dim>
class PoissonOperator
{
public:
  using VectorType       = LinearAlgebra::distributed::Vector<double>;
  using FECellIntegrator = FEEvaluation<dim, -1, 0, 1, double>;

  PoissonOperator(const MatrixFree<dim, double> &matrix_free,
                  const bool                     do_helmholtz)
    : matrix_free(matrix_free)
    , do_helmholtz(do_helmholtz)
  {}

  void
  initialize_dof_vector(VectorType &vec)
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  rhs(VectorType &vec) const
  {
    const int dummy = 0;

    matrix_free.template cell_loop<VectorType, int>(
      [&](const auto &data, auto &dst, const auto &, const auto range) {
        const auto       i = data.get_cell_active_fe_index(range);
        FECellIntegrator phi(matrix_free, 0, 0, 0, i, i);

        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            phi.reinit(cell);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_value(1.0, q);

            phi.integrate_scatter(true, false, dst);
          }
      },
      vec,
      dummy,
      true);
  }


  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    matrix_free.template cell_loop<VectorType, VectorType>(
      [&](const auto &data, auto &dst, const auto &src, const auto range) {
        const auto       i = data.get_cell_active_fe_index(range);
        FECellIntegrator phi(matrix_free, 0, 0, 0, i, i);

        for (unsigned int cell = range.first; cell < range.second; ++cell)
          {
            phi.reinit(cell);
            phi.gather_evaluate(src, do_helmholtz, true);

            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              {
                if (do_helmholtz)
                  phi.submit_value(phi.get_value(q), q);

                phi.submit_gradient(phi.get_gradient(q), q);
              }

            phi.integrate_scatter(do_helmholtz, true, dst);
          }
      },
      dst,
      src,
      true);
  }

private:
  const MatrixFree<dim, double> &matrix_free;
  const bool                     do_helmholtz;
};

template <int dim>
void
test(const unsigned version, const unsigned int degree, const bool do_helmholtz)
{
  Triangulation<dim> tria;

  const unsigned int subdivisions = dim == 2 ? 25 : 8;

  if (version == 0)
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, subdivisions);
  else if (version == 1)
    GridGenerator::subdivided_hyper_cube(tria, subdivisions);
  else if (version == 2)
    GridGenerator::subdivided_hyper_cube_with_simplices_mix(tria, subdivisions);

  Simplex::FE_P<dim>    fe1(degree);
  FE_Q<dim>             fe2(degree);
  hp::FECollection<dim> fes(fe1, fe2);

  Simplex::QGauss<dim> quad1(degree + 1);
  QGauss<dim>          quad2(degree + 1);
  hp::QCollection<dim> quads(quad1, quad2);

  MappingFE<dim>             mapping1(Simplex::FE_P<dim>(1));
  MappingQ<dim>              mapping2(1);
  hp::MappingCollection<dim> mappings(mapping1, mapping2);

  DoFHandler<dim> dof_handler(tria);

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->reference_cell_type() == ReferenceCell::Type::Tri ||
        cell->reference_cell_type() == ReferenceCell::Type::Tet)
      cell->set_active_fe_index(0);
    else
      cell->set_active_fe_index(1);

  dof_handler.distribute_dofs(fes);

  AffineConstraints<double> constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, 0, constraints);
  constraints.close();

  const auto solve_and_postprocess =
    [&](const auto &poisson_operator,
        auto &      x,
        auto &      b) -> std::tuple<unsigned int, double, double, double> {
    ReductionControl reduction_control(1000, 1e-10, 1e-4);
    SolverCG<typename std::remove_reference<decltype(x)>::type> solver(
      reduction_control);
    solver.solve(poisson_operator, x, b, PreconditionIdentity());

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      printf("Solved in %d iterations.\n", reduction_control.last_step());

    constraints.distribute(x);

#if 0
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    x.update_ghost_values();
    data_out.add_data_vector(dof_handler, x, "solution");
    data_out.build_patches(mappings, 2);
    data_out.write_vtu_with_pvtu_record("./", "result", 0, MPI_COMM_WORLD);
#endif

    Vector<double> difference(tria.n_active_cells());

    VectorTools::integrate_difference(mappings,
                                      dof_handler,
                                      x,
                                      Functions::ZeroFunction<dim>(),
                                      difference,
                                      quads,
                                      VectorTools::NormType::L2_norm);

    return {reduction_control.last_step(),
            reduction_control.last_value(),
            x.linfty_norm(),
            VectorTools::compute_global_error(tria,
                                              difference,
                                              VectorTools::NormType::L2_norm)};
  };

  const auto mf_algo = [&]() {
    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.mapping_update_flags = update_gradients | update_values;

    MatrixFree<dim, double> matrix_free;
    matrix_free.reinit(
      mappings, dof_handler, constraints, quads, additional_data);

    PoissonOperator<dim> poisson_operator(matrix_free, do_helmholtz);

    LinearAlgebra::distributed::Vector<double> x, b;
    poisson_operator.initialize_dof_vector(x);
    poisson_operator.initialize_dof_vector(b);

    poisson_operator.rhs(b);

    return solve_and_postprocess(poisson_operator, x, b);
  };

  const auto mb_algo = [&]() {
    Vector<double> x, b;

    x.reinit(dof_handler.n_dofs());
    b.reinit(dof_handler.n_dofs());

    SparseMatrix<double> A;

    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);

    SparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);
    A.reinit(sparsity_pattern);

    const auto flags = update_values | update_gradients | update_JxW_values;

    hp::FEValues<dim> hp_fe_values(mappings, fes, quads, flags);

    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned() == false)
          continue;

        hp_fe_values.reinit(cell);

        auto &fe_values = hp_fe_values.get_present_fe_values();

        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
        cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_rhs.reinit(dofs_per_cell);

        for (const auto q : fe_values.quadrature_point_indices())
          {
            for (const auto i : fe_values.dof_indices())
              for (const auto j : fe_values.dof_indices())
                cell_matrix(i, j) += (fe_values.shape_grad(i, q) *        //
                                        fe_values.shape_grad(j, q) +      //
                                      static_cast<double>(do_helmholtz) * //
                                        fe_values.shape_value(i, q) *     //
                                        fe_values.shape_value(j, q)) *    //
                                     fe_values.JxW(q);                    //

            for (const unsigned int i : fe_values.dof_indices())
              cell_rhs(i) += (fe_values.shape_value(i, q) * //
                              1. *                          //
                              fe_values.JxW(q));            //
          }

        local_dof_indices.resize(cell->get_fe().n_dofs_per_cell());
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, A, b);
      }

    return solve_and_postprocess(A, x, b);
  };

  const auto compare = [&](const auto result_mf, const auto result_mb) {
    AssertDimension(std::get<0>(result_mf), std::get<0>(result_mb));
    Assert(std::abs(std::get<1>(result_mf) - std::get<1>(result_mb)) < 1e-6,
           ExcNotImplemented());
    Assert(std::abs(std::get<2>(result_mf) - std::get<2>(result_mb)) < 1e-6,
           ExcNotImplemented());
    Assert(std::abs(std::get<3>(result_mf) - std::get<3>(result_mb)) < 1e-6,
           ExcNotImplemented());

    deallog << "mesh=";
    if (version == 0)
      deallog << "P";
    else if (version == 1)
      deallog << "Q";
    else if (version == 2)
      deallog << "M";
    deallog << " : ";

    deallog << "dim=" << dim << " ";
    deallog << "degree=" << degree << " ";
    deallog << "Type=";

    if (do_helmholtz)
      deallog << "Helmholtz";
    else
      deallog << "Possion  ";
    deallog << " : ";

    deallog << "Convergence step " << std::get<0>(result_mf) << " value "
            << std::get<1>(result_mf) << " max " << std::get<2>(result_mf)
            << " norm " << std::get<3>(result_mf) << "." << std::endl;
  };

  compare(mf_algo(), mb_algo());
}


int
main(int argc, char **argv)
{
  initlog();

  deallog.depth_file(1);

  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  for (unsigned int i = 0; i < 3; ++i)
    test<2>(i, /*degree=*/1, /*do_helmholtz*/ false);
  deallog << std::endl;

  for (unsigned int i = 0; i < 3; ++i)
    test<2>(i, /*degree=*/1, /*do_helmholtz*/ true);
  deallog << std::endl;

  for (unsigned int i = 0; i < 3; ++i)
    test<2>(i, /*degree=*/2, /*do_helmholtz*/ false);
  deallog << std::endl;

  for (unsigned int i = 0; i < 3; ++i)
    test<2>(i, /*degree=*/2, /*do_helmholtz*/ true);
  deallog << std::endl;
}
