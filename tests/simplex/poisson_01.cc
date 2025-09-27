// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Solve Poisson problem on a tet mesh and on a quad mesh with the same number
// of subdivisions.

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

#include "simplex_grids.h"


template <int dim>
struct Parameters
{
  unsigned int degree = 2;


  // GridGenerator
  bool                      use_grid_generator = true;
  std::vector<unsigned int> repetitions;
  Point<dim>                p1;
  Point<dim>                p2;

  // GridIn
  std::string file_name_in = "";

  // GridOut
  std::string file_name_out = "";
};

template <int dim, int spacedim>
MPI_Comm
get_communicator(const Triangulation<dim, spacedim> &tria)
{
  if (auto tria_ =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(&tria))
    return tria_->get_mpi_communicator();

  return MPI_COMM_SELF;
}

template <int dim, int spacedim = dim>
void
test(const Triangulation<dim, spacedim> &tria,
     const FiniteElement<dim, spacedim> &fe,
     const Quadrature<dim>              &quad,
     const hp::QCollection<dim - 1>     &face_quad,
     const Mapping<dim, spacedim>       &mapping,
     const double                        r_boundary,
     const bool                          do_use_fe_face_values = true)
{
  std::string label =
    (dynamic_cast<const parallel::shared::Triangulation<dim, spacedim> *>(
       &tria) ?
       "parallel::shared::Triangulation" :
       (dynamic_cast<
          const parallel::fullydistributed::Triangulation<dim, spacedim> *>(
          &tria) ?
          "parallel::fullydistributed::Triangulation" :
          (dynamic_cast<
             const parallel::distributed::Triangulation<dim, spacedim> *>(
             &tria) ?
             "parallel::distributed::Triangulation" :
             "Triangulation")));

  deallog << "   on " << label << std::endl;


  for (const auto &cell : tria.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary() &&
          (std::abs(face->center()[0] - r_boundary) < 1e-6))
        face->set_boundary_id(1);
      else if (face->at_boundary() && face->center()[1] == 0.0)
        face->set_boundary_id(2);
      else if (face->at_boundary() && face->center()[1] == 1.0)
        face->set_boundary_id(2);
      else if (dim == 3 && face->at_boundary() && face->center()[2] == 0.0)
        face->set_boundary_id(2);
      else if (dim == 3 && face->at_boundary() && face->center()[2] == 1.0)
        face->set_boundary_id(2);
      else if (face->at_boundary())
        face->set_boundary_id(0);


  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraint_matrix;
  DoFTools::make_zero_boundary_constraints(dof_handler, 0, constraint_matrix);
  constraint_matrix.close();

  // constraint_matrix.print(std::cout);

  const MPI_Comm comm = get_communicator(dof_handler.get_triangulation());

  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);


#ifdef DEAL_II_WITH_TRILINOS
  using VectorType = LinearAlgebra::distributed::Vector<double>;
  TrilinosWrappers::SparseMatrix system_matrix;
  VectorType                     solution;
  VectorType                     system_rhs;
#else
  using VectorType = Vector<double>;
  SparsityPattern        sparsity_pattern;
  SparseMatrix<double>   system_matrix;
  VectorType             solution(dof_handler.n_dofs());
  VectorType             system_rhs(dof_handler.n_dofs());
#endif


#ifdef DEAL_II_WITH_TRILINOS
  TrilinosWrappers::SparsityPattern dsp(dof_handler.locally_owned_dofs(), comm);
#else
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
#endif
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraint_matrix);
#ifdef DEAL_II_WITH_TRILINOS
  dsp.compress();
  system_matrix.reinit(dsp);


  solution.reinit(dof_handler.locally_owned_dofs(),
                  locally_relevant_dofs,
                  comm);
  system_rhs.reinit(dof_handler.locally_owned_dofs(),
                    locally_relevant_dofs,
                    comm);
#else
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
#endif

  const UpdateFlags flag = update_JxW_values | update_values |
                           update_gradients | update_quadrature_points;
  FEValues<dim, spacedim> fe_values(mapping, fe, quad, flag);

  std::shared_ptr<FEFaceValues<dim, spacedim>> fe_face_values;

  if (do_use_fe_face_values)
    fe_face_values.reset(
      new FEFaceValues<dim, spacedim>(mapping, fe, face_quad, flag));

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quad.size();

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  for (const auto &cell : dof_handler.cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            1.0 *                               // 1.0
                            fe_values.JxW(q_index));            // dx
          }

      if (fe_face_values)
        for (const auto &face : cell->face_iterators())
          if (face->at_boundary() && (face->boundary_id() == 1))
            {
              fe_face_values->reinit(cell, face);
              for (const auto q : fe_face_values->quadrature_point_indices())
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  cell_rhs(i) +=
                    (1.0 *                               // 1.0
                     fe_face_values->shape_value(i, q) * // phi_i(x_q)
                     fe_face_values->JxW(q));            // dx
            }

      cell->get_dof_indices(dof_indices);

      constraint_matrix.distribute_local_to_global(
        cell_matrix, cell_rhs, dof_indices, system_matrix, system_rhs);
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  SolverControl        solver_control(1000, 1e-12);
  SolverCG<VectorType> solver(solver_control);

  const std::vector<ReferenceCell> reference_cells = tria.get_reference_cells();
  Assert(reference_cells.size() == 1, ExcNotImplemented());
  unsigned int lower = 0;
  unsigned int upper = 0;
  if (reference_cells[0] == ReferenceCells::Triangle)
    {
      lower = 111;
      upper = 115;
    }
  else if (reference_cells[0] == ReferenceCells::Quadrilateral)
    {
      lower = 96;
      upper = 100;
    }
  else if (reference_cells[0] == ReferenceCells::Tetrahedron)
    {
      lower = 154;
      upper = 158;
    }
  else if (reference_cells[0] == ReferenceCells::Hexahedron)
    {
      lower = 132;
      upper = 136;
    }
  else if (reference_cells[0] == ReferenceCells::Wedge)
    {
      lower = 194;
      upper = 198;
    }
  else if (reference_cells[0] == ReferenceCells::Pyramid)
    {
      lower = 81;
      upper = 85;
    }
  else
    Assert(false, ExcInternalError());

  check_solver_within_range(
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()),
    solver_control.last_step(),
    lower,
    upper);

  // system_rhs.print(std::cout);
  // solution.print(std::cout);

  bool hex_mesh = true;

  for (const auto &cell : tria.active_cell_iterators())
    hex_mesh &= (cell->n_vertices() == GeometryInfo<dim>::vertices_per_cell);

  deallog << std::endl;


  if (false)
    {
      Vector<double> difference(tria.n_active_cells());

      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        Functions::ZeroFunction<dim>(),
                                        difference,
                                        quad,
                                        VectorTools::L2_norm);

      deallog << VectorTools::compute_global_error(tria,
                                                   difference,
                                                   VectorTools::L2_norm)
              << std::endl;
      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");

      data_out.build_patches(mapping);

      static unsigned int counter = 0;

      std::ofstream output("result" + std::to_string(counter++) + ".vtk");
      data_out.write_vtk(output);
    }
}

template <int dim, int spacedim = dim>
void
test_tet(const MPI_Comm comm, const Parameters<dim> &params)
{
  const unsigned int tria_type = 2;

  // 1) Create triangulation...
  Triangulation<dim, spacedim> *tria;

  // a) serial triangulation
  Triangulation<dim, spacedim> tr_1;

  // b) shared triangulation (with artificial cells)
  parallel::shared::Triangulation<dim> tr_2(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_custom_signal);

  tr_2.signals.create.connect([&]() {
    GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(comm),
                                       tr_2);
  });

  // c) distributed triangulation
  parallel::fullydistributed::Triangulation<dim> tr_3(comm);


  // ... choose the right triangulation
  if (tria_type == 0 || tria_type == 2)
    tria = &tr_1;
  else if (tria_type == 1)
    tria = &tr_2;

  // ... create triangulation
  if (params.use_grid_generator)
    {
      // ...via GridGenerator
      GridGenerator::subdivided_hyper_rectangle_with_simplices(
        *tria, params.repetitions, params.p1, params.p2, false);
    }
  else
    {
      // ...via GridIn
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(*tria);
      std::ifstream input_file(params.file_name_in);
      grid_in.read_ucd(input_file);
      // std::ifstream input_file("test_tet_geometry.unv");
      // grid_in.read_unv(input_file);
    }

  // ... partition serial triangulation and create distributed triangulation
  if (tria_type == 0 || tria_type == 2)
    {
      GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(comm),
                                         tr_1);

      auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(tr_1, comm);

      tr_3.create_triangulation(construction_data);

      tria = &tr_3;
    }

  // 2) Output generated triangulation via GridOut
  GridOut       grid_out;
  std::ofstream out(params.file_name_out + "." +
                    std::to_string(Utilities::MPI::this_mpi_process(comm)) +
                    ".vtk");
  grid_out.write_vtk(*tria, out);

  // 3) Select components
  FE_SimplexP<dim> fe(params.degree);

  QGaussSimplex<dim> quad(params.degree + 1);

  hp::QCollection<dim - 1> face_quad{QGaussSimplex<dim - 1>(params.degree + 1)};

  FE_SimplexP<dim> fe_mapping(1);
  MappingFE<dim>   mapping(fe_mapping);

  // 4) Perform test (independent of mesh type)
  test(*tria, fe, quad, face_quad, mapping, params.p2[0]);
}

template <int dim, int spacedim = dim>
void
test_hex(const MPI_Comm comm, const Parameters<dim> &params)
{
  // 1) Create triangulation...
  parallel::distributed::Triangulation<dim, spacedim> tria(comm);

  if (params.use_grid_generator)
    {
      // ...via GridGenerator
      GridGenerator::subdivided_hyper_rectangle(
        tria, params.repetitions, params.p1, params.p2, false);
    }
  else
    {
      // ...via GridIn
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(tria);
      std::ifstream input_file(params.file_name_in);
      grid_in.read_ucd(input_file);
    }

  // 2) Output generated triangulation via GridOut
  GridOut       grid_out;
  std::ofstream out(params.file_name_out + "." +
                    std::to_string(Utilities::MPI::this_mpi_process(comm)) +
                    ".vtk");
  grid_out.write_vtk(tria, out);

  // 3) Select components
  FE_Q<dim> fe(params.degree);

  QGauss<dim> quad(params.degree + 1);

  hp::QCollection<dim - 1> face_quad{QGauss<dim - 1>(params.degree + 1)};

  MappingQ<dim, spacedim> mapping(1);

  // 4) Perform test (independent of mesh type)
  test(tria, fe, quad, face_quad, mapping, params.p2[0]);
}

template <int dim, int spacedim = dim>
void
test_wedge(const MPI_Comm comm, const Parameters<dim> &params)
{
  const unsigned int tria_type = 2;

  // 1) Create triangulation...
  Triangulation<dim, spacedim> *tria;

  // a) serial triangulation
  Triangulation<dim, spacedim> tr_1;

  // b) shared triangulation (with artificial cells)
  parallel::shared::Triangulation<dim> tr_2(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_custom_signal);

  tr_2.signals.create.connect([&]() {
    GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(comm),
                                       tr_2);
  });

  // c) distributed triangulation
  parallel::fullydistributed::Triangulation<dim> tr_3(comm);


  // ... choose the right triangulation
  if (tria_type == 0 || tria_type == 2)
    tria = &tr_1;
  else if (tria_type == 1)
    tria = &tr_2;

  // ... create triangulation
  if (params.use_grid_generator)
    {
      // ...via GridGenerator
      GridGenerator::subdivided_hyper_rectangle_with_wedges(
        *tria, params.repetitions, params.p1, params.p2, false);
    }
  else
    {
      // ...via GridIn
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(*tria);
      std::ifstream input_file(params.file_name_in);
      grid_in.read_ucd(input_file);
      // std::ifstream input_file("test_tet_geometry.unv");
      // grid_in.read_unv(input_file);
    }

  // ... partition serial triangulation and create distributed triangulation
  if (tria_type == 0 || tria_type == 2)
    {
      GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(comm),
                                         tr_1);

      auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(tr_1, comm);

      tr_3.create_triangulation(construction_data);

      tria = &tr_3;
    }

  // 2) Output generated triangulation via GridOut
  GridOut       grid_out;
  std::ofstream out(params.file_name_out + "." +
                    std::to_string(Utilities::MPI::this_mpi_process(comm)) +
                    ".vtk");
  grid_out.write_vtk(*tria, out);

  // 3) Select components
  FE_WedgeP<dim> fe(params.degree);

  QGaussWedge<dim> quad(params.degree + 1);

  hp::QCollection<dim - 1> face_quad{QGaussSimplex<dim - 1>(params.degree + 1),
                                     QGaussSimplex<dim - 1>(params.degree + 1),
                                     QGauss<dim - 1>(params.degree + 1),
                                     QGauss<dim - 1>(params.degree + 1),
                                     QGauss<dim - 1>(params.degree + 1)};

  FE_WedgeP<dim> fe_mapping(1);
  MappingFE<dim> mapping(fe_mapping);

  // 4) Perform test (independent of mesh type)
  test(*tria, fe, quad, face_quad, mapping, params.p2[0], true);
}

template <int dim, int spacedim = dim>
void
test_pyramid(const MPI_Comm comm, const Parameters<dim> &params)
{
  const unsigned int tria_type = 2;

  // 1) Create triangulation...
  Triangulation<dim, spacedim> *tria;

  // a) serial triangulation
  Triangulation<dim, spacedim> tr_1;

  // b) shared triangulation (with artificial cells)
  parallel::shared::Triangulation<dim> tr_2(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_custom_signal);

  tr_2.signals.create.connect([&]() {
    GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(comm),
                                       tr_2);
  });

  // c) distributed triangulation
  parallel::fullydistributed::Triangulation<dim> tr_3(comm);


  // ... choose the right triangulation
  if (tria_type == 0 || tria_type == 2)
    tria = &tr_1;
  else if (tria_type == 1)
    tria = &tr_2;

  // ... create triangulation
  if (params.use_grid_generator)
    {
      // ...via GridGenerator
      GridGenerator::subdivided_hyper_rectangle_with_pyramids(
        *tria, params.repetitions, params.p1, params.p2, false);
    }
  else
    {
      // ...via GridIn
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(*tria);
      std::ifstream input_file(params.file_name_in);
      grid_in.read_ucd(input_file);
      // std::ifstream input_file("test_tet_geometry.unv");
      // grid_in.read_unv(input_file);
    }

  // ... partition serial triangulation and create distributed triangulation
  if (tria_type == 0 || tria_type == 2)
    {
      GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(comm),
                                         tr_1);

      auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(tr_1, comm);

      tr_3.create_triangulation(construction_data);

      tria = &tr_3;
    }

  // 2) Output generated triangulation via GridOut
  GridOut       grid_out;
  std::ofstream out(params.file_name_out + "." +
                    std::to_string(Utilities::MPI::this_mpi_process(comm)) +
                    ".vtk");
  grid_out.write_vtk(*tria, out);

  // 3) Select components
  FE_PyramidP<dim> fe(params.degree);

  QGaussPyramid<dim> quad(params.degree + 1);

  hp::QCollection<dim - 1> face_quad{QGauss<dim - 1>(params.degree + 1),
                                     QGaussSimplex<dim - 1>(params.degree + 1),
                                     QGaussSimplex<dim - 1>(params.degree + 1),
                                     QGaussSimplex<dim - 1>(params.degree + 1),
                                     QGaussSimplex<dim - 1>(params.degree + 1)};

  FE_PyramidP<dim> fe_mapping(1);
  MappingFE<dim>   mapping(fe_mapping);

  // 4) Perform test (independent of mesh type)
  test(*tria, fe, quad, face_quad, mapping, params.p2[0], true);
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  const MPI_Comm comm = MPI_COMM_WORLD;

  // 2D
  {
    Parameters<2> params;
    params.use_grid_generator = true;
    params.repetitions        = std::vector<unsigned int>{10, 10};

    // test TRI
    {
      deallog << "Solve problem on TRI mesh:" << std::endl;

      params.file_name_out = "mesh-tri";
      params.p1            = Point<2>(0, 0);
      params.p2            = Point<2>(1, 1);
      test_tet(comm, params);
    }

    // test QUAD
    {
      deallog << "Solve problem on QUAD mesh:" << std::endl;

      params.file_name_out = "mesh-quad";
      params.p1            = Point<2>(1.1, 0); // shift to the right for
      params.p2            = Point<2>(2.1, 1); // visualization purposes
      test_hex(comm, params);
    }
  }

  // 3D
  {
    Parameters<3> params;
    params.use_grid_generator = true;
    params.repetitions        = std::vector<unsigned int>{10, 10, 10};

    // test TET
    {
      deallog << "Solve problem on TET mesh:" << std::endl;

      params.file_name_out = "mesh-tet";
      params.p1            = Point<3>(0, 0, 0);
      params.p2            = Point<3>(1, 1, 1);
      test_tet(comm, params);
    }

    // test HEX
    {
      deallog << "Solve problem on HEX mesh:" << std::endl;

      params.file_name_out = "mesh-hex";
      params.p1            = Point<3>(1.1, 0, 0);
      params.p2            = Point<3>(2.1, 1, 1);
      test_hex(comm, params);
    }

    // test WEDGE
    {
      deallog << "Solve problem on WEDGE mesh:" << std::endl;

      params.file_name_out = "mesh-wedge";
      params.p1            = Point<3>(2.2, 0, 0);
      params.p2            = Point<3>(3.2, 1, 1);
      test_wedge(comm, params);
    }

    // test PYRAMID
    {
      deallog << "Solve problem on PYRAMID mesh:" << std::endl;

      params.file_name_out = "mesh-pyramid";
      params.repetitions   = std::vector<unsigned int>{10, 10, 10};
      params.p1            = Point<3>(3.3, 0, 0);
      params.p2            = Point<3>(4.3, 1, 1);

      params.degree = 1;

      test_pyramid(comm, params);
    }
  }
}
