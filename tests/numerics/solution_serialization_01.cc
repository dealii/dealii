// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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



// This test is used to make sure that FESeries::Fourier/Legendre
// work with non-primitive elements

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/solution_serialization.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    (void)component;

    return p[0];
  }
};

template <int dim>
void
create_mesh(Triangulation<dim> &                         tria,
            const std::vector<BoundingBox<dim, double>> &bbs)
{
  GridGenerator::subdivided_hyper_cube(tria, 2);

  for (unsigned int i = 0; i < bbs.size(); ++i)
    {
      for (const auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            if (bbs[i].point_inside(cell->center()))
              cell->set_refine_flag();
          }

      tria.execute_coarsening_and_refinement();
    }
}

template <int dim>
void
print_mesh(Triangulation<dim> &tria)
{
  for (unsigned int l = 0; l < tria.n_global_levels(); ++l)
    {
      deallog << "level " << std::to_string(l) << ":" << std::endl;
      for (const auto cell : tria.cell_iterators_on_level(l))
        deallog << " " << cell->center() << std::endl;

      deallog << std::endl << std::endl;
    }

  deallog << "active level:" << std::endl;
  for (const auto cell : tria.active_cell_iterators())
    deallog << " " << cell->center() << std::endl;

  deallog << std::endl << std::endl;
}


template <int dim>
void
test()
{
  const unsigned int no_refinement = 2;

  // 1) create triangulation by refining a quadrant over and over;
  // this leads to refining cells on coarser levels

  // ... specify the quadrant via a bounding box
  std::pair<Point<dim, double>, Point<dim, double>> points;
  points.first  = {0.0, 0.5};
  points.second = {0.5, 1.0};

  std::vector<BoundingBox<dim, double>> bbs1(no_refinement,
                                             BoundingBox<dim, double>(points));

  // ... perform refinement
  Triangulation<dim> tria1(
    Triangulation<dim>::MeshSmoothing::limit_level_difference_at_vertices);
  create_mesh(tria1, bbs1);

  // ... remove last level
  for (const auto cell : tria1.cell_iterators_on_level(bbs1.size()))
    cell->set_coarsen_flag();
  tria1.execute_coarsening_and_refinement();
  print_mesh(tria1);

  // 2) create triangulation by refining only cells on the
  // finest level

  // ... specify the cells to refine on each level by a boundig box
  std::vector<BoundingBox<dim, double>> bbs2;

  for (unsigned int l = 1; l < tria1.n_global_levels(); ++l)
    {
      std::vector<Point<dim, double>> points;

      for (const auto cell : tria1.cell_iterators_on_level(l))
        for (const auto v : cell->vertex_indices())
          points.push_back(cell->vertex(v));

      bbs2.emplace_back(points);
    }

  // ... run test
  Triangulation<dim> tria2(
    Triangulation<dim>::MeshSmoothing::limit_level_difference_at_vertices);
  create_mesh(tria2, bbs2);
  print_mesh(tria2);

  MappingQ1<dim> mapping;
  FE_Q<dim>      fe(2);
  QGauss<dim>    quad(3);

  DoFHandler<dim> dof_handler_1(tria1);
  dof_handler_1.distribute_dofs(fe);

  DoFHandler<dim> dof_handler_2(tria2);
  dof_handler_2.distribute_dofs(fe);

  using VectorType = LinearAlgebra::distributed::Vector<double>;

  VectorType vec1(dof_handler_1.n_dofs());
  VectorType vec2(dof_handler_2.n_dofs());

  AffineConstraints<double> dummy;
  dummy.close();

  VectorTools::project(
    mapping, dof_handler_1, dummy, quad, RightHandSideFunction<dim>(), vec1);

  SolutionSerialization<dim, VectorType> ss1(dof_handler_1);
  ss1.add_vector(vec1);
  ss1.save("temp");

  SolutionSerialization<dim, VectorType> ss2(dof_handler_2);
  ss2.add_vector(vec2);
  ss2.load("temp");

  Vector<float> norm_per_cell(tria2.n_active_cells());
  VectorTools::integrate_difference(dof_handler_2,
                                    vec2,
                                    RightHandSideFunction<dim>(),
                                    norm_per_cell,
                                    quad,
                                    VectorTools::L2_norm);
  const double error_L2_norm =
    VectorTools::compute_global_error(tria2,
                                      norm_per_cell,
                                      VectorTools::L2_norm);

  if (error_L2_norm < 1e-10)
    deallog << "OK!" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  test<2>();
}
