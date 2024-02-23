// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test FEFieldFunction for parallel computations. Specifically, we should be
// able to evaluate it in owned or ghost cells.

// This used to crash with VectorTools::ExcPointNotAvailableHere() because we
// didn't expect getting a ghost cell.

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class LinearFunction : public Function<dim>
{
public:
  double
  value(const Point<dim> &p, const unsigned int) const
  {
    return p[0] + 2;
  }
};


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  // the mesh needs to be irregular enough to throw off the "cell hint" and
  // GridTools::find_active_cell_around_point() to return ghost cells

  GridGenerator::quarter_hyper_shell(tr, Point<dim>(), 0.5, 1.0, 0, true);
  static const SphericalManifold<dim> spherical_manifold;
  tr.set_manifold(99, spherical_manifold);

  for (typename Triangulation<dim>::active_cell_iterator cell =
         tr.begin_active();
       cell != tr.end();
       ++cell)
    cell->set_all_manifold_ids(99);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tr.begin_active();
       cell != tr.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        cell->face(f)->set_all_manifold_ids(numbers::flat_manifold_id);
  //  static const SphericalManifold<dim> boundary_shell;
  //  tr.set_manifold (0, boundary_shell);
  //  tr.set_manifold (1, boundary_shell);

  tr.refine_global((dim == 2) ? 3 : 1);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  TrilinosWrappers::MPI::Vector interpolated(dofh.locally_owned_dofs(),
                                             MPI_COMM_WORLD);
  VectorTools::interpolate(dofh, LinearFunction<dim>(), interpolated);

  const IndexSet relevant_set = DoFTools::extract_locally_relevant_dofs(dofh);
  TrilinosWrappers::MPI::Vector x_rel(relevant_set, MPI_COMM_WORLD);
  x_rel = interpolated;

  Functions::FEFieldFunction<dim, TrilinosWrappers::MPI::Vector> field_function(
    dofh, x_rel);


  std::vector<Point<dim>> points;

  for (typename Triangulation<dim>::active_cell_iterator cell =
         tr.begin_active();
       cell != tr.end();
       ++cell)
    {
      if (cell->is_artificial())
        continue;

      points.push_back(cell->center());

      // throw off the current cell hint
      field_function.set_active_cell(dofh.end());

      {
        Point<dim> p = cell->center();
        Assert(std::abs(field_function.value(p) - (p[0] + 2.)) < 1e-10,
               ExcInternalError());

        Tensor<1, dim> gradient = field_function.gradient(p);
        Tensor<1, dim> exact_gradient;
        exact_gradient[0] = 1.0;
        Assert((gradient - exact_gradient).norm() < 1e-10, ExcInternalError());

        Assert(std::abs(field_function.laplacian(p)) < 1e-10,
               ExcInternalError());
      }

      if (cell->is_locally_owned())
        {
          // more evil: use a corner so we can end up in a different cell
          Point<dim> p = cell->vertex(0);
          Assert(std::abs(field_function.value(p) - (p[0] + 2.)) < 1e-10,
                 ExcInternalError());

          Tensor<1, dim> gradient = field_function.gradient(p);
          Tensor<1, dim> exact_gradient;
          exact_gradient[0] = 1.0;
          Assert((gradient - exact_gradient).norm() < 1e-10,
                 ExcInternalError());

          Assert(std::abs(field_function.laplacian(p)) < 1e-10,
                 ExcInternalError());
        }
    }

  // also hit the code path for list version (no need to cover gradients/laplace
  // here):
  std::vector<double> values(points.size());
  field_function.value_list(points, values, 0);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;

  tr.reset_manifold(99);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
