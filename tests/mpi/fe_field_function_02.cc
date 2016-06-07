// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// Test FEFieldFunction for parallel computations. Specifically, we should be able
// to evaluate it in owned or ghost cells.

// This used to crash with VectorTools::ExcPointNotAvailableHere() because we
// didn't expect getting a ghost cell.

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>


template <int dim>
class LinearFunction : public Function<dim>
{
public:
  double value (const Point<dim> &p,
                const unsigned int) const
  {
    return p[0] + 2;
  }
};


template<int dim>
void test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  // the mesh needs to be irregular enough to throw off the "cell hint" and
  // GridTools::find_active_cell_around_point() to return ghost cells

  GridGenerator::quarter_hyper_shell (tr,
                                      Point<dim>(),
                                      0.5,
                                      1.0,
                                      0,
                                      true);
  static const SphericalManifold<dim> spherical_manifold;
  tr.set_manifold (99, spherical_manifold);

  for (typename Triangulation<dim>::active_cell_iterator
       cell = tr.begin_active();
       cell != tr.end(); ++cell)
    cell->set_all_manifold_ids (99);
  for (typename Triangulation<dim>::active_cell_iterator
       cell = tr.begin_active();
       cell != tr.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->at_boundary(f))
        cell->face(f)->set_all_manifold_ids (numbers::invalid_manifold_id);
  static const HyperShellBoundary<dim> boundary_shell;
  tr.set_boundary (0, boundary_shell);
  tr.set_boundary (1, boundary_shell);

  tr.refine_global ((dim==2)?3:1);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs (fe);

  TrilinosWrappers::MPI::Vector interpolated(dofh.locally_owned_dofs(),
                                             MPI_COMM_WORLD);
  VectorTools::interpolate (dofh,
                            LinearFunction<dim>(),
                            interpolated);

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs (dofh, relevant_set);
  TrilinosWrappers::MPI::Vector x_rel(relevant_set, MPI_COMM_WORLD);
  x_rel = interpolated;

  Functions::FEFieldFunction<dim,DoFHandler<dim>,TrilinosWrappers::MPI::Vector>
  field_function (dofh, x_rel);

  std::vector<Point< dim > > points;

  for (typename Triangulation<dim>::active_cell_iterator cell=tr.begin_active();
       cell!=tr.end(); ++cell)
    {
      if (cell->is_artificial())
        continue;

      points.push_back(cell->center());

      // throw off the current cell hint
      field_function.set_active_cell(dofh.end());

      {
        Point<dim> p = cell->center();
        Assert(std::abs(field_function.value (p)-(p[0]+2.))<1e-10, ExcInternalError());

        Tensor<1,dim> gradient = field_function.gradient (p);
        Tensor<1,dim> exact_gradient;
        exact_gradient[0]=1.0;
        Assert((gradient-exact_gradient).norm()<1e-10, ExcInternalError());

        Assert(std::abs(field_function.laplacian(p))<1e-10, ExcInternalError());
      }

      if (cell->is_locally_owned())
        {
          // more evil: use a corner so we can end up in a different cell
          Point<dim> p = cell->vertex(0);
          Assert(std::abs(field_function.value (p)-(p[0]+2.))<1e-10, ExcInternalError());

          Tensor<1,dim> gradient = field_function.gradient (p);
          Tensor<1,dim> exact_gradient;
          exact_gradient[0]=1.0;
          Assert((gradient-exact_gradient).norm()<1e-10, ExcInternalError());

          Assert(std::abs(field_function.laplacian(p))<1e-10, ExcInternalError());
        }
    }

  // also hit the code path for list version (no need to cover gradients/laplace here):
  std::vector<double> values(points.size());
  field_function.value_list(points, values, 0);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
