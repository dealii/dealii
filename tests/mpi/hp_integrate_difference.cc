// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test VectorTools::integrate_difference for parallel computations
// with the hp::DoFHandler. This includes applying hanging node
// constraints and consequently verifies that we compute them
// correctly.
//
// The way this test works is this: We create a domain [-1,1]^d
// and use the following function f(x,y) (or f(x,y,z)):
// - in the lower left quadrant,  f(x,y)=xy
// - in the lower right quadrant, f(x,y)=xy+(xy)^2
// - in the upper left quadrant,  f(x,y)=xy+(xy)^3
// - in the upper right quadrant, f(x,y)=xy+(xy)^2+(xy)^3+(xy)*4
//
// We interpolate this function onto a finite element space that is
// chosen as follows:
// - in the lower left quadrant,  Q1
// - in the lower right quadrant, Q2
// - in the upper left quadrant,  Q3
// - in the upper right quadrant, Q4
// In other words, the function f(...) is in the space.
//
// We can then run two tests with it:
// - Interpolate the function onto the finite element space and compute
//   its L2 norm. This can be done analytically. In particular, the
//   area under the square of the functions above is, for the four
//   quadrants:
//   . 1/9
//   . 47/1800
//   . 2332/11025
//   . 1816349/3175200
//   This makes the total sum under the square of the function equal to
//   2923673/3175200 and the L2 norm under the function equal to
//   sqrt(5847346)/2520, which is about 0.9595748472.
//   (In 3d, we also integrate the same function over z=-1..1, so the
//   volume under f(...)^2 increases by a factor of 2, and the L2 norm
//   by a factor of sqrt(2); the numerical value is then 1.357043763.)
// - Interpolate the function onto the finite element space and compute
//   the L2 norm of the difference between the interpolated function
//   and the original function. This should be zero.
// This test does both.

#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
class CheckFunction : public Function<dim>
{
public:
  double
  value(const Point<dim> &p, const unsigned int) const
  {
    const double x = p[0];
    const double y = p[1];

    // function is bilinear everywhere
    double f = x * y;

    // on the right half of the domain, add something biquadratic that's
    // zero at x=0
    if (x >= 0)
      f += x * x * y * y;


    // in the top half of the domain, add something bicubic that's
    // zero at y=0
    if (y >= 0)
      f += y * y * y * x * x * x;

    // in the top right quadrant, add something biquartic that's
    // zero at x=0 and y=0
    if (x >= 0 && y >= 0)
      f += x * x * x * x * y * y * y * y;

    return f;
  }
};


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(3);

  const hp::FECollection<dim> fe(FE_Q<dim>(1),
                                 FE_Q<dim>(2),
                                 FE_Q<dim>(3),
                                 FE_Q<dim>(4));
  DoFHandler<dim>             dof_handler(tr);

  // set DoF indices as described at the top of the file
  for (auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const double x = cell->center()[0];
        const double y = cell->center()[1];

        if (x < 0 && y < 0)
          cell->set_active_fe_index(0);
        else if (x > 0 && y < 0)
          cell->set_active_fe_index(1);
        else if (x < 0 && y > 0)
          cell->set_active_fe_index(2);
        else if (x > 0 && y > 0)
          cell->set_active_fe_index(3);
      }

  dof_handler.distribute_dofs(fe);

  // interpolate the function above onto the finite element space
  TrilinosWrappers::MPI::Vector interpolated(dof_handler.locally_owned_dofs(),
                                             MPI_COMM_WORLD);
  VectorTools::interpolate(dof_handler, CheckFunction<dim>(), interpolated);

  // then also apply constraints
  const IndexSet relevant_set =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  AffineConstraints<double> hanging_node_constraints(
    dof_handler.locally_owned_dofs(), relevant_set);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          hanging_node_constraints);
  hanging_node_constraints.close();
  hanging_node_constraints.distribute(interpolated);

  // extract a vector that has ghost elements
  TrilinosWrappers::MPI::Vector x_rel(relevant_set, MPI_COMM_WORLD);
  x_rel = interpolated;

  // Create a sufficiently high order quadrature formula
  hp::QCollection<dim> quadrature(QGauss<dim>(3),
                                  QGauss<dim>(4),
                                  QGauss<dim>(5),
                                  QGauss<dim>(6));

  {
    // integrate the difference between the function above and
    // the zero function. for this case, we can compute the exact values
    // by hand. the ones printed in the output are correct
    Vector<float> results(tr.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      x_rel,
                                      Functions::ZeroFunction<dim>(),
                                      results,
                                      quadrature,
                                      VectorTools::L2_norm);
    const double global =
      VectorTools::compute_global_error(tr, results, VectorTools::L2_norm);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      deallog << "L2 norm = " << global << std::endl;

    Assert(std::fabs(global - std::sqrt(5847346.) / 2520. *
                                (dim == 3 ? std::sqrt(2) : 1)) < 1e-7,
           ExcInternalError());
  }


  {
    // Now also integrate the difference between the function above and
    // the its interpolant. This should then of course be zero
    Vector<float> results(tr.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      x_rel,
                                      CheckFunction<dim>(),
                                      results,
                                      quadrature,
                                      VectorTools::L2_norm);
    const double global =
      VectorTools::compute_global_error(tr, results, VectorTools::L2_norm);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      deallog << "L2 error = " << global << std::endl;

    Assert(std::fabs(global) < 1e-15, ExcInternalError());
  }
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      deallog.push("2d");
      test<2>();
      deallog.pop();

      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
}
