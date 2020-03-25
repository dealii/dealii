// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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


#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

// Test that, for a single thread, CellSimilarity calculations return 'none'
// when we are using curved cells. This is a regression test for a bug where,
// for this particular manifold, some cells would be marked as translations
// even though they were not (they varied in curvature of faces). This (small)
// error causes the Laplace solver to converge at order 2 instead of order 5,
// resulting in much higher L2 errors when the grid is refined.


static const unsigned int               fe_order          = 4;
static const dealii::types::boundary_id boundary_id       = 0;
static const dealii::types::manifold_id cubic_manifold_id = 1;
static const double                     pi                = numbers::PI;

// ----------------------------------------------------------------------------
// Manufactured solution and manufactured forcing
// ----------------------------------------------------------------------------

template <int dim>
class HardManufacturedSolution : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &point, const unsigned int) const
  {
    const double &x = point[0];
    const double &y = point[1];

    return (Utilities::fixed_power<3>(y) +
            std::exp(-Utilities::fixed_power<2>(y)) +
            std::sin(4.5 * Utilities::fixed_power<2>(y)) + std::sin(20 * y)) *
           (20 * std::cos(4 * pi * x) + 0.1 * std::sin(20 * pi * x) -
            80 * std::sin(6 * pi * x));
  }
};



template <int dim>
class HardManufacturedForcing : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &point, const unsigned int) const
  {
    const double &x = point[0];
    const double &y = point[1];
    return -40.0 *
             (Utilities::fixed_power<3>(y) +
              std::exp(-Utilities::fixed_power<2>(y)) +
              std::sin(4.5 * Utilities::fixed_power<2>(y)) + std::sin(20 * y)) *
             (-8.0 * pi * pi * std::cos(4 * pi * x) -
              1.0 * pi * pi * std::sin(20 * pi * x) +
              72.0 * pi * pi * std::sin(6 * pi * x)) -
           (4 * Utilities::fixed_power<2>(y) *
              std::exp(-Utilities::fixed_power<2>(y)) -
            81.0 * Utilities::fixed_power<2>(y) *
              std::sin(4.5 * Utilities::fixed_power<2>(y)) +
            6 * y + 9.0 * std::cos(4.5 * Utilities::fixed_power<2>(y)) -
            2 * std::exp(-Utilities::fixed_power<2>(y)) -
            400 * std::sin(20 * y)) *
             (20 * std::cos(4 * pi * x) + 0.1 * std::sin(20 * pi * x) -
              80 * std::sin(6 * pi * x));
  }
};


// ----------------------------------------------------------------------------
// Description of the curved geometry
// ----------------------------------------------------------------------------

template <int dim>
class PushForward : public Function<dim>
{
public:
  PushForward()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    switch (component)
      {
        case 0:
          return point[0];
        case 1:
          {
            const double &x = point[0];
            return point[1] + 0.25 * (2 * x - 1.0) * (x - 1.0) * x;
          }
        default:
          Assert(false, ExcNotImplemented());
      }
    return std::numeric_limits<double>::quiet_NaN();
  }
};

template <int dim>
class PullBack : public Function<dim>
{
public:
  PullBack()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    switch (component)
      {
        case 0:
          return point[0];
        case 1:
          {
            const double &x = point[0];
            return point[1] - 0.25 * (2 * x - 1.0) * (x - 1.0) * x;
          }
        default:
          Assert(false, ExcNotImplemented());
      }
    return std::numeric_limits<double>::quiet_NaN();
  }
};

/**
 * A struct containing the push forward and pull back functions. This is
 * only used with CubicRoofManifold: the only reason this struct exists is
 * so that CubicRoofManifold holds instances of these functions (which are
 * fed into FunctionManifold, so they must exist at that point in the
 * constructor call) and destroys them after destroying its FunctionManifold
 * part.
 */
template <int dim>
struct CubicRoofFunctions
{
  PushForward<dim> forward;
  PullBack<dim>    backward;
};

/**
 * A manifold describing the cubic roof.
 */
template <int dim>
class CubicRoofManifold : private CubicRoofFunctions<dim>,
                          public FunctionManifold<dim>
{
public:
  CubicRoofManifold()
    : FunctionManifold<dim>(this->forward, this->backward)
  {}
};

template <int dim>
std::shared_ptr<Manifold<dim>>
cubic_roof(Triangulation<dim> &triangulation)
{
  std::shared_ptr<Manifold<dim>> boundary(new CubicRoofManifold<dim>());
  GridGenerator::hyper_cube(triangulation);

  triangulation.set_all_manifold_ids(cubic_manifold_id);
  triangulation.set_manifold(cubic_manifold_id, *boundary);
  return boundary;
}

// ----------------------------------------------------------------------------
// The Laplace solver that has trouble with 1 thread, but works well with 2
// ----------------------------------------------------------------------------
template <int dim>
class JxWError
{
public:
  JxWError(const unsigned int n_global_refines);

  double
  run();

protected:
  std::shared_ptr<Function<dim>> manufactured_solution;
  std::shared_ptr<Function<dim>> manufactured_forcing;

  std::shared_ptr<Manifold<dim>> boundary_manifold;
  Triangulation<dim>             triangulation;
  FE_Q<dim>                      finite_element;
  DoFHandler<dim>                dof_handler;
  QGauss<dim>                    cell_quadrature;
  MappingQGeneric<dim>           cell_mapping;

  AffineConstraints<double> all_constraints;
  SparsityPattern           sparsity_pattern;
  SparseMatrix<double>      system_matrix;
  Vector<double>            system_rhs;
  Vector<double>            solution;

  void
  setup_dofs();
  void
  setup_matrices();
  double
  solve();
};



template <int dim>
JxWError<dim>::JxWError(const unsigned int n_global_refines)
  : manufactured_solution(new HardManufacturedSolution<dim>())
  , manufactured_forcing(new HardManufacturedForcing<dim>())
  , finite_element(fe_order)
  , dof_handler(triangulation)
  , cell_quadrature(fe_order + 1)
  , cell_mapping(fe_order)

{
  boundary_manifold = cubic_roof(triangulation);
  triangulation.refine_global(n_global_refines);
}



template <int dim>
void
JxWError<dim>::setup_dofs()
{
  dof_handler.distribute_dofs(finite_element);
  VectorTools::interpolate_boundary_values(cell_mapping,
                                           dof_handler,
                                           boundary_id,
                                           *manufactured_solution,
                                           all_constraints);
  all_constraints.close();

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dynamic_sparsity_pattern,
                                  all_constraints,
                                  /*keep_constrained_dofs=*/false);
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);
}


template <int dim>
void
JxWError<dim>::setup_matrices()
{
  system_matrix.reinit(sparsity_pattern);
  system_rhs.reinit(dof_handler.n_dofs());

  const UpdateFlags flags = update_values | update_gradients |
                            update_JxW_values | update_quadrature_points;
  FEValues<dim> fe_values(cell_mapping, finite_element, cell_quadrature, flags);

  const unsigned int dofs_per_cell = finite_element.dofs_per_cell;
  FullMatrix<double> cell_system(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      cell->get_dof_indices(local_dof_indices);
      cell_system = 0.0;
      cell_rhs    = 0.0;
      fe_values.reinit(cell);

      for (unsigned int q_point_n = 0;
           q_point_n < fe_values.n_quadrature_points;
           ++q_point_n)
        {
          const double point_forcing =
            manufactured_forcing->value(fe_values.quadrature_point(q_point_n));

          for (unsigned int test_n = 0; test_n < dofs_per_cell; ++test_n)
            {
              for (unsigned int trial_n = 0; trial_n < dofs_per_cell; ++trial_n)
                {
                  cell_system(test_n, trial_n) +=
                    fe_values.JxW(q_point_n) *
                    (fe_values.shape_grad(test_n, q_point_n) *
                     fe_values.shape_grad(trial_n, q_point_n));
                }

              cell_rhs[test_n] += fe_values.JxW(q_point_n) *
                                  fe_values.shape_value(test_n, q_point_n) *
                                  point_forcing;
            }
        }

      all_constraints.distribute_local_to_global(
        cell_system, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
}



template <int dim>
double
JxWError<dim>::solve()
{
  {
    SolverControl      solver_control(std::max(types::global_dof_index(100),
                                          static_cast<types::global_dof_index>(
                                            system_rhs.size())),
                                 1e-14 * system_rhs.l2_norm(),
                                 false,
                                 false);
    SolverCG<>         solver(solver_control);
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solution.reinit(system_rhs);
    solver.solve(system_matrix, solution, system_rhs, preconditioner);
    all_constraints.distribute(solution);
  }

  Vector<double> cell_l2_error(triangulation.n_cells());
  VectorTools::integrate_difference(
    cell_mapping,
    dof_handler,
    solution,
    *manufactured_solution,
    cell_l2_error,
    // use QIterated to avoid spurious superconvergence
    QIterated<dim>(QGauss<1>(finite_element.degree), 2),
    VectorTools::L2_norm);

  const double l2_error =
    VectorTools::compute_global_error(triangulation,
                                      cell_l2_error,
                                      VectorTools::L2_norm);
  return l2_error;
}



template <int dim>
double
JxWError<dim>::run()
{
  setup_dofs();
  setup_matrices();
  const double error = solve();
  return error;
}



int
main(int argc, char **argv)
{
  // Use exactly one thread so that the CellSimilarity checks are not disabled.
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  static const int dim = 2;

  initlog();
  deallog << std::setprecision(10);
  for (unsigned int n_global_refines = 3; n_global_refines < 6;
       ++n_global_refines)
    {
      JxWError<dim> solver(n_global_refines);
      deallog << "L2 error: " << solver.run() << std::endl;
    }
}
