// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


/**
 * Test MGTwoLevelTransfer::prolongate_and_add() and
 * MGTwoLevelTransfer::interpolate() between uniformly refined meshes for
 * systems of simplex elements (FE_SimplexP, FE_SimplexDGP) with one and with
 * dim components. A linear vector
 * field with different values in each component is contained in both the
 * coarse and the fine space, so both operations have to be exact.
 *
 * Like mg_transfer_a_05 and mg_transfer_a_06 but for FESystem.
 */

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <deal.II/numerics/vector_tools.h>

#include "mg_transfer_util.h"

using namespace dealii;

template <int dim>
class LinearFunction : public Function<dim>
{
public:
  LinearFunction(const unsigned int n_components)
    : Function<dim>(n_components)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    double result = 1. + component;
    for (unsigned int d = 0; d < dim; ++d)
      result += (1. + d + dim * component) * p[d];
    return result;
  }
};

template <int dim, typename Number>
void
do_test(const FiniteElement<dim> &fe, const unsigned int n_components)
{
  // create coarse and fine grid
  Triangulation<dim> tria_coarse;
  Triangulation<dim> tria_fine;
  GridGenerator::subdivided_hyper_cube_with_simplices(tria_coarse, 1);
  GridGenerator::subdivided_hyper_cube_with_simplices(tria_fine, 1);
  tria_coarse.refine_global(1);
  tria_fine.refine_global(2);

  const FESystem<dim> fe_system(fe, n_components);

  const MappingFE<dim> mapping(FE_SimplexP<dim>(1));

  // setup dof-handlers
  DoFHandler<dim> dof_handler_fine(tria_fine);
  dof_handler_fine.distribute_dofs(fe_system);

  DoFHandler<dim> dof_handler_coarse(tria_coarse);
  dof_handler_coarse.distribute_dofs(fe_system);

  // setup constraint matrix
  AffineConstraints<Number> constraint_coarse;
  constraint_coarse.close();

  AffineConstraints<Number> constraint_fine;
  constraint_fine.close();

  // setup transfer operator
  MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>> transfer;
  transfer.reinit(dof_handler_fine,
                  dof_handler_coarse,
                  constraint_fine,
                  constraint_coarse);

  const LinearFunction<dim> function(n_components);

  LinearAlgebra::distributed::Vector<Number> vec_fine, vec_coarse, vec_fine_ref,
    vec_coarse_ref;
  initialize_dof_vector(vec_fine,
                        dof_handler_fine,
                        numbers::invalid_unsigned_int);
  initialize_dof_vector(vec_fine_ref,
                        dof_handler_fine,
                        numbers::invalid_unsigned_int);
  initialize_dof_vector(vec_coarse,
                        dof_handler_coarse,
                        numbers::invalid_unsigned_int);
  initialize_dof_vector(vec_coarse_ref,
                        dof_handler_coarse,
                        numbers::invalid_unsigned_int);

  VectorTools::interpolate(mapping, dof_handler_fine, function, vec_fine_ref);
  VectorTools::interpolate(mapping,
                           dof_handler_coarse,
                           function,
                           vec_coarse_ref);

  const auto print_error =
    [](const std::string                                &label,
       LinearAlgebra::distributed::Vector<Number>       &vec,
       const LinearAlgebra::distributed::Vector<Number> &vec_ref) {
      vec -= vec_ref;
      const Number error = vec.linfty_norm() / vec_ref.linfty_norm();
      deallog << label << " error: " << (error < 1e-12 ? 0. : error)
              << std::endl;
    };

  // test prolongation: coarse interpolant -> fine interpolant
  vec_fine = 0.;
  transfer.prolongate_and_add(vec_fine, vec_coarse_ref);
  print_error("prolongate ", vec_fine, vec_fine_ref);

  // test interpolation: fine interpolant -> coarse interpolant
  vec_coarse = 0.;
  transfer.interpolate(vec_coarse, vec_fine_ref);
  print_error("interpolate", vec_coarse, vec_coarse_ref);
}

template <int dim, typename Number>
void
test(const FiniteElement<dim> &fe)
{
  for (const unsigned int n_components : {1u, static_cast<unsigned int>(dim)})
    {
      deallog.push(fe.get_name() + "^" + std::to_string(n_components));
      do_test<dim, Number>(fe, n_components);
      deallog.pop();
    }
}

int
main()
{
  initlog();

  deallog.precision(8);

  for (unsigned int degree = 1; degree <= 3; ++degree)
    {
      test<2, double>(FE_SimplexP<2>(degree));
      test<2, double>(FE_SimplexDGP<2>(degree));
    }

  for (unsigned int degree = 1; degree <= 3; ++degree)
    {
      test<3, double>(FE_SimplexP<3>(degree));
      test<3, double>(FE_SimplexDGP<3>(degree));
    }
}
