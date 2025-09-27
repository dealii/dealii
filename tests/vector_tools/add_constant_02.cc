// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test VectorTools::add_constant() in parallel

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



// 0 (dim components)
// 1
// 2
// 3
template <int dim>
class Reference : public Function<dim>
{
public:
  Reference()
    : Function<dim>(dim + 3)
  {}

  double
  value(const Point<dim> &p, const unsigned int c) const
  {
    if (c == dim)
      return 1.0;
    if (c == dim + 1)
      return 2.0;
    if (c == dim + 2)
      return 3.0;
    return 0.0;
  }
};



template <int dim, class VectorType>
void
test()
{
  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FESystem<dim>   fe(FE_RaviartThomas<dim>(1),
                   1,
                   FE_Q<dim>(2),
                   1,
                   FE_DGQ<dim>(2),
                   1,
                   FE_DGP<dim>(1),
                   1);
  DoFHandler<dim> dofhandler(tria);
  dofhandler.distribute_dofs(fe);

  VectorType vec(dofhandler.locally_owned_dofs(), MPI_COMM_WORLD);
  vec = 0.0;

  for (unsigned int c = dim; c < dim + 3; ++c)
    VectorTools::add_constant(vec, dofhandler, c, 1.0 * (c - dim + 1));

  Vector<double> cellwise_errors(tria.n_active_cells());
  QIterated<dim> quadrature(QTrapezoid<1>(), 5);

  const dealii::Function<dim, double> *w = nullptr;
  const IndexSet                       relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dofhandler);
  VectorType ghosted(dofhandler.locally_owned_dofs(),
                     relevant_dofs,
                     MPI_COMM_WORLD);
  ghosted = vec;
  VectorTools::integrate_difference(dofhandler,
                                    ghosted,
                                    Reference<dim>(),
                                    cellwise_errors,
                                    quadrature,
                                    VectorTools::L2_norm);

  const double error = VectorTools::compute_global_error(tria,
                                                         cellwise_errors,
                                                         VectorTools::L2_norm);

  AssertThrow(error < 1e-10, ExcMessage("Error in integrate_difference"));
}

template <int dim, class VectorType>
void
test_simplex()
{
  parallel::shared::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube_with_simplices(tria, 4);

  FESystem<dim>   fe(FE_SimplexP<dim>(1),
                   dim,
                   FE_SimplexP<dim>(1),
                   1,
                   FE_SimplexP<dim>(1),
                   1,
                   FE_SimplexDGP<dim>(1),
                   1);
  DoFHandler<dim> dofhandler(tria);
  dofhandler.distribute_dofs(fe);

  VectorType vec(dofhandler.locally_owned_dofs(), MPI_COMM_WORLD);
  vec = 0.0;
  for (unsigned int c = dim; c < dim + 3; ++c)
    VectorTools::add_constant(vec, dofhandler, c, 1.0 * (c - dim + 1));

  Vector<double>     cellwise_errors(tria.n_active_cells());
  QGaussSimplex<dim> quadrature(2);

  const dealii::Function<dim, double> *w = nullptr;
  MappingFE<dim>                       mapping(FE_SimplexP<dim>(1));

  const IndexSet relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dofhandler);
  VectorType ghosted(dofhandler.locally_owned_dofs(),
                     relevant_dofs,
                     MPI_COMM_WORLD);
  ghosted = vec;
  VectorTools::integrate_difference(mapping,
                                    dofhandler,
                                    ghosted,
                                    Reference<dim>(),
                                    cellwise_errors,
                                    quadrature,
                                    VectorTools::L2_norm);

  const double error = VectorTools::compute_global_error(tria,
                                                         cellwise_errors,
                                                         VectorTools::L2_norm);

  AssertThrow(error < 1e-10, ExcMessage("Error in integrate_difference"));
}


template <int dim>
void
big()
{
  test<dim, LinearAlgebra::distributed::Vector<double>>();
  test_simplex<dim, LinearAlgebra::distributed::Vector<double>>();

#ifdef DEAL_II_WITH_TRILINOS
  test<dim, TrilinosWrappers::MPI::Vector>();
  test_simplex<dim, TrilinosWrappers::MPI::Vector>();
#endif
#ifdef DEAL_II_WITH_PETSC
  test<dim, PETScWrappers::MPI::Vector>();
  test_simplex<dim, PETScWrappers::MPI::Vector>();
#endif
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  big<2>();
  big<3>();

  deallog << "OK" << std::endl;
}
