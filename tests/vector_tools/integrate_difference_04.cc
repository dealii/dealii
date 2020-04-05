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


// Test integrate_difference and compute_global_error in parallel
// see integrate_difference_02.cc for the serial version

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



// x+y+z, x^2+y^2, z+xy
// div = 1+2y+1
template <int dim>
class Ref : public Function<dim>
{
public:
  Ref()
    : Function<dim>(dim)
  {}

  double
  value(const Point<dim> &p, const unsigned int c) const
  {
    if (c == 0)
      return p[0] + p[1] + ((dim == 3) ? p[2] : 0.0);
    if (c == 1)
      return p[0] * p[0] + p[1] * p[1];
    if (c == 2)
      return p[2] + p[0] * p[1];
    return 0.0;
  }
};



template <int dim>
void
test(VectorTools::NormType norm, double value, double exp = 2.0)
{
  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_metis);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);

  FESystem<dim>   fe(FE_Q<dim>(4), dim);
  DoFHandler<dim> dofh(tria);
  dofh.distribute_dofs(fe);

  TrilinosWrappers::MPI::Vector interpolated(dofh.locally_owned_dofs(),
                                             MPI_COMM_WORLD);

  VectorTools::interpolate(dofh, Ref<dim>(), interpolated);

  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs(dofh, relevant_set);
  TrilinosWrappers::MPI::Vector solution(relevant_set, MPI_COMM_WORLD);
  solution = interpolated;

  Vector<double> cellwise_errors(tria.n_active_cells());
  QIterated<dim> quadrature(QTrapez<1>(), 5);

  const dealii::Function<dim, double> *w = nullptr;
  VectorTools::integrate_difference(dofh,
                                    solution,
                                    Functions::ZeroFunction<dim>(dim),
                                    cellwise_errors,
                                    quadrature,
                                    norm,
                                    w,
                                    exp);

  const double error =
    VectorTools::compute_global_error(tria, cellwise_errors, norm, exp);

  const double difference = std::abs(error - value);
  deallog << "computed: " << error << " expected: " << value
          << " difference: " << difference << std::endl;
  Assert(difference < 2e-3, ExcMessage("Error in integrate_difference"));
}

template <int dim>
void
test()
{
  deallog << "Hdiv_seminorm:" << std::endl;
  // sqrt(\int (div f)^2 = sqrt(\int (1+2y+1)^2)
  test<dim>(VectorTools::Hdiv_seminorm, 2.0 * std::sqrt(7.0 / 3.0));

  deallog << "L2_norm:" << std::endl;
  // sqrt(\int_\Omega f^2) = sqrt(\int (x+y+z)^2+(x^2+y^2)^2+(z+xy)^2)
  test<dim>(VectorTools::L2_norm, std::sqrt(229.0 / 60.0));

  deallog << "H1_seminorm:" << std::endl;
  // sqrt( \int sum_k | d/dxi f_k |_0^2 )
  // = sqrt( \int 3+ (2x)^2+(2y)^2 + y^2+x^2+1
  // = sqrt( 22/3  )
  test<dim>(VectorTools::H1_seminorm, std::sqrt(22.0 / 3.0));

  deallog << "H1_norm:" << std::endl;
  test<dim>(VectorTools::H1_norm, std::sqrt(229.0 / 60.0 + 22.0 / 3.0));

  deallog << "L1_norm:" << std::endl;
  test<dim>(VectorTools::L1_norm, 35.0 / 12.0);

  deallog << "Linfty_norm:" << std::endl;
  test<dim>(VectorTools::Linfty_norm, 3.0);

  deallog << "mean:" << std::endl;
  // int -(x+y+z + x^2+y^2 + z+xy)
  test<dim>(VectorTools::mean, -35.0 / 12.0);

  deallog << "Lp_norm:" << std::endl;
  // (int (x+y+z)^p+(x^2+y^2)^p+(z+xy)^p) ) ^ 1./p
  // = std::pow(9937.0/1680.0, 1.0/3.0)
  test<dim>(VectorTools::Lp_norm, std::pow(9937.0 / 1680.0, 1.0 / 3.0), 3.0);

  deallog << "W1p_seminorm:" << std::endl;
  // ( \int_K sum_k | d/dxi f_k |_0^2^(p/2) )^1/p
  // ( integrate 3^(3/2) + ((2x)^2+(2y)^2)^(3/2) + (y^2+x^2+1)^(3/2) ) ^1/p
  // = (12.4164) ^1/3 = 2.31560
  test<dim>(VectorTools::W1p_seminorm, 2.31560, 3.0);

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  deallog << std::setprecision(10);
  test<3>();
}
