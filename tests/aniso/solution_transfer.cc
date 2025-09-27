// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    : Function<dim>(){};

  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    double ret_value = sin(p[0] * 4) * cos(p[1] * 4);
    if (dim == 3)
      ret_value *= sin(5 * p[2] + 1);
    return ret_value;
  };
};


template <int dim>
void
transfer(std::ostream &out)
{
  MyFunction<dim>    function;
  Triangulation<dim> tria(Triangulation<dim>::allow_anisotropic_smoothing);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  FE_DGQ<dim>     fe(1);
  DoFHandler<dim> dof_handler(tria);
  Vector<double>  solution;
  MappingQ<dim>   mapping(1);
  DataOut<dim>    data_out;

  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());

  VectorTools::interpolate(mapping, dof_handler, function, solution);

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  deallog << "Initial solution" << std::endl << std::endl;
  data_out.write_gnuplot(out);

  SolutionTransfer<dim> soltrans(dof_handler);

  // test a): pure refinement
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  for (; cell != endc; ++cell)
    cell->set_refine_flag(RefinementCase<dim>::cut_x);

  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_coarsening_and_refinement(solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe);

  Vector<double> new_solution(dof_handler.n_dofs());
  soltrans.interpolate(new_solution);
  solution.reinit(dof_handler.n_dofs());
  solution = new_solution;

  data_out.clear_data_vectors();
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  deallog << "Interpolated/transferred solution after pure refinement"
          << std::endl
          << std::endl;
  data_out.write_gnuplot(out);

  // test b): with coarsening
  SolutionTransfer<dim> soltrans2(dof_handler);
  cell = tria.begin_active(tria.n_levels() - 1);
  endc = tria.end(tria.n_levels() - 1);
  for (; cell != endc; ++cell)
    cell->set_coarsen_flag();
  Vector<double> old_solution = solution;
  tria.prepare_coarsening_and_refinement();
  soltrans2.prepare_for_coarsening_and_refinement(old_solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe);
  solution.reinit(dof_handler.n_dofs());
  soltrans2.interpolate(solution);

  data_out.clear_data_vectors();
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  deallog << "Interpolated/transferred solution after coarsening" << std::endl
          << std::endl;
  data_out.write_gnuplot(out);
}


int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(4);

  transfer<2>(deallog.get_file_stream());
  transfer<3>(deallog.get_file_stream());
}
