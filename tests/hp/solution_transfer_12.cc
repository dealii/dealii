// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
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

#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <vector>

#include "../tests.h"

// a linear function that should be transferred exactly with Q1 and Q2
// elements
// This is a copy-paste of solution_transfer.cc test,
// with the two exceptions:
// - FE_Q_Hierarchical
// - manual interpolation as there is no general interpolation implemented
//   for this FE currently, thus VectorTools::interpolate can't be used.
template <int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    : Function<dim>(){};

  virtual double
  value(const Point<dim> &p, const unsigned int) const
  {
    double f = 0.25 + 2 * p[0];
    if (dim > 1)
      f += 0.772 * p[1];
    if (dim > 2)
      f -= 3.112 * p[2];
    return f;
  };
};


template <int dim>
void
transfer(std::ostream &out)
{
  MyFunction<dim>    func;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(5 - dim);
  const unsigned int    max_degree = 6 - dim;
  hp::FECollection<dim> fe_q;
  for (unsigned int deg = 1; deg <= max_degree; ++deg)
    {
      fe_q.push_back(FE_Q_Hierarchical<dim>(deg));
    }
  DoFHandler<dim> q_dof_handler(tria);
  Vector<double>  q_solution;
  MappingQ<dim>   mapping(1);

  // refine a few cells
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  ++cell;
  ++cell;
  for (; cell != endc; ++cell)
    cell->set_refine_flag();
  tria.prepare_coarsening_and_refinement();
  tria.execute_coarsening_and_refinement();

  // randomly assign FE orders
  unsigned int counter = 0;
  {
    typename DoFHandler<dim>::active_cell_iterator cell = q_dof_handler
                                                            .begin_active(),
                                                   endc = q_dof_handler.end();

    for (; cell != endc; ++cell, ++counter)
      {
        if (counter < 15)
          cell->set_active_fe_index(1);
        else
          cell->set_active_fe_index(Testing::rand() % max_degree);
      }
  }

  q_dof_handler.distribute_dofs(fe_q);
  q_solution.reinit(q_dof_handler.n_dofs());

  AffineConstraints<double> cm;
  cm.close();

  // interpolation is not straight forward for FE_Q_Hierarchical.
  // Since we know that we do linear function, just get those DoFs
  // which correspond to vertices and assign a value of shape function there
  // while keeping other modes zero.
  {
    q_solution = 0.;

    hp::QCollection<dim> quad;
    quad.push_back(QTrapezoid<dim>());
    hp::FEValues<dim>                              hp_fe_val(fe_q,
                                quad,
                                update_values | update_quadrature_points);
    std::vector<types::global_dof_index>           local_dof_indices;
    typename DoFHandler<dim>::active_cell_iterator cell = q_dof_handler
                                                            .begin_active(),
                                                   endc = q_dof_handler.end();

    for (; cell != endc; ++cell)
      {
        hp_fe_val.reinit(cell, 0);
        const FEValues<dim> &fe_val = hp_fe_val.get_present_fe_values();
        local_dof_indices.resize(fe_val.dofs_per_cell);

        cell->get_dof_indices(local_dof_indices);

        Assert(local_dof_indices.size() >= fe_val.n_quadrature_points,
               ExcInternalError());

        for (unsigned int q = 0; q < fe_val.n_quadrature_points; ++q)
          {
            q_solution[local_dof_indices[q]] =
              func.value(fe_val.quadrature_point(q), 0);
          }
      }
  }

  SolutionTransfer<dim, Vector<double>> q_soltrans(q_dof_handler);



  // test b): do some coarsening and
  // refinement
  q_soltrans.clear();

  counter = 0;
  cell    = tria.begin_active();
  endc    = tria.end();
  for (; cell != endc; ++cell, ++counter)
    {
      if (counter > 120)
        cell->set_coarsen_flag();
      else if (Testing::rand() % 3 == 0)
        cell->set_refine_flag();
      else if (Testing::rand() % 3 == 3)
        cell->set_coarsen_flag();
    }

  counter = 0;
  {
    typename DoFHandler<dim>::active_cell_iterator cell = q_dof_handler
                                                            .begin_active(),
                                                   endc = q_dof_handler.end();

    for (; cell != endc; ++cell, ++counter)
      {
        if (counter > 20 && counter < 90)
          cell->set_future_fe_index(0);
        else
          cell->set_future_fe_index(Testing::rand() % max_degree);
      }
  }

  Vector<double> q_old_solution = q_solution;
  tria.prepare_coarsening_and_refinement();
  q_soltrans.prepare_for_coarsening_and_refinement(q_old_solution);
  tria.execute_coarsening_and_refinement();

  q_dof_handler.distribute_dofs(fe_q);
  q_solution.reinit(q_dof_handler.n_dofs());
  q_soltrans.interpolate(q_solution);

  // check correctness by comparing the values
  // on points of QGauss of order 2.
  {
    double                                         error = 0;
    const hp::QCollection<dim>                     quad(QGauss<dim>(2));
    hp::FEValues<dim>                              hp_fe_val(fe_q,
                                quad,
                                update_values | update_quadrature_points);
    std::vector<double>                            vals(quad[0].size());
    typename DoFHandler<dim>::active_cell_iterator cell = q_dof_handler
                                                            .begin_active(),
                                                   endc = q_dof_handler.end();

    for (; cell != endc; ++cell)
      {
        hp_fe_val.reinit(cell, 0);
        const FEValues<dim> &fe_val = hp_fe_val.get_present_fe_values();
        fe_val.get_function_values(q_solution, vals);
        for (unsigned int q = 0; q < fe_val.n_quadrature_points; ++q)
          {
            error +=
              std::fabs(func.value(fe_val.quadrature_point(q), 0) - vals[q]);
          }
      }
    deallog << "Error in interpolating hp-FE_Q_Hierarchical: " << error
            << std::endl;
  }
}


int
main()
{
  initlog();

  deallog << "   1D solution transfer" << std::endl;
  transfer<1>(deallog.get_file_stream());

  deallog << "   2D solution transfer" << std::endl;
  transfer<2>(deallog.get_file_stream());

  deallog << "   3D solution transfer" << std::endl;
  transfer<3>(deallog.get_file_stream());
}
