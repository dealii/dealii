// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// check intergrid transfer

#include "../tests.h"
#include <deal.II/base/function_derivative.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out_stack.h>
#include <fstream>
#include <iomanip>

#include <deal.II/grid/intergrid_map.h>
#include <deal.II/dofs/dof_handler.h>

namespace DoFToolsEx
{
  /// transfers solution between differently refined grids
  template <int dim, class InVector, class OutVector>
  void transfer(const DoFHandler<dim> &source_dof, const InVector &source_vector, const DoFHandler<dim> &target_dof, OutVector &target_vector);
}

/**
 * Detailed desc.
 * Does it make sense to have InVector and OutVector templates? Shouldn't just one Vector be enough?
 *
 * @param
 *
 */
template <int dim, class InVector, class OutVector>
void
DoFToolsEx::transfer(const DoFHandler<dim> &source_dof, const InVector &source_vector, const DoFHandler<dim> &target_dof, OutVector &target_vector)
{
  // any sanity tests? Trias derived from same coarse grid?

  // build mappings between the cells
  InterGridMap<DoFHandler<dim> > source2target, target2source;
  source2target.make_mapping(source_dof, target_dof);
  target2source.make_mapping(target_dof, source_dof);

  // setup temporary vector
  InVector local_dofs(source_dof.get_fe().dofs_per_cell);

  // iterate over all active source cells
  typedef typename DoFHandler<dim>::active_cell_iterator cell_iterator;
  cell_iterator cell = source_dof.begin_active(), endc = source_dof.end();
  for (; cell != endc; ++cell)
    {
      if (cell->level() == source2target[cell]->level())
        {
          // cell has not been coarsend, but possibly been refined
          cell->get_dof_values(source_vector, local_dofs);
          source2target[cell]->set_dof_values_by_interpolation(local_dofs, target_vector);
        }
      else
        {
          // the source cell has been coarsened
          Assert(cell->level() > source2target[cell]->level(), ExcInternalError());
          target2source[ source2target[cell] ]->get_interpolated_dof_values(source_vector, local_dofs);
          source2target[cell]->set_dof_values(local_dofs, target_vector);
        }
    }

  // handle hanging node constraints inside or outside of this function?
}


class TestFunction : public Function<2>
{
public:
  TestFunction() : Function<2>() {}
  virtual double value(const Point<2> &p,
                       const unsigned int) const
  {
    return std::sin(3.14159*p(0))*std::sin(3.14159*p(1));
  }
};

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);


  // build test-case trias
  Triangulation<2> tria;
  Triangulation<2> refined_tria;
  Triangulation<2> coarse_tria;
  Triangulation<2> both_tria;

  // make grids
  GridGenerator::hyper_cube(tria, 0, 1);
  tria.refine_global(2);
  GridGenerator::hyper_cube(refined_tria, 0, 1);
  refined_tria.refine_global(3);
  GridGenerator::hyper_cube(coarse_tria, 0, 1);
  coarse_tria.refine_global(1);
  GridGenerator::hyper_cube(both_tria, 0, 1);
  both_tria.refine_global(2);
  Triangulation<2>::active_cell_iterator it = both_tria.begin_active();
  it->set_refine_flag();
  (++it)->set_refine_flag();
  for (int i = 1; i < 8; ++i, ++it)
    ;
  it->set_coarsen_flag();
  (++it)->set_coarsen_flag();
  (++it)->set_coarsen_flag();
  (++it)->set_coarsen_flag();
  both_tria.execute_coarsening_and_refinement();

  // finite element
  FE_Q<2> fe(1);

  // dof handlers
  DoFHandler<2> dof(tria), refined_dof(refined_tria), coarse_dof(coarse_tria), both_dof(both_tria);
  dof.distribute_dofs(fe);
  refined_dof.distribute_dofs(fe);
  coarse_dof.distribute_dofs(fe);
  both_dof.distribute_dofs(fe);

  // interpolate test function onto tria
  Vector<double> sol(dof.n_dofs());
  TestFunction test_function;
  VectorTools::interpolate(dof, test_function, sol);

  // setup solution vectors
  Vector<double> refined_sol(refined_dof.n_dofs()), coarse_sol(coarse_dof.n_dofs()), both_sol(both_dof.n_dofs());

  // transfer solutions
  DoFToolsEx::transfer(dof, sol, refined_dof, refined_sol);
  DoFToolsEx::transfer(dof, sol, coarse_dof, coarse_sol);
  DoFToolsEx::transfer(dof, sol, both_dof, both_sol);

  // handle hanging nodes
  ConstraintMatrix both_constraints;
  DoFTools::make_hanging_node_constraints(both_dof, both_constraints);
  both_constraints.close();
  both_constraints.distribute(both_sol);

  // reference output using DataOut to create seperate files
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof);
  data_out.add_data_vector(sol, "base");
  data_out.build_patches();
  data_out.write_gnuplot(logfile);

  data_out.clear();
  data_out.attach_dof_handler(refined_dof);
  data_out.add_data_vector(refined_sol, "refined");
  data_out.build_patches();
  data_out.write_gnuplot(logfile);

  data_out.clear();
  data_out.attach_dof_handler(coarse_dof);
  data_out.add_data_vector(coarse_sol, "coarse");
  data_out.build_patches();
  data_out.write_gnuplot(logfile);

  data_out.clear();
  data_out.attach_dof_handler(both_dof);
  data_out.add_data_vector(both_sol, "both");
  data_out.build_patches();
  data_out.write_gnuplot(logfile);

  // test output using DataOutStack
  DataOutStack<2> data_out_stack;
  data_out_stack.declare_data_vector("dof", DataOutStack<2>::dof_vector);
  data_out_stack.new_parameter_value(0.,0.);
  data_out_stack.attach_dof_handler(dof);
  data_out_stack.add_data_vector(sol, "dof");
  data_out_stack.build_patches();
  data_out_stack.finish_parameter_value();
  data_out_stack.write_gnuplot(logfile);
}
