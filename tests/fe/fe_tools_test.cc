// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

/*
 * Author: Ralf Hartmann, University of Heidelberg, 1999
 */


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

// This is needed for C++ output:
#include <fstream>
#include <cmath>

const double pi=std::acos(-1.);

class TestFunction: public Function<2>
{
public:
  TestFunction() {}
  virtual ~TestFunction() {}

  virtual double value(const Point<2> &p,
                       const unsigned int component) const;
};


double
TestFunction::value(const Point<2> &p,
                    const unsigned int component) const
{
  Assert(component==0, ExcInternalError());
  return std::sin(pi*p(0))*std::cos(pi*p(1));
}




void make_grid (Triangulation<2> &triangulation)
{
  GridGenerator::hyper_cube (triangulation);

  const Point<2> p0=triangulation.begin_active()->vertex(0);

  triangulation.refine_global(1);

  for (unsigned int step=0; step<1; ++step)
    {
      Triangulation<2>::active_cell_iterator cell, endc;
      cell = triangulation.begin_active();
      endc = triangulation.end();

      for (; cell!=endc; ++cell)
        if (cell->vertex(0)==p0)
          {
            cell->set_refine_flag ();
            break;
          }

      triangulation.execute_coarsening_and_refinement ();
    };
}


template<int dim>
void test(const Triangulation<dim> &tria,
          const Mapping<dim> &mapping,
          const FiniteElement<dim> &fe1,
          const std::string &fe_string1,
          const FiniteElement<dim> &fe2,
          const std::string &fe_string2,
          const unsigned int testcase)
{
  DoFHandler<dim> dof_handler1 (tria);
  DoFHandler<dim> dof_handler2 (tria);

  dof_handler1.distribute_dofs (fe1);
  dof_handler2.distribute_dofs (fe2);

  Vector<double> function1(dof_handler1.n_dofs());
  Vector<double> function1_back(dof_handler1.n_dofs());
  Vector<double> function2(dof_handler2.n_dofs());

  ConstraintMatrix constraints1;
  DoFTools::make_hanging_node_constraints (dof_handler1, constraints1);
  constraints1.close();

  ConstraintMatrix constraints2;
  DoFTools::make_hanging_node_constraints (dof_handler2, constraints2);
  constraints2.close();

  QGauss<dim> quadrature(4);
  TestFunction function;
  VectorTools::project(mapping,
                       dof_handler1,
                       constraints1,
                       quadrature,
                       function,
                       function1);

  switch (testcase)
    {
    case 1:
      FETools::interpolate(dof_handler1, function1,
                           dof_handler2, constraints2, function2);
      break;
    case 2:
      FETools::back_interpolate(dof_handler1, constraints1, function1,
                                dof_handler2, constraints2, function1_back);
      break;
    case 3:
      FETools::interpolation_difference(dof_handler1, constraints1, function1,
                                        dof_handler2, constraints2, function1_back);
      break;
    default:
      Assert(false, ExcNotImplemented());
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler1);
  data_out.add_data_vector (function1, fe_string1);
  data_out.build_patches (2);
//    std::string file1_name=fe_string1+"_function.gnuplot";
//    std::ofstream file1(file1_name.c_str());
//    data_out.write_gnuplot(file1);
//    file1.close();
  data_out.clear();

  std::string file2_name=fe_string1+"_";
  switch (testcase)
    {
    case 1:
      data_out.attach_dof_handler (dof_handler2);
      data_out.add_data_vector (function2, fe_string2);
      file2_name+=fe_string2+"_interpolation.gnuplot";
      break;
    case 2:
      data_out.attach_dof_handler (dof_handler1);
      data_out.add_data_vector (function1_back, fe_string1);
      file2_name+=fe_string2+"_back_interpolation.gnuplot";
      break;
    case 3:
      data_out.attach_dof_handler (dof_handler1);
      data_out.add_data_vector (function1_back, fe_string1);
      file2_name+=fe_string2+"_interpolation_diff.gnuplot";
      break;
    default:
      Assert(false, ExcNotImplemented());
    }
  deallog << file2_name << std::endl;

  data_out.build_patches (2);
  data_out.write_gnuplot(deallog.get_file_stream());
//  std::ofstream file2(file2_name.c_str());
//  file2.setf(std::ios::showpoint | std::ios::left);
//  data_out.write_gnuplot(file2);
//  file2.close();
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tria;
  MappingQ1<2> mapping;

  make_grid (tria);

  deallog.push("dg2dg1_int");
  test(tria, mapping, FE_DGQ<2>(2), "dg2", FE_DGQ<2>(1), "dg1", 1);
  deallog.pop();
  deallog.push("cg2dg1_back");
  test(tria, mapping, FE_Q<2>(2),   "cg2", FE_DGQ<2>(1), "dg1", 2);
  deallog.pop();
  deallog.push("dg2dg1_back");
  test(tria, mapping, FE_DGQ<2>(2), "dg2", FE_DGQ<2>(1), "dg1", 2);
  deallog.pop();
  deallog.push("cg2dg1_diff");
  test(tria, mapping, FE_Q<2>(2),   "cg2", FE_DGQ<2>(1), "dg1", 3);
  deallog.pop();
  deallog.push("dg2dg1_diff");
  test(tria, mapping, FE_DGQ<2>(2), "dg2", FE_DGQ<2>(1), "dg1", 3);
  deallog.pop();
}
