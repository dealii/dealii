// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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


// compute some convergence results from computing pi on a mesh that
// is deformed to represent a quarter of a ring

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <fstream>
#include <iostream>

#include <deal.II/fe/mapping_q_eulerian.h>


// .... IMPOSED DISPLACEMENT

template <int dim>
class ImposedDisplacement : public Function<dim>
{
public:
  ImposedDisplacement() : Function<dim> (dim) { }
  virtual void vector_value(const Point<dim> &p,
                            Vector<double> &value) const;
};

template <>
void ImposedDisplacement<2>::vector_value(const Point<2> &p,
                                          Vector<double> &value) const
{
  double radius = 1 + (sqrt(5)-1)*p(0);
  double angle  = 0.5*numbers::PI*(1-p(1));
  value(0) = radius*sin(angle)-p(0);
  value(1) = radius*cos(angle)-p(1);
}


// .... MAPPING TEST CLASS

template <int dim>
class MappingTest
{
public:
  MappingTest (unsigned int degree);
  ~MappingTest ();

  void run_test();
  void graphical_output();

private:
  double compute_area();
  void explicitly_move_mesh();
  void write_tria_to_eps(std::string id);

  Triangulation<dim>     triangulation;
  DoFHandler<dim>        dof_handler;
  FESystem<dim>          fe;

  unsigned int           degree;

  ImposedDisplacement<dim> imposed_displacement;
  Vector<double>           displacements;
};


// .... CONSTRUCTOR

template <int dim>
MappingTest<dim>::MappingTest (unsigned int degree)
  :
  dof_handler (triangulation),
  fe (FE_Q<dim>(degree),dim),
  degree(degree)
{ }


// .... DESTRUCTOR

template <int dim>
MappingTest<dim>::~MappingTest ()
{
  dof_handler.clear ();
}


// .... COMPUTE AREA

template <int dim>
double MappingTest<dim>::compute_area ()
{
  QGauss<dim>  quadrature_formula(degree+1);

  MappingQEulerian<dim> mapping(degree,displacements,dof_handler);

  FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                           update_JxW_values);

  const unsigned int   n_q_points = quadrature_formula.size();

  long double area = 0.;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      for (unsigned int q=0; q<n_q_points; ++q) area += fe_values.JxW(q);
    }

  return area;
}


// .... RUN TEST

template <int dim>
void MappingTest<dim>::run_test ()
{
  GridGenerator::hyper_cube (triangulation,0, 1);

  ConvergenceTable table;

  for (unsigned int ref_level = 0;
       ref_level < 5;
       ++ref_level, triangulation.refine_global(1))
    {

      dof_handler.distribute_dofs (fe);
      displacements.reinit (dof_handler.n_dofs());

      VectorTools::interpolate(MappingQ1<dim>(),dof_handler,
                               imposed_displacement,displacements);


      table.add_value("cells",triangulation.n_active_cells());
      table.add_value("dofs",dof_handler.n_dofs());

      long double area  = compute_area();
      long double error = std::fabs(numbers::PI-area)/numbers::PI;

      table.add_value("area",  static_cast<double> (area));
      table.add_value("error", static_cast<double> (error));
    }

  table.set_precision("area", 8);
  table.set_precision("error", 4);
  table.set_scientific("error", true);
  table.evaluate_convergence_rates("error",
                                   ConvergenceTable::reduction_rate_log2);
  table.write_text(deallog.get_file_stream());
  deallog << std::endl;

}


// .... EXPLICITLY MOVE MESH

template <int dim>
void MappingTest<dim>::explicitly_move_mesh ()
{
  std::vector<bool> moved (triangulation.n_vertices(),false);
  unsigned int vpc = GeometryInfo<dim>::vertices_per_cell;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active (),
  endc = dof_handler.end();

  for (; cell != endc; cell++)
    {
      for (unsigned int v=0; v < vpc; v++)
        {
          if (moved[cell->vertex_index(v)] == false)
            {
              moved[cell->vertex_index(v)] =  true;
              Point<dim> vertex_disp;
              for (unsigned int d=0; d<dim; d++)
                {
                  vertex_disp[d] = displacements(cell->vertex_dof_index(v,d));
                }
              cell->vertex(v) += vertex_disp;
            }
        }
    }
}



// .... GRAPHICAL OUTPUT

template <int dim>
void MappingTest<dim>::graphical_output ()
{
  GridGenerator::hyper_cube (triangulation,0, 1);
  triangulation.refine_global(4);

  dof_handler.distribute_dofs (fe);
  displacements.reinit (dof_handler.n_dofs());

  VectorTools::interpolate(MappingQ1<dim>(),dof_handler,
                           imposed_displacement,displacements);

  explicitly_move_mesh();
}


// .... MAIN

int main ()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // convergence studies

  for (unsigned int degree = 1; degree <=4; ++degree)
    {
      deallog << ".... Q" << degree << " Mapping ...." << std::endl;
      MappingTest<2> test_one(degree);
      test_one.run_test();
    }

  // graphical output

  MappingTest<2> test_two(1);
  test_two.graphical_output();

  return 0;
}

