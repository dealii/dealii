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



// a un-hp-ified version of hp/step-10


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <fstream>
std::ofstream logfile("output");


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/fe/mapping_q.h>

#include <iomanip>
#include <fstream>
#include <cmath>

const long double pi = 3.141592653589793238462643;



template <int dim>
void gnuplot_output()
{
  deallog << "Output of grids into gnuplot files:" << std::endl
          << "===================================" << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball (triangulation);
  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);

  for (unsigned int refinement=0; refinement<2;
       ++refinement, triangulation.refine_global(1))
    {
      deallog << "Refinement level: " << refinement << std::endl;

      std::string filename_base = "ball";
      filename_base += '0'+refinement;

      for (unsigned int degree=1; degree<4; ++degree)
        {
          deallog << "Degree = " << degree << std::endl;

          const MappingQ<dim> mapping (degree);


          GridOut grid_out;
          GridOutFlags::Gnuplot gnuplot_flags(false, 30);
          grid_out.set_flags(gnuplot_flags);

          grid_out.write_gnuplot (triangulation, deallog.get_file_stream(), &mapping);
        }
      deallog << std::endl;
    }
}

template <int dim>
void compute_pi_by_area ()
{
  deallog << "Computation of Pi by the area:" << std::endl
          << "==============================" << std::endl;

  const QGauss<dim> quadrature (4);

  for (unsigned int degree=1; degree<5; ++degree)
    {
      deallog << "Degree = " << degree << std::endl;

      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball (triangulation);

      static const HyperBallBoundary<dim> boundary;
      triangulation.set_boundary (0, boundary);

      MappingQ<dim> mapping(degree);

      const FE_Q<dim>     dummy_fe (1);

      DoFHandler<dim> dof_handler (triangulation);

      FEValues<dim> x_fe_values (mapping, dummy_fe, quadrature,
                                 update_JxW_values);

      ConvergenceTable table;


      for (unsigned int refinement=0; refinement < (degree != 4 ? 6 : 4);
           ++refinement, triangulation.refine_global (1))
        {
          table.add_value("cells", triangulation.n_active_cells());

          dof_handler.distribute_dofs (dummy_fe);

          long double area = 0;

          typename DoFHandler<dim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
          for (; cell!=endc; ++cell)
            {
              x_fe_values.reinit (cell);
              const FEValues<dim> &fe_values = x_fe_values.get_present_fe_values();
              for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
                area += fe_values.JxW (i);
            };

          table.add_value("eval.pi", static_cast<double> (area));
          table.add_value("error",   static_cast<double> (std::fabs(area-pi)));
        };

      table.omit_column_from_convergence_rate_evaluation("cells");
      table.omit_column_from_convergence_rate_evaluation("eval.pi");
      table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);

      table.set_precision("eval.pi", 16);
      table.set_scientific("error", true);

      table.write_text(deallog.get_file_stream());

      deallog << std::endl;
    };
}


template <int dim>
void compute_pi_by_perimeter ()
{
  deallog << "Computation of Pi by the perimeter:" << std::endl
          << "===================================" << std::endl;

  const QGauss<dim-1> quadrature (4);

  for (unsigned int degree=1; degree<5; ++degree)
    {
      deallog << "Degree = " << degree << std::endl;
      Triangulation<dim> triangulation;
      GridGenerator::hyper_ball (triangulation);

      static const HyperBallBoundary<dim> boundary;
      triangulation.set_boundary (0, boundary);

      const MappingQ<dim> mapping (degree);
      const FE_Q<dim>     fe (1);

      DoFHandler<dim> dof_handler (triangulation);

      FEFaceValues<dim> x_fe_face_values (mapping, fe, quadrature,
                                          update_JxW_values);
      ConvergenceTable table;

      for (unsigned int refinement=0; refinement < (degree != 4 ? 6 : 4);
           ++refinement, triangulation.refine_global (1))
        {
          table.add_value("cells", triangulation.n_active_cells());

          dof_handler.distribute_dofs (fe);

          typename DoFHandler<dim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
          long double perimeter = 0;
          for (; cell!=endc; ++cell)
            for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
              if (cell->face(face_no)->at_boundary())
                {
                  x_fe_face_values.reinit (cell, face_no);
                  const FEFaceValues<dim> &fe_face_values
                    = x_fe_face_values.get_present_fe_values ();

                  for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
                    perimeter += fe_face_values.JxW (i);
                };
          table.add_value("eval.pi", static_cast<double> (perimeter/2.));
          table.add_value("error",   static_cast<double> (std::fabs(perimeter/2.-pi)));
        };

      table.omit_column_from_convergence_rate_evaluation("cells");
      table.omit_column_from_convergence_rate_evaluation("eval.pi");
      table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);

      table.set_precision("eval.pi", 16);
      table.set_scientific("error", true);

      table.write_text(deallog.get_file_stream());

      deallog << std::endl;
    };
}


int main ()
{
  deallog << std::setprecision(6);
  logfile << std::setprecision(6);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  gnuplot_output<2>();

  compute_pi_by_area<2> ();
  compute_pi_by_perimeter<2> ();

  return 0;
}
