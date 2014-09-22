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



// step-10 with MappingQ(8) that shows pi with 14 digits already at the first
// iteration for the volume and MappingQ(20) for the boundary part that shows
// that things still work for very high order


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
void compute_pi_by_area ()
{
  deallog << "Computation of Pi by the area:" << std::endl
          << "==============================" << std::endl;

  const unsigned int degree = 10-dim;
  deallog << "Degree = " << degree << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball (triangulation);

  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);

  MappingQ<dim>     mapping(degree);
  const FE_Q<dim>   dummy_fe (1);
  const QGauss<dim> quadrature (degree+3);

  DoFHandler<dim> dof_handler (triangulation);

  FEValues<dim> x_fe_values (mapping, dummy_fe, quadrature,
                             update_JxW_values);

  // in 3D, we obtain many digits only on a finer mesh
  if (dim == 3)
    triangulation.refine_global(1);
  for (int refinement=0; refinement<4-dim;
       ++refinement, triangulation.refine_global (1))
    {
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

      // As a variation from step-10, no convergence table here because we
      // should be in the regime of roundoff errors where comparing rates does
      // not make sense as it is not necessarily stable between
      // systems. Rather, check that we have at least 14 digits of pi correct
      // on any level.
      if (dim == 2)
        {
          deallog << "Evaluation of pi on " << triangulation.n_active_cells()
                  << " cells: " << area << std::endl;
          // assert accuracy because numdiff might cut off digits from output
          Assert(std::abs(area - pi) < 1e-14,
                 ExcMessage("Calculation not accurate"));
        }
      else
        {
          area *= 0.75;
          deallog << "Evaluation of pi on in 3D " << triangulation.n_active_cells()
                  << " cells: " << area << std::endl;
          //Assert(std::abs(area - pi) < 1e-12,
          //       ExcMessage("Calculation not accurate"));
        }
    };
  deallog << std::endl;
}



template <int dim>
void compute_pi_by_perimeter ()
{
  deallog << "Computation of Pi by the perimeter:" << std::endl
          << "===================================" << std::endl;


  const unsigned int degree = 20;
  deallog << "Degree = " << degree << std::endl;
  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball (triangulation);

  static const HyperBallBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);

  const MappingQ<dim> mapping (degree);
  const FE_Q<dim>     fe (1);
  const QGauss<dim-1> quadrature (degree);

  DoFHandler<dim> dof_handler (triangulation);

  FEFaceValues<dim> x_fe_face_values (mapping, fe, quadrature,
                                      update_JxW_values);
  for (unsigned int refinement=0; refinement<2;
       ++refinement, triangulation.refine_global (1))
    {
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
      deallog << "Evaluation of pi on " << triangulation.n_active_cells()
              << " cells: " << perimeter/2. << std::endl;
      Assert(std::abs(perimeter/2. - pi) < 1e-14,
             ExcMessage("Calculation not accurate"));
    };

  deallog << std::endl;
}


int main ()
{
  deallog << std::setprecision(16);
  logfile << std::setprecision(16);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  compute_pi_by_area<2> ();
  compute_pi_by_perimeter<2> ();

  compute_pi_by_area<3> ();

  return 0;
}
