// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2015 by the deal.II authors
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


// Show the Jacobians, inverse Jacobians, and Jacobian gradients of
// FEFaceValues and FESubfaceValues on a Cartesian mesh with MappingCartesian
// with one quadrature point

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

template<int dim>
void test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  MappingCartesian<dim> mapping;
  FE_Nothing<dim> dummy;
  // choose a point that is not right in the middle of the cell so that the
  // Jacobian contains many nonzero entries
  Point<dim-1> quad_p;
  for (int d=0; d<dim-1; ++d)
    quad_p(d) = 0.42 + 0.11 * d;
  Quadrature<dim-1> quad(quad_p);
  FEFaceValues<dim> fe_val (mapping, dummy, quad,
                            update_jacobians | update_inverse_jacobians |
                            update_jacobian_grads);
  FESubfaceValues<dim> fe_sub_val (mapping, dummy, quad,
                                   update_jacobians | update_inverse_jacobians |
                                   update_jacobian_grads);
  deallog << dim << "d Jacobians:" << std::endl;
  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active(), endc = tria.end();
  for ( ; cell != endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        fe_val.reinit (cell, f);

        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int e=0; e<dim; ++e)
            deallog << fe_val.jacobian(0)[d][e] << " ";
        deallog << std::endl;

        // Also check the Jacobian with FESubfaceValues
        if (cell->at_boundary(f) == false &&
            cell->neighbor(f)->level() < cell->level())
          {
            fe_sub_val.reinit(cell->neighbor(f), cell->neighbor_face_no(f),
                              cell->neighbor_of_coarser_neighbor(f).second);

            for (unsigned int d=0; d<dim; ++d)
              for (unsigned int e=0; e<dim; ++e)
                deallog << fe_sub_val.jacobian(0)[d][e] << " ";
            deallog << std::endl;
          }
      }
  deallog << std::endl;

  deallog << dim << "d inverse Jacobians:" << std::endl;
  cell = tria.begin_active();
  endc = tria.end();
  for ( ; cell != endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        fe_val.reinit (cell, f);

        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int e=0; e<dim; ++e)
            deallog << fe_val.inverse_jacobian(0)[d][e] << " ";
        deallog << std::endl;

        // Also check the inverse Jacobian with FESubfaceValues
        if (cell->at_boundary(f) == false &&
            cell->neighbor(f)->level() < cell->level())
          {
            fe_sub_val.reinit(cell->neighbor(f), cell->neighbor_face_no(f),
                              cell->neighbor_of_coarser_neighbor(f).second);

            for (unsigned int d=0; d<dim; ++d)
              for (unsigned int e=0; e<dim; ++e)
                deallog << fe_sub_val.inverse_jacobian(0)[d][e] << " ";
            deallog << std::endl;
          }

    }
  deallog << std::endl;

  deallog << dim << "d Jacobian gradients:" << std::endl;
  cell = tria.begin_active();
  endc = tria.end();
  for ( ; cell != endc; ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        fe_val.reinit (cell, f);

        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int e=0; e<dim; ++e)
            for (unsigned int f=0; f<dim; ++f)
              deallog << fe_val.jacobian_grad(0)[d][e][f] << " ";
        deallog << std::endl;

        // Also check the Jacobian with FESubfaceValues
        if (cell->at_boundary(f) == false &&
            cell->neighbor(f)->level() < cell->level())
          {
            fe_sub_val.reinit(cell->neighbor(f), cell->neighbor_face_no(f),
                              cell->neighbor_of_coarser_neighbor(f).second);

            for (unsigned int d=0; d<dim; ++d)
              for (unsigned int e=0; e<dim; ++e)
                for (unsigned int f=0; f<dim; ++f)
                  deallog << fe_sub_val.jacobian_grad(0)[d][e][f] << " ";
            deallog << std::endl;
          }
      }
  deallog << std::endl;
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(8) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();

  return 0;
}
