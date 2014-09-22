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


// Computes measure, center and barycenter on a variety of cells

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <fstream>
#include <iomanip>

#define PRECISION 5


template<int dim>
void create_triangulation(const unsigned int,
                          Triangulation<dim> &)
{
  Assert(false, ExcNotImplemented());
}


template<>
void create_triangulation(const unsigned int case_no,
                          Triangulation<2> &tria)
{
  switch (case_no)
    {
    case 0:
      GridGenerator::hyper_cube(tria, 1., 3.);
      break;
    case 1:
    {
      GridGenerator::hyper_cube(tria, 1., 3.);
      Point<2> &v0=tria.begin_active()->vertex(0);
      v0(0) = 0.;
      Point<2> &v2=tria.begin_active()->vertex(3);
      v2(0) = 5.;
      v2(1) = 4.;
//      exact_areas.push_back(7.);
      break;
    }
    default:
      Assert(false, ExcNotImplemented());
    };
}


template<>
void create_triangulation(const unsigned int case_no,
                          Triangulation<3> &tria)
{
  switch (case_no)
    {
    case 0:
      GridGenerator::hyper_cube(tria, 1., 3.);
      break;
    case 1:
    {
      GridGenerator::hyper_cube(tria, 1., 3.);
      Point<3> &v0=tria.begin_active()->vertex(0);
      v0(0) = 0.;
      break;
    }
    default:
      Assert(false, ExcNotImplemented());
    };
}


template<int dim>
void test()
{
  Triangulation<dim> tria;
  for (unsigned int case_no=0; case_no<2; ++case_no)
    {
      create_triangulation(case_no, tria);
      deallog << "dim" << dim << ":case" << case_no << ":measure="
              << tria.begin_active()->measure() << std::endl;
      deallog << "dim" << dim << ":case" << case_no << ":center="
              << tria.begin_active()->center() << std::endl;
      deallog << "dim" << dim << ":case" << case_no << ":barycenter="
              << tria.begin_active()->barycenter() << std::endl;
      tria.clear();
    }
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (PRECISION);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
}

