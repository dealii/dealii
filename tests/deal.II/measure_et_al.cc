//----------------------------  project_to_surface_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  project_to_surface_01.cc  ---------------------------

// Computes measure, center and barycenter on a variety of cells

#include "../tests.h"
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
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
//	    exact_areas.push_back(7.);
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
  std::ofstream logfile ("measure_et_al/output");
  deallog << std::setprecision (PRECISION);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
}

