//----------------------------  measure_et_al.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  measure_et_al.cc  ---------------------------

// Computes diameter, extent_in_direction, and minimum_vertex_distance on a variety of cells

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
	break;
	}
		case 2:
      {
	GridGenerator::hyper_cube(tria, 1., 3.);
	Point<2> &v0=tria.begin_active()->vertex(0);
	v0(0) = 0.;
 	Point<2> &v3=tria.begin_active()->vertex(3);
 	v3(0) = 4.;
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
		case 2:// like case 1
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
  for (unsigned int case_no=0; case_no<3; ++case_no)
    {
      create_triangulation(case_no, tria);
      deallog << "dim" << dim << ":case" << case_no << ":diameter="
	      << tria.begin_active()->diameter() << std::endl;
      deallog << "dim" << dim << ":case" << case_no << ":extent_in_direction="
	      << tria.begin_active()->extent_in_direction(0) << std::endl;
      deallog << "dim" << dim << ":case" << case_no << ":minimum_vertex_distance="
	      << tria.begin_active()->minimum_vertex_distance() << std::endl;
      tria.clear();
    }
}


int main()
{
  std::ofstream logfile ("measure_et_al_02/output");
  deallog << std::setprecision (PRECISION);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
}

