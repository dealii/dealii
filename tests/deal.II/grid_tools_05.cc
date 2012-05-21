//----------------------------  grid_tools.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_tools.cc  ---------------------------


// check GridTools::collect_periodic_cell_pairs for some simple cases


#include "../tests.h"
#include <deal.II/base/logstream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <fstream>

std::ofstream logfile("grid_tools_05/output");

using namespace dealii;

template<int dim>
void test1 ()
{
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation, 0., 1.);
  triangulation.refine_global(2);

  deallog << std::endl << std::endl << "Hyper cube, dimension " << dim << ":" << std::endl;

  typedef typename Triangulation<dim>::cell_iterator CELL_IT;
  typedef typename std::map<CELL_IT, CELL_IT> MAP;

  MAP matched_cells = GridTools::collect_periodic_cell_pairs(triangulation, 0, /* direction: */1);

  deallog << "Matching:" << std::endl;
  for (typename MAP::const_iterator it = matched_cells.begin(); it != matched_cells.end(); ++it) {
    deallog << it->first << "  -  " << it->second << std::endl;
  }
}

template<int dim>
void test2 ()
{
  /* Unfortunately, there is no parallelogram for 3D.. */
  Assert(dim == 2,
         ExcNotImplemented());

  Triangulation<dim> triangulation;

  Tensor<2,dim> vectors;
  vectors[0][0] = 1.;
  vectors[0][1] = 0.;
  vectors[1][0] = 0.5;
  vectors[1][1] = 1.;

  Tensor<1,dim> offset;
  offset[0] = 0.5;
  offset[1] = 1.;

  GridGenerator::parallelogram(triangulation, vectors);
  triangulation.refine_global(3);

  deallog << std::endl << std::endl << "Parallelogram, dimension " << dim << ":" << std::endl;


  typedef typename Triangulation<dim>::cell_iterator CELL_IT;
  typedef typename std::map<CELL_IT, CELL_IT> MAP;

  MAP matched_cells = GridTools::collect_periodic_cell_pairs(triangulation, 0, /* direction: */0);
  deallog << "Matching, direction 0:" << std::endl;
  for (typename MAP::const_iterator it = matched_cells.begin(); it != matched_cells.end(); ++it) {
    deallog << it->first << "  -  " << it->second << std::endl;
  }

  MAP matched_cells_2 = GridTools::collect_periodic_cell_pairs(triangulation, 0, /* direction: */1, offset);
  deallog << "Matching, direction 1:" << std::endl;
  for (typename MAP::const_iterator it = matched_cells_2.begin(); it != matched_cells_2.end(); ++it) {
    deallog << it->first << "  -  " << it->second << std::endl;
  }
}

int main()
{
  deallog << std::setprecision(4);
  logfile << std::setprecision(4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  /* Let's try to match the hyper_cube in 3D:*/
  test1<3> ();

  /* Let's try to match a parallelogramm in 2D */
  test2<2> ();

  return 0;
}
