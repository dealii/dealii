//----------------------------  distorted_cells_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  distorted_cells_01.cc  ---------------------------


// check that indeed Triangulation::create_triangulation throws an
// exception if we have distorted cells

#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <fstream>


// create a (i) pinched cell (where two vertices coincide), or (ii)
// twisted cell (where two vertices are swapped)
template <int dim>
void check (const unsigned int testcase)
{
  std::vector<Point<dim> > vertices;
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
    vertices.push_back (GeometryInfo<dim>::unit_cell_vertex(v));

  switch (testcase)
    {
      case 1:
	    deallog << "Pinched cell in " << dim << "d" << std::endl;
	    vertices[0] = vertices[1];
	    break;
      case 2:
	    deallog << "Twisted cell in " << dim << "d" << std::endl;
	    std::swap (vertices[0], vertices[1]);
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    }
  

  std::vector<CellData<dim> > cells;
  {
    CellData<dim> cell;
    for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_cell; ++j)
      cell.vertices[j]   = j;
    cells.push_back (cell);
  }

  Triangulation<dim> coarse_grid;

  bool flag = false;
  try
    {
      coarse_grid.create_triangulation (vertices, cells, SubCellData());
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      flag = true;
      
      deallog << dcv.distorted_cells.size() << " distorted cells"
	      << std::endl;
      Assert (dcv.distorted_cells.front() == coarse_grid.begin(0),
	      ExcInternalError());
    }

  Assert (flag == true, ExcInternalError());
}


int main () 
{
  std::ofstream logfile("distorted_cells_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int testcase=1; testcase<=2; ++testcase)
    {
      check<1> (testcase);
      check<2> (testcase);
      check<3> (testcase);
    }
}

  
  
