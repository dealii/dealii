//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007 by the deal.II authors and Luca Heltai
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test output and input of meshes, to see if the meshes are different
// after a couple of global refinements

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <fstream>
#include <iomanip>


template<int dim>
void test(std::ostream& out)
{
  GridOut go;
  go.set_flags(GridOutFlags::Ucd(false, true, true));
  go.set_flags(GridOutFlags::Msh(true, true));
  Triangulation<dim> tr;
  
  GridGenerator::hyper_cube_with_cylindrical_hole(tr, .3, .4, 1, 1, false);
  CylinderBoundary<dim> boundary(.3, 2);
  tr.set_boundary(1, boundary);
  {
    std::ofstream grid_file("grid_output_input/coarse_grid.inp");
    go.write_ucd(tr, grid_file);
    grid_file.close();
  }
  tr.refine_global(1);
  deallog << "Writing refined from constructor" << std::endl;
  go.write_ucd(tr, out);
  
  tr.set_boundary(1);
  tr.clear();
  
  GridIn<dim> gi;
  gi.attach_triangulation(tr);
  {
    deallog << "Read coarse grid" << std::endl;
    std::ifstream grid_file("grid_output_input/coarse_grid.inp");
    gi.read_ucd(grid_file);
    grid_file.close();
  }

  tr.set_boundary(1, boundary);
  tr.refine_global(1);
  deallog << "Writing refined from file" << std::endl;
  go.write_ucd(tr, out);
  
  {
    std::ofstream grid_file("grid_output_input/grid.msh");
    go.write_msh(tr, grid_file);
    grid_file.close();
  }
  
  tr.set_boundary(1);
  tr.clear();
}

int main()
{
  std::ofstream logfile("grid_output_input/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();
}
