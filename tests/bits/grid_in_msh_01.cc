//----------------------------  grid_in_msh_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors and Luca Heltai
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_in_msh_01.cc  ---------------------------


// check whether we can read in with the gmsh format

#include "../tests.h"
#include <base/logstream.h>

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_in.h>
  
#include <fstream>
#include <cmath>
  

template <int dim>
void gmsh_grid (const char *name)
{
  Triangulation<dim> tria;
  GridIn<dim> grid_in;
  grid_in.attach_triangulation (tria);
  std::ifstream input_file(name);
  grid_in.read_msh(input_file);
  
  deallog << "  " << tria.n_active_cells() << " active cells" << std::endl;

  int hash = 0;
  int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator c=tria.begin_active();
       c!=tria.end(); ++c, ++index)
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells()+1);
  deallog << "  hash=" << hash << std::endl;
}


int main () 
{
  std::ofstream logfile("grid_in_msh_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
				
  try
    {
      gmsh_grid<2> ("grid_in_msh_01.2d.msh");
      gmsh_grid<2> ("grid_in_msh_01.2da.msh");
      gmsh_grid<3> ("grid_in_msh_01.3d.msh");
      gmsh_grid<3> ("grid_in_msh_01.3da.msh");
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
  
  return 0;
}
