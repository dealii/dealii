//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2017 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read a file in iges format, and write it out again in the same
// format.

#include "../tests.h"
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

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
  initlog();
  std::ofstream geo("file.geo");

  geo << "Lx = 25.0;" << std::endl
      << "Ly = 1.0;" << std::endl
      << "Point(1) = {0, 0, 0, Lx};" << std::endl
      << "Point(2) = {Lx, 0, 0, Lx};" << std::endl
      << "Point(3) = {Lx, Ly, 0, Lx};" << std::endl
      << "Point(4) = {0, Ly, 0, Lx};" << std::endl
      << "Line(1) = {1, 2};" << std::endl
      << "Line(2) = {2, 3};" << std::endl
      << "Line(3) = {3, 4};" << std::endl
      << "Line(4) = {4, 1};" << std::endl
      << "Line Loop(5) = {3, 4, 1, 2};" << std::endl
      << "Plane Surface(0) = {5};" << std::endl
      << "Transfinite Surface{0} = {5};" << std::endl
      << "Recombine Surface {0};" << std::endl;

  geo.close();

  const int ierr = std::system(DEAL_II_GMSH_EXECUTABLE_PATH " -2 file.geo 1>file.log 2>file_warn.log");
  Assert(ierr==0, ExcInternalError());
  gmsh_grid<2>("file.msh");
  cat_file("file.msh");

  return 0;
}
