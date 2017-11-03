//----------------------------  manifold_id_01.cc  ---------------------------
//    Copyright (C) 2011 - 2015 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  flat_manifold_01.cc  ---------------------------


// Test that the flat manifold does what it should

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  std::vector<Point<spacedim> > vertices (GeometryInfo<dim>::vertices_per_cell);
  std::vector<CellData<dim> > cells (1);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    cells[0].vertices[i] = i;
  cells[0].material_id = 0;

  vertices[0] = Point<dim>(0,0,0);
  vertices[1] = Point<dim>(1,0,0);
  vertices[2] = Point<dim>(0.5,0.4,0);
  vertices[3] = Point<dim>(1.5,0.4,0);

  vertices[4] = Point<dim>(0,0,1);
  vertices[5] = Point<dim>(1,0,1);
  vertices[6] = Point<dim>(0.5,0.4,1);
  vertices[7] = Point<dim>(1.5,0.4,1);

  Triangulation<dim,spacedim> tria;
  tria.create_triangulation (vertices, cells, SubCellData());

  typename Triangulation<dim,spacedim>::active_cell_iterator
  cell = tria.begin_active();

  Point<dim> p1 (0.5,0,0);
  deallog << "Normal vector of face 4: " << cell->get_manifold().normal_vector(cell->face(4),p1) << std::endl;
  deallog << "Center of face 4: " << cell->face(4)->center() << std::endl;
}

int main ()
{
  initlog();

  test<3,3>();

  return 0;
}

