#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// Regression test for the face numbers of mixed meshes. A recent commit
// (f35f629b6cd) exposed an inconsistency between the way this was done via
// connectivity.h vs. ReferenceCell. This test prints out the faces (via
// GridOut) which in turn verifies that we consistently (in the sense of
// versions of deal.II, not in the sense of orientations) order the faces of
// wedges and pyramids.

int
main(int argc, char **argv)
{
  initlog();

  std::vector<Point<3>> vertices(13);
  vertices[0]  = Point<3>(0, 0, 0);
  vertices[1]  = Point<3>(3, 0, 0);
  vertices[2]  = Point<3>(3, 3, 0);
  vertices[3]  = Point<3>(0, 3, 0);
  vertices[4]  = Point<3>(0, 0, 3);
  vertices[5]  = Point<3>(3, 0, 3);
  vertices[6]  = Point<3>(3, 3, 3);
  vertices[7]  = Point<3>(0, 3, 3);
  vertices[8]  = Point<3>(1, 1, 5);
  vertices[9]  = Point<3>(3, 1, 5);
  vertices[10] = Point<3>(8, 0, 3);
  vertices[11] = Point<3>(8, 3, 3);
  vertices[12] = Point<3>(8, 1, 5);

  std::vector<CellData<3>> cells(4);
  cells[0].vertices = {0, 1, 3, 2, 4, 5, 7, 6};
  cells[1].vertices = {4, 5, 7, 6, 8};
  cells[2].vertices = {5, 9, 6, 8};
  cells[3].vertices = {11, 10, 12, 6, 5, 9};

  std::vector<CellData<2>> quads(14);
  quads[0].boundary_id  = 2;
  quads[0].vertices     = {0, 1, 3, 2};
  quads[1].boundary_id  = 5;
  quads[1].vertices     = {0, 3, 4, 7};
  quads[2].boundary_id  = 1;
  quads[2].vertices     = {0, 4, 1, 5};
  quads[3].boundary_id  = 3;
  quads[3].vertices     = {1, 2, 5, 6};
  quads[4].boundary_id  = 4;
  quads[4].vertices     = {3, 7, 2, 6};
  quads[5].boundary_id  = 5;
  quads[5].vertices     = {4, 7, 8};
  quads[6].boundary_id  = 1;
  quads[6].vertices     = {5, 4, 8};
  quads[7].boundary_id  = 3;
  quads[7].vertices     = {6, 9, 8};
  quads[8].boundary_id  = 4;
  quads[8].vertices     = {7, 6, 8};
  quads[9].boundary_id  = 1;
  quads[9].vertices     = {9, 5, 8};
  quads[10].boundary_id = 6;
  quads[10].vertices    = {10, 11, 12};
  quads[11].boundary_id = 1;
  quads[11].vertices    = {10, 12, 5, 9};
  quads[12].boundary_id = 4;
  quads[12].vertices    = {11, 10, 6, 5};
  quads[13].boundary_id = 3;
  quads[13].vertices    = {12, 11, 9, 6};

  SubCellData subcell_data;
  subcell_data.boundary_quads = quads;

  Triangulation<3> tria;
  tria.create_triangulation(vertices, cells, subcell_data);
  GridOut().write_vtk(tria, deallog.get_file_stream());
}
