#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/geometry_info.h>

#include <fstream>

#include <bitset>

using namespace dealii;

//
// Test GeometryInfo<dim>::face_to_cell_vertices
// for correct behaviour under face_orientation face_flip and face_rotation
//


template<int dim>
void test_vertices()
{
  deallog << dim << "D:" << std::endl;

  for(unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i) {

    deallog << "face " << i << ":" << std::endl;

    for(unsigned int o = 0; o < 8; ++o) {
      const std::bitset<3> orientation = o;

      deallog << "orientation " << orientation[0]
              << ", flip " << orientation[1]
              << ", rotation " << orientation[2]
              << ":" << std::endl << "    ";

      for(unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_face; ++j) {
        deallog << " (" << j << " -> "
                << GeometryInfo<dim>::face_to_cell_vertices(i, j, orientation[0], orientation[1], orientation[2])
                << " )";
      }
      deallog << std::endl;
    }
  }
}


template<int dim>
void test_lines()
{
  deallog << dim << "D:" << std::endl;

  for(unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i) {

    deallog << "face " << i << ":" << std::endl;

    for(unsigned int o = 0; o < 8; ++o) {
      const std::bitset<3> orientation = o;

      deallog << "orientation " << orientation[0]
              << ", flip " << orientation[1]
              << ", rotation " << orientation[2]
              << ":" << std::endl << "    ";

      for(unsigned int j = 0; j < GeometryInfo<dim>::lines_per_face; ++j) {
        deallog << " (" << j << " -> "
                << GeometryInfo<dim>::face_to_cell_lines(i, j, orientation[0], orientation[1], orientation[2])
                << " )";
      }
      deallog << std::endl;
    }
  }
}


int main()
{
  std::ofstream logfile("geometry_info_8/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog << std::endl << "GeometryInfo<dim>::face_to_cell_vertices:" << std::endl;

  test_vertices<1>();
  test_vertices<2>();
  test_vertices<3>();

  deallog << std::endl << std::endl << "GeometryInfo<dim>::face_to_cell_lines:" << std::endl;

  test_lines<2>();
  test_lines<3>();

  return 0;
}

