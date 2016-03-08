//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2015 by the deal.II authors
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read goteborg.iges and dump its topological structure to the
// logfile.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/opencascade/utilities.h>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>
#include <Standard_Stream.hxx>

using namespace OpenCASCADE;

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  std::vector<TopoDS_Face> faces;
  std::vector<TopoDS_Edge> edges;
  std::vector<TopoDS_Vertex> vertices;

  TopoDS_Shape sh = read_STEP(SOURCE_DIR "/step_files/goteborg.step");
  extract_geometrical_shapes(sh, faces, edges, vertices);

  dealii::OpenCASCADE::NURBSPatchManifold<2,3> manifold(faces[0]);

  auto bounds = manifold.get_uv_bounds();
  auto u_min = std::get<0>(bounds);
  auto u_max = std::get<1>(bounds);
  auto v_min = std::get<2>(bounds);
  auto v_max = std::get<3>(bounds);  
  
  deallog << " u_min = " << std::get<0>(bounds) << std::endl;
  deallog << " u_max = " << std::get<1>(bounds) << std::endl;
  deallog << " v_min = " << std::get<2>(bounds) << std::endl;
  deallog << " v_max = " << std::get<3>(bounds) << std::endl;

  deallog << "=======================================" << std::endl;
  int len = 10;
  for(unsigned int i = 0; i <= 10; ++i)
  {
    double step = ((double)i)/10;
    deallog << " pos = " <<  step << std::endl;
    double u_pos = step * u_min + (1-step) * u_max;
    double v_pos = step * v_min + (1-step) * v_max;
    Point<2> uv(u_pos, v_pos);
    deallog << " uv = " << uv << std::endl;
    auto q = manifold.push_forward(uv);
    deallog << " q = " << q << std::endl;
    auto uv_ = manifold.pull_back(q);
    deallog << " uv = " << uv << std::endl;
    deallog << "=======================================" << std::endl;
  }

  double u_pos =  (u_min + u_max)/2;
  double v_pos =  (v_min + v_max)/2;
  Point<2> uv(u_pos, v_pos);
  auto D = manifold.push_forward_gradient(uv);
  deallog << "=======================================" << std::endl;
  deallog << " | " <<  D[1][0] << " | " <<  D[0][0] <<" | " << std::endl;
  deallog << " | " <<  D[1][1] << " | " <<  D[0][1] <<" | " << std::endl;
  deallog << " | " <<  D[1][2] << " | " <<  D[0][2] <<" | " << std::endl;
  deallog << "=======================================" << std::endl;

  return 0;
}

