// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Read goteborg.iges and dump its topological structure to the
// logfile.

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>

#include "../tests.h"

using namespace OpenCASCADE;

int
main()
{
  initlog();
  deallog.precision(10);

  std::vector<TopoDS_Face>   faces;
  std::vector<TopoDS_Edge>   edges;
  std::vector<TopoDS_Vertex> vertices;

  TopoDS_Shape sh = read_STEP(SOURCE_DIR "/step_files/goteborg.step");
  extract_geometrical_shapes(sh, faces, edges, vertices);

  dealii::OpenCASCADE::NURBSPatchManifold<2, 3> manifold(faces[0]);

  const double u_min = 1.96115;
  const double u_max = 10.0000;
  const double v_min = 0.00000;
  const double v_max = 11.0000;

  deallog << "=======================================" << std::endl;
  int len = 10;
  for (unsigned int i = 0; i <= 10; ++i)
    {
      double step = ((double)i) / 10;
      deallog << " pos = " << step << std::endl;
      double   u_pos = step * u_min + (1 - step) * u_max;
      double   v_pos = step * v_min + (1 - step) * v_max;
      Point<2> uv(u_pos, v_pos);
      deallog << " uv = " << uv << std::endl;
      Point<3> q = manifold.push_forward(uv);
      deallog << " q = " << q << std::endl;
      Point<2> uv_ = manifold.pull_back(q);
      deallog << " uv = " << uv << std::endl;
      deallog << "=======================================" << std::endl;
    }

  double                  u_pos = (u_min + u_max) / 2;
  double                  v_pos = (v_min + v_max) / 2;
  Point<2>                uv(u_pos, v_pos);
  DerivativeForm<1, 2, 3> D = manifold.push_forward_gradient(uv);
  deallog << "=======================================" << std::endl;
  deallog << " | " << D[0][1] << " | " << D[0][0] << " | " << std::endl;
  deallog << " | " << D[1][1] << " | " << D[1][0] << " | " << std::endl;
  deallog << " | " << D[2][1] << " | " << D[2][0] << " | " << std::endl;
  deallog << "=======================================" << std::endl;

  return 0;
}
