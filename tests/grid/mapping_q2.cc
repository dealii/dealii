// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// read a file in the MSH format with quadratic elements
// created by the GMSH program. Use the additional support
// points to create a quadratic mapping.

#include "../../include/deal.II/fe/mapping_q2.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_manifold.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <istream>
#include <string>

#include "../tests.h"

template <int dim, int spacedim = dim>
void
check_file(const std::string file_name)
{
  // create triangulation
  Triangulation<dim, spacedim> tria;
  // create GridIn object and attach triangulation
  GridIn<dim, spacedim> gi;
  gi.attach_triangulation(tria);

  // vector storing the support points of the quadratic elements
  std::vector<std::vector<Point<spacedim>>> support_points;

  // read input file, extract and store support points
  std::filebuf fb;
  if (fb.open(file_name, std::ios::in))
    {
      std::istream is(&fb);
      std::string  file_ending =
        file_name.substr(file_name.find_last_of(".") + 1);
      if (file_ending == "msh")
        gi.read_msh(is, support_points);
      else
        AssertThrow(false, ExcNotImplemented());
    }

  // create different mappings for comparison
  // MappingGeneric with polynomial degree 1
  MappingQGeneric<dim, spacedim> mapping_1(1);
  // New MappingQ2 (always with polynomial degree 2) and
  // additional input parameter support_points
  MappingQ2<dim, spacedim> mapping_2(support_points);

  // create quadrature rules. Use GaussLobatto since the quadrature
  // points correspond to the support points
  QGaussLobatto<dim> quad_1(2);
  QGaussLobatto<dim> quad_2(3);

  // create finite elements
  FE_Q<dim, spacedim> fe_1(1);
  FE_Q<dim, spacedim> fe_2(2);

  // set update flags
  UpdateFlags flags = update_quadrature_points;

  // create different FEValues objects with the above mappings,
  // finite elements, and quadrature rules
  FEValues<dim, spacedim> fe_values_1(mapping_1, fe_1, quad_1, flags);
  FEValues<dim, spacedim> fe_values_2(mapping_1, fe_2, quad_2, flags);
  FEValues<dim, spacedim> fe_values_3(mapping_2, fe_2, quad_2, flags);

  for (auto cell : tria.active_cell_iterators())
    {
      // print vertices
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        deallog << cell->vertex(v) << std::endl;
      deallog << std::endl;

      // print quadrature points, should be the same results as
      // above from print vertices
      deallog << "MappingQ(1)" << std::endl;
      fe_values_1.reinit(cell);
      for (auto q : fe_values_1.get_quadrature_points())
        deallog << q << std::endl;
      deallog << std::endl;

      // print quadrature points, since quadrature rule has degree 3
      // also support points are computed and printed.
      deallog << "MappingQ(2)" << std::endl;
      fe_values_2.reinit(cell);
      for (auto q : fe_values_2.get_quadrature_points())
        deallog << q << std::endl;
      deallog << std::endl;

      // print quadrature points, this time the support points should
      // be at the same location as specified in the input file, the
      // ordering will be different though.
      deallog << "MappingQ2()" << std::endl;
      fe_values_3.reinit(cell);
      for (auto q : fe_values_3.get_quadrature_points())
        deallog << q << std::endl;
      deallog << std::endl;
    }
  deallog << "END OF TEST" << std::endl;
  deallog << std::endl;
}

int
main()
{
  initlog();
  // QUAD9 ELEMENTS
  // one quadrangle with 9 nodes, dim = spacedim = 2, Gmsh file version 2
  deallog.push("quad9_1ele_dim2_spacedim2");
  check_file<2, 2>(std::string(SOURCE_DIR "/mapping_q2/quad9_1ele.msh"));
  deallog.pop();

  // one quadrangle with 9 nodes, dim = 2 and spacedim = 3
  deallog.push("quad9_1ele_dim2_spacedim3");
  check_file<2, 3>(std::string(SOURCE_DIR "/mapping_q2/quad9_1ele.msh"));
  deallog.pop();

  // four 9-noded quadrangles, dim = spacedim = 2, Gmsh file version 4
  deallog.push("quad9_4ele_dim2_spacedim2_v4");
  check_file<2, 2>(std::string(SOURCE_DIR "/mapping_q2/quad9_4ele_v4.msh"));
  deallog.pop();

  // four 9-noded quadrangles, dim = 2 and spacedim = 3, Gmsh file version 4
  deallog.push("quad9_4ele_dim2_spacedim3_v4");
  check_file<2, 3>(std::string(SOURCE_DIR "/mapping_q2/quad9_4ele_v4.msh"));
  deallog.pop();



  // HEX27 ELEMENTS
  // one hexadron with 27 nodes, dim = spacedim = 3, Gmsh file version 2
  deallog.push("hex27_1ele");
  check_file<3, 3>(std::string(SOURCE_DIR "/mapping_q2/hex27_1ele.msh"));
  deallog.pop();

  // 8 hexadra with 27 nodes, dim = spacedim = 3, Gmsh file version 2
  deallog.push("hex27_8ele_v2");
  check_file<3, 3>(std::string(SOURCE_DIR "/mapping_q2/hex27_8ele_v2.msh"));
  deallog.pop();
}
