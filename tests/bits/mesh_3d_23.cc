//----------------------------  mesh_3d_23.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2008 by Timo Heister and the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_23.cc  ---------------------------


// like mesh_3d_22, but further reduced: when creating output with
// DataOut using a MappingQ(3) on a mesh with flipped cells, we get
// bogus output at the interior face.

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <fstream>
#include <cmath>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <vector>


int main ()
{
  std::ofstream logfile ("mesh_3d_23/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.0e-10);

  Triangulation<3> triangulation;
  GridIn<3> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream inputStream("mesh_3d_22/mesh.msh");
  grid_in.read_msh (inputStream);

  MappingQ<3> mapping(3);
  FE_Q<3> fe(3);
  DoFHandler<3> dofh(triangulation);

  dofh.distribute_dofs(fe);

  Vector<double> x(dofh.n_dofs());
  DataOut<3> data_out;

  data_out.attach_dof_handler(dofh);
  data_out.add_data_vector(x, "u");
  data_out.build_patches(mapping, 3);

  data_out.write_gnuplot (deallog.get_file_stream());
}

