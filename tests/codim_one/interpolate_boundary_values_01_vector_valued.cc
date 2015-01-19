// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// test VectorTools::interpolate_boundary_values for codim=1. like
// _01_vector_valued, but for vector-valued functions

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>

#include <string>

std::ofstream logfile("output");

template <int dim>
class X : public Function<dim>
{
public:
  X() : Function<dim>(dim) {}

  double value (const Point<dim> &p,
                const unsigned int component) const
  {
    return p[component];
  }
};


template <int dim, int spacedim>
void test(std::string filename)
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  deallog << tria.n_active_cells() << " active cells" << std::endl;

  FESystem<dim,spacedim> fe(FE_Q<dim,spacedim> (2), spacedim);
  DoFHandler<dim,spacedim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;

  std::map<types::global_dof_index, double> bv;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            X<spacedim>(),
                                            bv);
  deallog << bv.size() << " boundary degrees of freedom" << std::endl;

  for (std::map<types::global_dof_index, double>::const_iterator i = bv.begin();
       i != bv.end(); ++i)
    deallog << i->first << ' ' << i->second << std::endl;

  for (typename DoFHandler<dim,spacedim>::active_cell_iterator
       cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->at_boundary(f))
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
          for (unsigned int i=0; i<fe.dofs_per_vertex; ++i)
            {
              Assert (bv.find(cell->face(f)->vertex_dof_index(v,i))
                      != bv.end(),
                      ExcInternalError());
              Assert (bv[cell->face(f)->vertex_dof_index(v,i)]
                      ==
                      X<spacedim>()
                      .value(cell->face(f)->vertex(v),i),
                      ExcInternalError());
            }
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2,3>(SOURCE_DIR "/grids/square.inp");
  test<2,3>(SOURCE_DIR "/grids/sphere_1.inp");

  return 0;
}

