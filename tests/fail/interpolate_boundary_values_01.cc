
//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// test VectorTools::interpolate_boundary_values for codim=1

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>
#include <base/function_lib.h>
#include <grid/tria.h>
#include <grid/grid_in.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <numerics/vectors.h>

#include <string>

std::ofstream logfile("interpolate_boundary_values_01/output");

template <int dim, int spacedim>
void test(std::string filename) {
    Triangulation<dim, spacedim> tria;
    GridIn<dim, spacedim> gi;
    gi.attach_triangulation (tria);
    std::ifstream in (filename.c_str());
    gi.read_ucd (in);

    deallog << tria.n_active_cells() << " active cells" << std::endl;

    FE_Q<dim,spacedim> fe(2);
    DoFHandler<dim,spacedim> dof_handler (tria);
    dof_handler.distribute_dofs (fe);

    deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;

    std::map<unsigned int, double> bv;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      Functions::SquareFunction<spacedim>(),
					      bv);
    deallog << bv.size() << " boundary degrees of freedom" << std::endl;

    for (std::map<unsigned int, double>::const_iterator i = bv.begin();
	 i != bv.end(); ++i)
      deallog << i->first << ' ' << i->second << std::endl;
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2,3>("grids/square.inp");
  test<2,3>("grids/sphere_1.inp");

  return 0;
}

