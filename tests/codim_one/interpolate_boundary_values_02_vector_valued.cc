
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


// test VectorTools::interpolate_boundary_values for codim=1. like _02
// but for vector-valued elements

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>
#include <base/function_lib.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <numerics/vectors.h>

#include <string>

std::ofstream logfile("interpolate_boundary_values_02_vector_valued/output");

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

void test() {
  const int dim = 1;
  const int spacedim = 2;

    Triangulation<dim, spacedim> tria;
    GridGenerator::hyper_cube (tria);
    deallog << tria.n_active_cells() << " active cells" << std::endl;

    FESystem<dim,spacedim> fe(FE_Q<dim,spacedim>(2), spacedim);
    DoFHandler<dim,spacedim> dof_handler (tria);
    dof_handler.distribute_dofs (fe);

    deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;

				     // test left and right boundary
				     // separatel
    for (unsigned int boundary_id=0; boundary_id<2; ++boundary_id)
      {
	std::map<unsigned int, double> bv;
	VectorTools::interpolate_boundary_values (dof_handler,
						  boundary_id,
						  X<spacedim>(),
						  bv);
	deallog << bv.size() << " boundary degrees of freedom" << std::endl;

	for (std::map<unsigned int, double>::const_iterator i = bv.begin();
	     i != bv.end(); ++i)
	  deallog << i->first << ' ' << i->second << std::endl;

	for (DoFHandler<dim,spacedim>::active_cell_iterator
	       cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
	  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	    if (cell->at_boundary(f) &&
		(cell->face(f)->boundary_indicator() == boundary_id))
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
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test();

  return 0;
}

