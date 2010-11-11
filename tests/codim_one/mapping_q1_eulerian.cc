
//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2009, 2010 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// computes points in real space starting from some quadrature points on the unit element

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>

// all include files you need here

#include <grid/tria.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <fe/mapping.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <fe/mapping_q1_eulerian.h>
#include <base/quadrature_lib.h>

#include <fstream>
#include <string>

std::ofstream logfile("mapping_q1_eulerian/output");

template <int dim, int spacedim>
void test(std::string filename) {
  Triangulation<dim, spacedim> tria;
  FE_Q<dim, spacedim> base_fe(1);
  FESystem<dim, spacedim> fe(base_fe, spacedim);

  DoFHandler<dim, spacedim> shift_dh(tria);
  
  GridIn<dim, spacedim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  shift_dh.distribute_dofs(fe);
  
  Vector<double> shift(shift_dh.n_dofs());

  shift.add(+1);

  GridOut grid_out;
  grid_out.set_flags (GridOutFlags::Ucd(true));
  grid_out.write_ucd (tria, logfile);

  QTrapez<dim> quad;
  MappingQ1Eulerian<dim,Vector<double>,spacedim> mapping(shift, shift_dh);
  
  typename Triangulation<dim,spacedim>::active_cell_iterator cell=tria.begin_active(), 
    endc=tria.end() ;
  Point<spacedim> real;
  Point<dim> unit;
  for(;cell!=endc;++cell)
    {
      deallog<<cell<< std::endl;	
      for(unsigned int q=0; q<quad.size(); ++q)
	{
	  real = mapping.transform_unit_to_real_cell(cell, quad.point(q));
	  // unit = mapping.transform_real_to_unit_cell(cell, real);
	  deallog<<quad.point(q)<< " -> " << real << std::endl;
	}
    }	   
    
}

int main () 
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  test<1,2>("grids/circle_1.inp");
  test<2,3>("grids/square.inp");
  test<2,3>("grids/sphere_1.inp");

  return 0;
}
                  
