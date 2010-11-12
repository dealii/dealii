//----------------------------  vertex_as_face_07.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vertex_as_face_07.cc  ---------------------------

// verify that we can do things like cell->face() in 1d as well. here:
// test cell->face(0)->get_dof_indices()
// compared to _06, we now test for an hp DoFHandler

#include "../tests.h"
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>

#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <hp/dof_handler.h>
#include <hp/fe_collection.h>
#include <dofs/dof_accessor.h>

#include <fstream>


template <int spacedim>
void print_dofs (const typename hp::DoFHandler<1,spacedim>::face_iterator &i,
		 const unsigned int fe_index,
		 const unsigned int n)
{
  std::vector<unsigned int> dof_indices (n);
  i->get_dof_indices (dof_indices, fe_index);
  for (unsigned int i=0; i<n; ++i)
    deallog << dof_indices[i] << ' ';
  deallog << std::endl;
}



template <int spacedim>
void print_dofs (const typename hp::DoFHandler<1,spacedim>::cell_iterator &i,
		 const unsigned int n)
{
  std::vector<unsigned int> dof_indices (n);
  i->get_dof_indices (dof_indices);
  for (unsigned int i=0; i<n; ++i)
    deallog << dof_indices[i] << ' ';
  deallog << std::endl;
}



template <int spacedim>
void test ()
{
  Triangulation<1,spacedim> tria;
  GridGenerator::hyper_cube (tria);

  FESystem<1,spacedim> fe1(FE_Q<1,spacedim>(2),1,
			   FE_Q<1,spacedim>(1),1);
  FESystem<1,spacedim> fe2(FE_Q<1,spacedim>(3),1,
			   FE_Q<1,spacedim>(2),1);
  hp::FECollection<1,spacedim> fe_collection;
  fe_collection.push_back (fe1);
  fe_collection.push_back (fe2);

  hp::DoFHandler<1,spacedim> dof_handler (tria);
  dof_handler.begin_active()->set_active_fe_index (0);
  dof_handler.distribute_dofs (fe_collection);

  deallog << "Coarse mesh:" << std::endl;
  print_dofs<spacedim> (dof_handler.begin_active()->face(0), 0, fe1.dofs_per_face);
  print_dofs<spacedim> (dof_handler.begin_active()->face(1), 0, fe1.dofs_per_face);

  tria.refine_global (2);
  {
    unsigned int index = 0;
    for (typename hp::DoFHandler<1,spacedim>::active_cell_iterator
	   cell = dof_handler.begin_active();
	 cell != dof_handler.end(); ++cell, index = (index + 1) % fe_collection.size())
      cell->set_active_fe_index (index);
  }
  dof_handler.distribute_dofs (fe_collection);

  for (typename hp::DoFHandler<1,spacedim>::active_cell_iterator
	 cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      deallog << "Cell: " << cell
	      << ", active_fe_index=" << cell->active_fe_index()
	      << std::endl;

      print_dofs<spacedim> (cell,
			    cell->get_fe().dofs_per_cell);
      print_dofs<spacedim> (cell->face(0),
			    cell->active_fe_index(),
			    cell->get_fe().dofs_per_face);
      print_dofs<spacedim> (cell->face(1),
			    cell->active_fe_index(),
			    cell->get_fe().dofs_per_face);
    }
}



int main ()
{
  std::ofstream logfile("vertex_as_face_07/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();

  return 0;
}
