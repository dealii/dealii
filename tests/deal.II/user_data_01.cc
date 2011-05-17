//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// Check the handling of user pointers and user indices

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <string>




/**
 * Check if user pointers are the same as entered below
 */
template<int dim>
void
check_user_pointers(Triangulation<dim>& tr)
{
				   // Check if values are the same as below
  Triangulation<dim>* p = &tr;
  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end();++it)
    {
      deallog << '.';
      if (it->user_pointer() != p++)
	deallog << "Error" << std::endl;
    }
  deallog << "cells" << std::endl;

  if (dim > 1)
    for (typename Triangulation<dim>::face_iterator it = tr.begin_face();
	 it != tr.end_face();++it)
    {
      deallog << '.';
      if (it->user_pointer() != p++)
	deallog << "Error" << std::endl;
    }
  deallog << "faces" << std::endl;

  if (dim > 2)
    for (typename Triangulation<dim>::line_iterator it = tr.begin_line();
	 it != tr.end_line();++it)
    {
      deallog << '.';
      if (it->user_pointer() != p++)
	deallog << "Error" << std::endl;
    }
  deallog << "lines" << std::endl;
}

/**
 * Check if user indices are the same as entered below
 */
template<int dim>
void
check_user_indices(Triangulation<dim>& tr)
{
				   // Check if values are the same as below
  unsigned int p=1;
  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end();++it)
    {
      deallog << '.';
      if (it->user_index() != p++)
	deallog << "Error" << std::endl;
    }
  deallog << "cells" << std::endl;

  if (dim > 1)
    for (typename Triangulation<dim>::face_iterator it = tr.begin_face();
	 it != tr.end_face();++it)
    {
      deallog << '.';
      if (it->user_index() != p++)
	deallog << "Error" << std::endl;
    }
  deallog << "faces" << std::endl;

  if (dim > 2)
    for (typename Triangulation<dim>::line_iterator it = tr.begin_line();
	 it != tr.end_line();++it)
    {
      deallog << '.';
      if (it->user_index() != p++)
	deallog << "Error" << std::endl;
    }
  deallog << "lines" << std::endl;
}


template<int dim>
void
user_pointers(Triangulation<dim>& tr)
{
  deallog << "Pointers" << dim << 'D' << std::endl;
  
				   // Fill user pointers with some nonsense
  Triangulation<dim>* p = &tr;
  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end();++it)
    it->set_user_pointer(p++);
  
  if (dim > 1)
    for (typename Triangulation<dim>::face_iterator it = tr.begin_face();
	 it != tr.end_face();++it)
      it->set_user_pointer(p++);
  
  if (dim > 2)
    for (typename Triangulation<dim>::line_iterator it = tr.begin_line();
	 it != tr.end_line();++it)
      it->set_user_pointer(p++);

				   // Check if they are still the same
  check_user_pointers(tr);
				   // Create two pointer index clashes here
  tr.begin()->user_index();
  tr.begin()->user_pointer();
  
  
				   // Check if save and load work
  std::vector<void*> cell_pointers(tr.n_cells());
  deallog << "Save" << dim << 'D' << std::endl;
  tr.save_user_pointers(cell_pointers);
  tr.clear_user_data();
  deallog << "Load" << dim << 'D' << std::endl;
  tr.load_user_pointers(cell_pointers);
  check_user_pointers(tr);  
}


template<int dim>
void
user_indices(Triangulation<dim>& tr)
{
  deallog << "Indices" << dim << 'D' << std::endl;
  
				   // Fill user pointers with some nonsense
  unsigned int p=1;
  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end();++it)
    it->set_user_index(p++);
  
  if (dim > 1)
    for (typename Triangulation<dim>::face_iterator it = tr.begin_face();
	 it != tr.end_face();++it)
      it->set_user_index(p++);
  
  if (dim > 2)
    for (typename Triangulation<dim>::line_iterator it = tr.begin_line();
	 it != tr.end_line();++it)
      it->set_user_index(p++);

				   // Check if they are still the same
  check_user_indices(tr);
				   // Create two pointer index clashes here
  tr.begin()->user_pointer();
  tr.begin()->user_index();
  
  
				   // Check if save and load work
  std::vector<unsigned int> indices(tr.n_cells());
  deallog << "Save" << dim << 'D' << std::endl;
  tr.save_user_indices(indices);
  tr.clear_user_data();
  deallog << "Load" << dim << 'D' << std::endl;
  tr.load_user_indices(indices);
  check_user_indices(tr);  
}


template <int dim>
void check()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  
  user_pointers(tr);
  tr.clear_user_pointers();
  user_indices(tr);
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();
  std::ofstream logfile("user_data_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2>();
  check<3>();  
}
