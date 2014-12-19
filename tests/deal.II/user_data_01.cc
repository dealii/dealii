// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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
check_user_pointers(Triangulation<dim> &tr)
{
  // Check if values are the same as below
  Triangulation<dim> *p = &tr;
  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end(); ++it)
    {
      deallog << '.';
      if (it->user_pointer() != p++)
        deallog << "Error" << std::endl;
    }
  deallog << "cells" << std::endl;

  // now for faces and lines. since
  // we visit them multiple times,
  // check only the first time around
  tr.clear_user_flags ();
  if (dim > 1)
    for (typename Triangulation<dim>::cell_iterator it = tr.begin();
         it != tr.end(); ++it)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (it->face(f)->user_flag_set() == false)
          {
            deallog << '.';
            if (it->face(f)->user_pointer() != p++)
              deallog << "Error" << std::endl;
            it->face(f)->set_user_flag();
          }
  deallog << "faces" << std::endl;

  if (dim > 2)
    for (typename Triangulation<dim>::cell_iterator it = tr.begin();
         it != tr.end(); ++it)
      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
        if (it->line(l)->user_flag_set() == false)
          {
            deallog << '.';
            if (it->line(l)->user_pointer() != p++)
              deallog << "Error" << std::endl;
            it->line(l)->set_user_flag();
          }
  deallog << "lines" << std::endl;
}

/**
 * Check if user indices are the same as entered below
 */
template<int dim>
void
check_user_indices(Triangulation<dim> &tr)
{
  // Check if values are the same as below
  unsigned int p=1;
  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end(); ++it)
    {
      deallog << '.';
      if (it->user_index() != p++)
        deallog << "Error" << std::endl;
    }
  deallog << "cells" << std::endl;

  // now for faces and lines. since
  // we visit them multiple times,
  // check only the first time around
  tr.clear_user_flags ();
  if (dim > 1)
    for (typename Triangulation<dim>::cell_iterator it = tr.begin();
         it != tr.end(); ++it)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (it->face(f)->user_flag_set() == false)
          {
            deallog << '.';
            if (it->face(f)->user_index() != p++)
              deallog << "Error" << std::endl;
            it->face(f)->set_user_flag();
          }
  deallog << "faces" << std::endl;

  if (dim > 2)
    for (typename Triangulation<dim>::cell_iterator it = tr.begin();
         it != tr.end(); ++it)
      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
        if (it->line(l)->user_flag_set() == false)
          {
            deallog << '.';
            if (it->line(l)->user_index() != p++)
              deallog << "Error" << std::endl;
            it->line(l)->set_user_flag();
          }
  deallog << "lines" << std::endl;
}


template<int dim>
void
user_pointers(Triangulation<dim> &tr)
{
  deallog << "Pointers" << dim << 'D' << std::endl;

  // Fill user pointers with some
  // nonsense. clear them first
  Triangulation<dim> *p = &tr;
  tr.clear_user_pointers();
  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end(); ++it)
    it->set_user_pointer(p++);

  // we hit faces and lines more than
  // once, possibly. only set them
  // the first time around
  tr.clear_user_flags ();
  if (dim > 1)
    for (typename Triangulation<dim>::cell_iterator it = tr.begin();
         it != tr.end(); ++it)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (it->face(f)->user_flag_set() == false)
          {
            it->face(f)->set_user_pointer(p++);
            it->face(f)->set_user_flag();
          }

  if (dim > 2)
    for (typename Triangulation<dim>::cell_iterator it = tr.begin();
         it != tr.end(); ++it)
      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
        if (it->line(l)->user_flag_set() == false)
          {
            it->line(l)->set_user_pointer(p++);
            it->line(l)->set_user_flag();
          }

  // Check if they are still the same
  check_user_pointers(tr);
  // Create two pointer index clashes here
  try
    {
      tr.begin()->user_index();
    }
  catch (...)
    {
    }
  tr.begin()->user_pointer();


  // Check if save and load work
  std::vector<void *> cell_pointers(tr.n_cells());
  deallog << "Save" << dim << 'D' << std::endl;
  tr.save_user_pointers(cell_pointers);
  tr.clear_user_data();
  deallog << "Load" << dim << 'D' << std::endl;
  tr.load_user_pointers(cell_pointers);
  check_user_pointers(tr);
}


template<int dim>
void
user_indices(Triangulation<dim> &tr)
{
  deallog << "Indices" << dim << 'D' << std::endl;

  // Fill user indices with some
  // nonsense. clear them first

  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end(); ++it)
    {
      it->clear_user_index();

      if (dim > 1)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          it->face(f)->clear_user_index();

      if (dim > 2)
        for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
          it->line(l)->clear_user_index();
    }


  unsigned int p=1;
  for (typename Triangulation<dim>::cell_iterator it = tr.begin();
       it != tr.end(); ++it)
    it->set_user_index(p++);

  // we hit faces and lines more than
  // once, possibly. only set them
  // the first time around
  tr.clear_user_flags ();
  if (dim > 1)
    for (typename Triangulation<dim>::cell_iterator it = tr.begin();
         it != tr.end(); ++it)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (it->face(f)->user_flag_set() == false)
          {
            it->face(f)->set_user_index(p++);
            it->face(f)->set_user_flag();
          }

  if (dim > 2)
    for (typename Triangulation<dim>::cell_iterator it = tr.begin();
         it != tr.end(); ++it)
      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
        if (it->line(l)->user_flag_set() == false)
          {
            it->line(l)->set_user_index(p++);
            it->line(l)->set_user_flag();
          }

  // Check if they are still the same
  check_user_indices(tr);
  // Create two pointer index clashes here
  try
    {
      tr.begin()->user_pointer();
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
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

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2>();
  check<3>();
}
