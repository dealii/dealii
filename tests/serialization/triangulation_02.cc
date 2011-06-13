//----------------------------  triangulation_02.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  triangulation_02.cc  ---------------------------

// check serialization for Triangulation<1,dim>. do the same as in the
// _01 test but for an adaptively refined and coarsened triangulation

#include "serialization.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include "../../include/deal.II/grid/tria.h"


template <int dim, int spacedim>
bool operator == (const Triangulation<dim,spacedim> &t1,
		  const Triangulation<dim,spacedim> &t2)
{
				   // test a few attributes, though we can't
				   // test everything unfortunately...
  if (t1.n_active_cells() != t2.n_active_cells())
    return false;
  
  if (t1.n_cells() != t2.n_cells())
    return false;

  if (t1.n_faces() != t2.n_faces())
    return false;

  typename Triangulation<dim,spacedim>::cell_iterator
    c1 = t1.begin(),
    c2 = t2.begin();
  for (; (c1 != t1.end()) && (c2 != t2.end()); ++c1, ++c2)
    {
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
	{
	  if (c1->vertex(v) != c2->vertex(v))
	    return false;
	  if (c1->vertex_index(v) != c2->vertex_index(v))
	    return false;
	}

      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	{
	  if (c1->face(f)->at_boundary() != c2->face(f)->at_boundary())
	    return false;

	  if (c1->face(f)->at_boundary())
	    {
	      if (c1->face(f)->boundary_indicator() !=
		  c2->face(f)->boundary_indicator())
		return false;
	    }
	  else
	    {
	      if (c1->neighbor(f)->level() != c2->neighbor(f)->level())
		return false;
	      if (c1->neighbor(f)->index() != c2->neighbor(f)->index())
		return false;
	    }    
	}
      
      if (c1->subdomain_id() != c2->subdomain_id())
	return false;
      
      if (c1->material_id() != c2->material_id())
	return false;

      if (c1->user_index() != c2->user_index())
	return false;

      if (c1->user_flag_set() != c2->user_flag_set())
	return false;
    }

  // also check the order of raw iterators as they contain
  // something about the history of the triangulation
  typename Triangulation<dim,spacedim>::cell_iterator
    r1 = t1.begin_raw(),
    r2 = t2.begin_raw();
  for (; (r1 != t1.end()) && (r2 != t2.end()); ++r1, ++r2)
  {
    if (r1->level() != r2->level())
      return false;
    if (r1->index() != r2->index())
      return false;
  }

  return true;
}

template <int dim, int spacedim>
void do_boundary (Triangulation<dim,spacedim> &t1)
{
  typename Triangulation<dim,spacedim>::cell_iterator
    c1 = t1.begin();
  for (; c1 != t1.end(); ++c1)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (c1->at_boundary(f))
	c1->face(f)->set_boundary_indicator (42);
}


template <int spacedim>
void do_boundary (Triangulation<1,spacedim> &)
{}


template <int dim, int spacedim>
void test ()
{
  Triangulation<dim,spacedim> tria_1, tria_2;

  GridGenerator::hyper_cube (tria_1);
  tria_1.refine_global (2);
  // coarsen again as this takes away the finest level but may leave
  // around some of this level's cells
  for (typename Triangulation<dim,spacedim>::active_cell_iterator
    cell = tria_1.begin_active(2); cell != tria_1.end(); ++cell)
    cell->set_coarsen_flag ();
  tria_1.execute_coarsening_and_refinement();
  // now add one cell again
  tria_1.begin_active()->set_refine_flag ();
  tria_1.execute_coarsening_and_refinement();
  
  tria_1.begin_active()->set_subdomain_id (1);
  tria_1.begin_active()->set_material_id (2);
  tria_1.begin_active()->set_user_index (3);
  tria_1.begin_active()->set_user_flag ();
  tria_1.begin_active()->set_refine_flag (RefinementCase<dim>::cut_x);

  do_boundary (tria_1);
  
  verify (tria_1, tria_2);
}


int main ()
{
  std::ofstream logfile("triangulation_02/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1,1> ();
  test<1,2> ();
  test<2,2> ();
  test<2,3> ();
  test<3,3> ();
    
  deallog << "OK" << std::endl;
}
