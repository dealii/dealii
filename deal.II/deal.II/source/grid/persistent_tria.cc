//----------------------------  persistent_tria.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  persistent_tria.cc  ---------------------------


#include <grid/persistent_tria.h>
#include <grid/magic_numbers.h>
#include <iostream>


template <int dim>
PersistentTriangulation<dim>::
PersistentTriangulation (const Triangulation<dim> &coarse_grid) :
		coarse_grid (&coarse_grid) 
{};


template <int dim>
PersistentTriangulation<dim>::
PersistentTriangulation (const PersistentTriangulation<dim> &old_tria) :
						 // default initialize
						 // tria, i.e. it will be
						 // empty on first use
		Triangulation<dim> (),
		coarse_grid (old_tria.coarse_grid),
		refine_flags (old_tria.refine_flags),
		coarsen_flags (old_tria.coarsen_flags)
{
  Assert (old_tria.n_levels() == 0, ExcTriaNotEmpty ());
};


template <int dim>
PersistentTriangulation<dim>::~PersistentTriangulation () 
{};


template <int dim>
void
PersistentTriangulation<dim>::execute_coarsening_and_refinement () 
{
				   // first save flags
  refine_flags.push_back (vector<bool>());
  coarsen_flags.push_back (vector<bool>());
  save_refine_flags (refine_flags.back());
  save_coarsen_flags (coarsen_flags.back());

				   // then refine triangulation
  Triangulation<dim>::execute_coarsening_and_refinement ();
};


template <int dim>
void
PersistentTriangulation<dim>::restore () {
				   // copy the old triangulation.
				   // this will yield an error if the
				   // underlying triangulation was not
				   // empty
  Triangulation<dim>::copy_triangulation (*coarse_grid);

				   // for each of the previous refinement
				   // sweeps
  for (unsigned int i=0; i<refine_flags.size(); ++i) 
    {
				       // get flags
      load_refine_flags  (refine_flags[i]);
      load_coarsen_flags (coarsen_flags[i]);

      Triangulation<dim>::execute_coarsening_and_refinement ();
    };
};


template <int dim>
void
PersistentTriangulation<dim>::copy_triangulation (const Triangulation<dim> &old_grid) 
{
  clear ();
  coarse_grid  = &old_grid;
  refine_flags.clear ();
  coarsen_flags.clear ();
};


template <int dim>
void
PersistentTriangulation<dim>::create_triangulation (const vector<Point<dim> >    &,
						    const vector<CellData<dim> > &,
						    const SubCellData            &)
{
  Assert (false, ExcFunctionNotUseful());
};


template <int dim>
void
PersistentTriangulation<dim>::write_flags(ostream &out) const
{
  const unsigned int n_flag_levels=refine_flags.size();
  
  AssertThrow (out, ExcIO());

  out << mn_persistent_tria_flags_begin << ' ' << n_flag_levels << endl;

  for (unsigned int i=0; i<n_flag_levels; ++i)
    {
      write_bool_vector (mn_tria_refine_flags_begin, refine_flags[i],
			 mn_tria_refine_flags_end, out);
      write_bool_vector (mn_tria_coarsen_flags_begin, coarsen_flags[i],
			 mn_tria_coarsen_flags_end, out);
    }
  
  out << mn_persistent_tria_flags_end << endl;

  AssertThrow (out, ExcIO());
}


template <int dim>
void
PersistentTriangulation<dim>::read_flags(istream &in)
{
  Assert(refine_flags.size()==0 && coarsen_flags.size()==0, ExcTriaNotEmpty());

  AssertThrow (in, ExcIO());

  unsigned int magic_number;
  in >> magic_number;
  AssertThrow(magic_number==mn_persistent_tria_flags_begin, ExcGridReadError());

  unsigned int n_flag_levels;
  in >> n_flag_levels;
  for (unsigned int i=0; i<n_flag_levels; ++i)
    {
      refine_flags.push_back (vector<bool>());
      coarsen_flags.push_back (vector<bool>());
      read_bool_vector (mn_tria_refine_flags_begin, refine_flags.back(),
			mn_tria_refine_flags_end, in);
      read_bool_vector (mn_tria_coarsen_flags_begin, coarsen_flags.back(),
			mn_tria_coarsen_flags_end, in);
    }
  
  in >> magic_number;
  AssertThrow(magic_number==mn_persistent_tria_flags_end, ExcGridReadError());

  AssertThrow (in, ExcIO());
}


// explicit instantiations
template class PersistentTriangulation<deal_II_dimension>;

