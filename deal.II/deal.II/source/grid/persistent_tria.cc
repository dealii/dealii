/* $Id$ */

#include <grid/persistent_tria.h>


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
PersistentTriangulation<dim>::block_write (ostream &) const 
{
  Assert (false, ExcFunctionNotUseful());
};



template <int dim>
void
PersistentTriangulation<dim>::block_read (istream &)
{
  Assert (false, ExcFunctionNotUseful());
};



template <int dim>
void
PersistentTriangulation<dim>::create_triangulation (const vector<Point<dim> >    &,
						    const vector<CellData<dim> > &,
						    const SubCellData            &)
{
  Assert (false, ExcFunctionNotUseful());
};




// explicit instantiations
template class PersistentTriangulation<deal_II_dimension>;

