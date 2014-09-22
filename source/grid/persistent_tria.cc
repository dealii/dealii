// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/grid/persistent_tria.h>
#include <deal.II/grid/magic_numbers.h>
#include <iostream>

DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
const unsigned int
PersistentTriangulation<dim,spacedim>::dimension;

template <int dim, int spacedim>
const unsigned int
PersistentTriangulation<dim,spacedim>::spacedimension;


template <int dim, int spacedim>
PersistentTriangulation<dim,spacedim>::
PersistentTriangulation (const Triangulation<dim,spacedim> &coarse_grid)
  :
  coarse_grid (&coarse_grid, typeid(*this).name())
{}



template <int dim, int spacedim>
PersistentTriangulation<dim,spacedim>::
PersistentTriangulation (const PersistentTriangulation<dim,spacedim> &old_tria)
  :
  // default initialize
  // tria, i.e. it will be
  // empty on first use
  Triangulation<dim,spacedim> (),
  coarse_grid (old_tria.coarse_grid),
  refine_flags (old_tria.refine_flags),
  coarsen_flags (old_tria.coarsen_flags)
{
  Assert (old_tria.n_levels() == 0, ExcTriaNotEmpty ());
}



template <int dim, int spacedim>
PersistentTriangulation<dim,spacedim>::~PersistentTriangulation ()
{}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::execute_coarsening_and_refinement ()
{
  // first save flags
  refine_flags.push_back (std::vector<bool>());
  coarsen_flags.push_back (std::vector<bool>());
  this->save_refine_flags (refine_flags.back());
  this->save_coarsen_flags (coarsen_flags.back());

  // then refine triangulation. if
  // this function throws an
  // exception, that's fine since it
  // is the last call here
  Triangulation<dim,spacedim>::execute_coarsening_and_refinement ();
}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::restore ()
{
  // for each of the previous
  // refinement sweeps
  for (unsigned int i=0; i<refine_flags.size()+1; ++i)
    restore(i);
}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::restore (const unsigned int step)
{

  if (step==0)
    // copy the old triangulation.
    // this will yield an error if
    // the underlying triangulation
    // was not empty
    Triangulation<dim,spacedim>::copy_triangulation (*coarse_grid);
  else
    // for each of the previous
    // refinement sweeps
    {
      Assert(step<refine_flags.size()+1,
             ExcDimensionMismatch(step, refine_flags.size()+1));

      this->load_refine_flags  (refine_flags[step-1]);
      this->load_coarsen_flags (coarsen_flags[step-1]);

      Triangulation<dim,spacedim>::execute_coarsening_and_refinement ();
    }
}



template <int dim, int spacedim>
unsigned int
PersistentTriangulation<dim,spacedim>::n_refinement_steps() const
{
  return refine_flags.size();
}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::copy_triangulation (const Triangulation<dim,spacedim> &old_grid)
{
  this->clear ();
  coarse_grid  = &old_grid;
  refine_flags.clear ();
  coarsen_flags.clear ();
}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::create_triangulation (const std::vector<Point<spacedim> > &,
    const std::vector<CellData<dim> > &,
    const SubCellData &)
{
  Assert (false, ExcImpossibleInDim(dim));
}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::create_triangulation_compatibility (
  const std::vector<Point<spacedim> > &,
  const std::vector<CellData<dim> > &,
  const SubCellData &)
{
  Assert (false, ExcImpossibleInDim(dim));
}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::write_flags(std::ostream &out) const
{
  const unsigned int n_flag_levels=refine_flags.size();

  AssertThrow (out, ExcIO());

  out << mn_persistent_tria_flags_begin << ' ' << n_flag_levels << std::endl;

  for (unsigned int i=0; i<n_flag_levels; ++i)
    {
      this->write_bool_vector (mn_tria_refine_flags_begin, refine_flags[i],
                               mn_tria_refine_flags_end, out);
      this->write_bool_vector (mn_tria_coarsen_flags_begin, coarsen_flags[i],
                               mn_tria_coarsen_flags_end, out);
    }

  out << mn_persistent_tria_flags_end << std::endl;

  AssertThrow (out, ExcIO());
}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::read_flags(std::istream &in)
{
  Assert(refine_flags.size()==0 && coarsen_flags.size()==0,
         ExcFlagsNotCleared());
  AssertThrow (in, ExcIO());

  unsigned int magic_number;
  in >> magic_number;
  AssertThrow(magic_number==mn_persistent_tria_flags_begin,
              typename Triangulation<dim>::ExcGridReadError());

  unsigned int n_flag_levels;
  in >> n_flag_levels;
  for (unsigned int i=0; i<n_flag_levels; ++i)
    {
      refine_flags.push_back (std::vector<bool>());
      coarsen_flags.push_back (std::vector<bool>());
      this->read_bool_vector (mn_tria_refine_flags_begin, refine_flags.back(),
                              mn_tria_refine_flags_end, in);
      this->read_bool_vector (mn_tria_coarsen_flags_begin, coarsen_flags.back(),
                              mn_tria_coarsen_flags_end, in);
    }

  in >> magic_number;
  AssertThrow(magic_number==mn_persistent_tria_flags_end,
              typename Triangulation<dim>::ExcGridReadError());

  AssertThrow (in, ExcIO());
}



template <int dim, int spacedim>
void
PersistentTriangulation<dim,spacedim>::clear_flags()
{
  refine_flags.clear();
  coarsen_flags.clear();
}



template <int dim, int spacedim>
std::size_t
PersistentTriangulation<dim,spacedim>::memory_consumption () const
{
  return (Triangulation<dim,spacedim>::memory_consumption () +
          MemoryConsumption::memory_consumption (coarse_grid) +
          MemoryConsumption::memory_consumption (refine_flags) +
          MemoryConsumption::memory_consumption (coarsen_flags));
}


// explicit instantiations
template class PersistentTriangulation<1>;
template class PersistentTriangulation<2>;
template class PersistentTriangulation<3>;
template class PersistentTriangulation<1,2>;
template class PersistentTriangulation<2,3>;

DEAL_II_NAMESPACE_CLOSE

