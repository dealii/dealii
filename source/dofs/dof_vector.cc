// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
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

#include <deal.II/dofs/dof_vector.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  template <class DH, class VECTOR>
  void reinit_dof_vector (const DH &, VECTOR &)
  {
    Assert(false, ExcNotImplemented());
  }

  template <int dim, int spacedim, typename number>
  void reinit_dof_vector (const DoFHandler<dim, spacedim> &dh,
                          Vector<number> &v)
  {
    v.reinit(dh.n_dofs());
  }

  template <int dim, int spacedim, typename number>
  void reinit_dof_vector (const DoFHandler<dim, spacedim> &dh,
                          MGLevelObject<Vector<number> > &v)
  {
    const unsigned int nlev = dh.get_tria().n_levels();
    v.resize(0, nlev-1);
    for (unsigned int l=0; l<nlev; ++l)
      v[l].reinit(dh.n_dofs(l));
  }

  template <int dim, int spacedim, typename number>
  void reinit_dof_vector (const DoFHandler<dim, spacedim> &dh,
                          BlockVector<number> &v)
  {
    v.reinit(dh.block_info().global());
  }

  template <int dim, int spacedim, typename number>
  void reinit_dof_vector (const DoFHandler<dim, spacedim> &dh,
                          MGLevelObject<BlockVector<number> > &v)
  {
    const unsigned int nlev = dh.get_tria().n_levels();
    v.resize(0, nlev-1);
    for (unsigned int l=0; l<nlev; ++l)
      v[l].reinit(dh.block_info().level(l));
  }

  template <int dim, int spacedim, typename number>
  void reinit_dof_vector (const hp::DoFHandler<dim, spacedim> &dh,
                          Vector<number> &v)
  {
    v.reinit(dh.n_dofs());
  }

}

//----------------------------------------------------------------------//

template <class DH, class VECTOR>
const ConstraintMatrix DoFVector<DH, VECTOR>::no_constraints;


template <class DH, class VECTOR>
DoFVector<DH, VECTOR>::DoFVector (const DH &dh, const ConstraintMatrix &co)
  : dh(&dh),
    constraints(&co),
    my_data(0),
    other_data(0)
{
  my_data = mem.alloc();
  const_data = my_data;
}


template <class DH, class VECTOR>
DoFVector<DH, VECTOR>::DoFVector (const DH &dh)
  : dh(&dh),
    constraints(&no_constraints),
    my_data(0),
    other_data(0)
{
  my_data = mem.alloc();
  const_data = my_data;
}


template <class DH, class VECTOR>
DoFVector<DH, VECTOR>::DoFVector (const DH &dh, const ConstraintMatrix &co, VECTOR &v)
  : dh(&dh),
    constraints(&co),
    my_data(0),
    other_data(&v),
    const_data(&v)
{}


template <class DH, class VECTOR>
DoFVector<DH, VECTOR>::DoFVector (const DH &dh, const ConstraintMatrix &co, const VECTOR &v)
  : dh(&dh),
    constraints(&co),
    my_data(0),
    other_data(0),
    const_data(&v)
{}


template <class DH, class VECTOR>
DoFVector<DH, VECTOR>::DoFVector (const DH &dh, VECTOR &v)
  : dh(&dh),
    constraints(&no_constraints),
    my_data(0),
    other_data(&v),
    const_data(&v)
{}


template <class DH, class VECTOR>
DoFVector<DH, VECTOR>::DoFVector (const DH &dh, const VECTOR &v)
  : dh(&dh),
    constraints(&no_constraints),
    my_data(0),
    other_data(0),
    const_data(&v)
{}


template <class DH, class VECTOR>
DoFVector<DH, VECTOR> &
DoFVector<DH, VECTOR>::operator= (const DoFVector<DH, VECTOR> &other)
{
  Assert(other_data == 0, ExcNotOwner());

  dh = &other.get_dof_handler();
  constraints = &other.get_constraint_matrix();
  *my_data = other.get_data();
  return *this;
}


template <class DH, class VECTOR>
DoFVector<DH, VECTOR>::~DoFVector ()
{
  if (my_data != 0)
    mem.free(my_data);
}


template <class DH, class VECTOR>
void
DoFVector<DH, VECTOR>::sync ()
{
  Assert(other_data == 0, ExcNotOwner());
  reinit_dof_vector(get_dof_handler(), get_data());
}


#include "dof_vector.inst"


DEAL_II_NAMESPACE_CLOSE
