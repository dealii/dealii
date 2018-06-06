// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// Test that all the iterator types are defined.

#include <deal.II/base/index_set.h>
#include <deal.II/base/iterator_range.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_vector_base.h>

#include <deal.II/particles/particle_iterator.h>

#include "../tests.h"

int
main()
{
  initlog();

  typedef std::iterator_traits<TriaRawIterator<TriaAccessor<2, 2, 2>>>::
    iterator_category tria_raw_iterator_category;
  typedef std::iterator_traits<TriaIterator<TriaAccessor<2, 2, 2>>>::
    iterator_category tria_iterator_category;
  typedef std::iterator_traits<TriaActiveIterator<TriaAccessor<2, 2, 2>>>::
    iterator_category tria_active_iterator_category;

  typedef std::iterator_traits<
    Particles::ParticleIterator<2>>::iterator_category particle_category;

  typedef std::iterator_traits<IteratorRange<
    Particles::ParticleIterator<2>>::IteratorOverIterators>::iterator_category
    iterator_over_iterator_category;

  typedef std::iterator_traits<
    internal::BlockVectorIterators::Iterator<BlockVector<double>, false>>::
    iterator_category block_vector_base_iterator_category;

  typedef std::iterator_traits<IndexSet::IntervalIterator>::iterator_category
    intervall_iterator_category;
  typedef std::iterator_traits<IndexSet::ElementIterator>::iterator_category
    element_iterator_category;

  deallog << "OK" << std::endl;
}
