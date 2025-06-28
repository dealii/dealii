// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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

  using tria_raw_iterator_category = std::iterator_traits<
    TriaRawIterator<TriaAccessor<2, 2, 2>>>::iterator_category;
  using tria_iterator_category = std::iterator_traits<
    TriaIterator<TriaAccessor<2, 2, 2>>>::iterator_category;
  using tria_active_iterator_category = std::iterator_traits<
    TriaActiveIterator<TriaAccessor<2, 2, 2>>>::iterator_category;

  using particle_category =
    std::iterator_traits<Particles::ParticleIterator<2>>::iterator_category;

  using iterator_over_iterator_category = std::iterator_traits<IteratorRange<
    Particles::ParticleIterator<2>>::IteratorOverIterators>::iterator_category;

  using block_vector_base_iterator_category = std::iterator_traits<
    internal::BlockVectorIterators::Iterator<BlockVector<double>,
                                             false>>::iterator_category;

  using interval_iterator_category =
    std::iterator_traits<IndexSet::IntervalIterator>::iterator_category;
  using element_iterator_category =
    std::iterator_traits<IndexSet::ElementIterator>::iterator_category;

  deallog << "OK" << std::endl;
}
