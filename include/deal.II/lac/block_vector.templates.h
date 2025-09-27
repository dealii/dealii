// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_block_vector_templates_h
#define dealii_block_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <algorithm>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <typename Number>
BlockVector<Number>::BlockVector(const unsigned int n_blocks,
                                 const size_type    block_size)
{
  reinit(n_blocks, block_size);
}



template <typename Number>
BlockVector<Number>::BlockVector(const std::vector<size_type> &block_sizes)
{
  reinit(block_sizes, false);
}


template <typename Number>
BlockVector<Number>::BlockVector(const BlockIndices &n)
{
  reinit(n, false);
}


template <typename Number>
BlockVector<Number>::BlockVector(const BlockVector<Number> &v)
  : BlockVectorBase<Vector<Number>>()
{
  this->components.resize(v.n_blocks());
  this->block_indices = v.block_indices;

  for (size_type i = 0; i < this->n_blocks(); ++i)
    this->components[i] = v.components[i];
}



template <typename Number>
template <typename OtherNumber>
BlockVector<Number>::BlockVector(const BlockVector<OtherNumber> &v)
{
  reinit(v, true);
  *this = v;
}



#ifdef DEAL_II_WITH_TRILINOS

template <typename Number>
BlockVector<Number>::BlockVector(const TrilinosWrappers::MPI::BlockVector &v)
{
  this->block_indices = v.get_block_indices();
  this->components.resize(this->n_blocks());

  for (size_type i = 0; i < this->n_blocks(); ++i)
    this->components[i] = v.block(i);

  BaseClass::collect_sizes();
}

#endif


template <typename Number>
void
BlockVector<Number>::reinit(const unsigned int n_blocks,
                            const size_type    block_size,
                            const bool         omit_zeroing_entries)
{
  std::vector<size_type> block_sizes(n_blocks, block_size);
  reinit(block_sizes, omit_zeroing_entries);
}


template <typename Number>
void
BlockVector<Number>::reinit(const std::vector<size_type> &block_sizes,
                            const bool                    omit_zeroing_entries)
{
  this->block_indices.reinit(block_sizes);
  if (this->components.size() != this->n_blocks())
    this->components.resize(this->n_blocks());

  for (size_type i = 0; i < this->n_blocks(); ++i)
    this->components[i].reinit(block_sizes[i], omit_zeroing_entries);
}


template <typename Number>
void
BlockVector<Number>::reinit(const BlockIndices &n,
                            const bool          omit_zeroing_entries)
{
  this->block_indices = n;
  if (this->components.size() != this->n_blocks())
    this->components.resize(this->n_blocks());

  for (size_type i = 0; i < this->n_blocks(); ++i)
    this->components[i].reinit(n.block_size(i), omit_zeroing_entries);
}


template <typename Number>
template <typename Number2>
void
BlockVector<Number>::reinit(const BlockVector<Number2> &v,
                            const bool                  omit_zeroing_entries)
{
  this->block_indices = v.get_block_indices();
  if (this->components.size() != this->n_blocks())
    this->components.resize(this->n_blocks());

  for (size_type i = 0; i < this->n_blocks(); ++i)
    this->block(i).reinit(v.block(i), omit_zeroing_entries);
}



#ifdef DEAL_II_WITH_TRILINOS
template <typename Number>
inline BlockVector<Number> &
BlockVector<Number>::operator=(const TrilinosWrappers::MPI::BlockVector &v)
{
  BaseClass::operator=(v);
  return *this;
}
#endif


template <typename Number>
void
BlockVector<Number>::swap(BlockVector<Number> &v) noexcept
{
  std::swap(this->components, v.components);

  dealii::swap(this->block_indices, v.block_indices);
}



template <typename Number>
void
BlockVector<Number>::print(std::ostream      &out,
                           const unsigned int precision,
                           const bool         scientific,
                           const bool         across) const
{
  for (size_type i = 0; i < this->n_blocks(); ++i)
    {
      if (across)
        out << 'C' << i << ':';
      else
        out << "Component " << i << std::endl;
      this->components[i].print(out, precision, scientific, across);
    }
}



template <typename Number>
void
BlockVector<Number>::block_write(std::ostream &out) const
{
  for (size_type i = 0; i < this->n_blocks(); ++i)
    this->components[i].block_write(out);
}



template <typename Number>
void
BlockVector<Number>::block_read(std::istream &in)
{
  for (size_type i = 0; i < this->n_blocks(); ++i)
    this->components[i].block_read(in);
}



DEAL_II_NAMESPACE_CLOSE

#endif
