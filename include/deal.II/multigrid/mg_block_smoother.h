// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_mg_block_smoother_h
#define dealii_mg_block_smoother_h


#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_smoother.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/*
 * MGSmootherBase is defined in mg_base.h
 */

/*!@addtogroup mg */
/*@{*/

/**
 * General smoother class for block vectors. This class gives complete freedom
 * to the choice of a block smoother by being initialized with a matrix and a
 * smoother object. Therefore, the smoother object for each level must be
 * constructed by hand.
 *
 * @author Guido Kanschat, 2005
 */
template <typename MatrixType, class RelaxationType, typename number>
class MGSmootherBlock : public MGSmoother<BlockVector<number>>
{
public:
  /**
   * Constructor.
   */
  MGSmootherBlock(const unsigned int steps     = 1,
                  const bool         variable  = false,
                  const bool         symmetric = false,
                  const bool         transpose = false,
                  const bool         reverse   = false);

  /**
   * Initialize for matrices. The parameter <tt>matrices</tt> can be any
   * object having functions <tt>get_minlevel()</tt> and
   * <tt>get_maxlevel()</tt> as well as an <tt>operator[]</tt> returning a
   * reference to @p MatrixType.
   *
   * The same convention is used for the parameter <tt>smoothers</tt>, such
   * that <tt>operator[]</tt> returns the object doing the block-smoothing on
   * a single level.
   *
   * This function stores pointers to the level matrices and smoothing
   * operator for each level.
   */
  template <class MGMatrixType, class MGRelaxationType>
  void
  initialize(const MGMatrixType &matrices, const MGRelaxationType &smoothers);

  /**
   * Empty all vectors.
   */
  void
  clear();

  /**
   * Switch on/off reversed. This is mutually exclusive with transpose().
   */
  void
  set_reverse(const bool);

  /**
   * Implementation of the interface for @p Multigrid. This function does
   * nothing, which by comparison with the definition of this function means
   * that the smoothing operator equals the null operator.
   */
  virtual void
  smooth(const unsigned int         level,
         BlockVector<number> &      u,
         const BlockVector<number> &rhs) const;

  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * Pointer to the matrices.
   */
  MGLevelObject<LinearOperator<BlockVector<number>>> matrices;

  /**
   * Pointer to the matrices.
   */
  MGLevelObject<LinearOperator<BlockVector<number>>> smoothers;

  /**
   * Reverse?
   */
  bool reverse;

  /**
   * Memory for auxiliary vectors.
   */
  SmartPointer<VectorMemory<BlockVector<number>>,
               MGSmootherBlock<MatrixType, RelaxationType, number>>
    mem;
};

/**@}*/

//---------------------------------------------------------------------------

#ifndef DOXYGEN

template <typename MatrixType, class RelaxationType, typename number>
inline MGSmootherBlock<MatrixType, RelaxationType, number>::MGSmootherBlock(
  const unsigned int steps,
  const bool         variable,
  const bool         symmetric,
  const bool         transpose,
  const bool         reverse)
  : MGSmoother<BlockVector<number>>(steps, variable, symmetric, transpose)
  , reverse(reverse)
  , mem(&this->vector_memory)
{}


template <typename MatrixType, class RelaxationType, typename number>
inline void
MGSmootherBlock<MatrixType, RelaxationType, number>::clear()
{
  unsigned int i = matrices.min_level(), max_level = matrices.max_level();
  for (; i <= max_level; ++i)
    {
      smoothers[i] = LinearOperator<BlockVector<number>>();
      matrices[i]  = LinearOperator<BlockVector<number>>();
    }
}


template <typename MatrixType, class RelaxationType, typename number>
template <class MGMatrixType, class MGRelaxationType>
inline void
MGSmootherBlock<MatrixType, RelaxationType, number>::initialize(
  const MGMatrixType &    m,
  const MGRelaxationType &s)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      // Workaround: Unfortunately, not every "m[i]" object has a
      // rich enough interface to populate reinit_(domain|range)_vector.
      // Thus, apply an empty LinearOperator exemplar.
      matrices[i] = linear_operator<BlockVector<number>>(
        LinearOperator<BlockVector<number>>(), m[i]);
      smoothers[i] = linear_operator<BlockVector<number>>(matrices[i], s[i]);
    }
}


template <typename MatrixType, class RelaxationType, typename number>
inline void
MGSmootherBlock<MatrixType, RelaxationType, number>::set_reverse(
  const bool flag)
{
  reverse = flag;
}


template <typename MatrixType, class RelaxationType, typename number>
inline std::size_t
MGSmootherBlock<MatrixType, RelaxationType, number>::memory_consumption() const
{
  return sizeof(*this) + matrices.memory_consumption() +
         smoothers.memory_consumption() +
         this->vector_memory.memory_consumption();
}


template <typename MatrixType, class RelaxationType, typename number>
inline void
MGSmootherBlock<MatrixType, RelaxationType, number>::smooth(
  const unsigned int         level,
  BlockVector<number> &      u,
  const BlockVector<number> &rhs) const
{
  LogStream::Prefix prefix("Smooth");

  unsigned int maxlevel = matrices.max_level();
  unsigned int steps2   = this->steps;

  if (this->variable)
    steps2 *= (1 << (maxlevel - level));

  typename VectorMemory<BlockVector<number>>::Pointer r(*this->mem);
  typename VectorMemory<BlockVector<number>>::Pointer d(*this->mem);
  r->reinit(u);
  d->reinit(u);

  bool T = this->transpose;
  if (this->symmetric && (steps2 % 2 == 0))
    T = false;

  for (unsigned int i = 0; i < steps2; ++i)
    {
      if (T)
        {
          matrices[level].vmult(*r, u);
          r->sadd(-1., 1., rhs);
          smoothers[level].Tvmult(*d, *r);
        }
      else
        {
          matrices[level].vmult(*r, u);
          r->sadd(-1., 1., rhs);
          smoothers[level].vmult(*d, *r);
        }
      u += *d;
      if (this->symmetric)
        T = !T;
    }
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
