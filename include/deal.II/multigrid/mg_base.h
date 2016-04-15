// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2016 by the deal.II authors
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

#ifndef dealii__mg_base_h
#define dealii__mg_base_h

/*
 * This file contains some abstract base classes
 * used by Multigrid.
 */

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mg */
/*@{*/


/**
 * Multilevel matrix base. This class sets up the interface needed by
 * multilevel algorithms. It has no relation to the actual matrix type and
 * takes the vector class as only template argument.
 *
 * Usually, the derived class MGMatrix, operating on an MGLevelObject of
 * matrices will be sufficient for applications.
 *
 * @author Guido Kanschat, 2002
 */
template <typename VectorType>
class MGMatrixBase : public Subscriptor
{
public:
  /*
   * Virtual destructor.
   */
  virtual ~MGMatrixBase();

  /**
   * Matrix-vector-multiplication on a certain level.
   */
  virtual void vmult (const unsigned int level,
                      VectorType         &dst,
                      const VectorType   &src) const = 0;

  /**
   * Adding matrix-vector-multiplication on a certain level.
   */
  virtual void vmult_add (const unsigned int level,
                          VectorType         &dst,
                          const VectorType   &src) const = 0;

  /**
   * Transpose matrix-vector-multiplication on a certain level.
   */
  virtual void Tvmult (const unsigned int level,
                       VectorType         &dst,
                       const VectorType   &src) const = 0;

  /**
   * Adding transpose matrix-vector-multiplication on a certain level.
   */
  virtual void Tvmult_add (const unsigned int level,
                           VectorType         &dst,
                           const VectorType   &src) const = 0;
};


/**
 * Base class for coarse grid solvers.  This defines the virtual parenthesis
 * operator, being the interface used by multigrid methods. Any implementation
 * will be done by derived classes.
 *
 * @author Guido Kanschat, 2002
 */
template <typename VectorType>
class MGCoarseGridBase : public Subscriptor
{
public:
  /**
   * Virtual destructor.
   */
  virtual ~MGCoarseGridBase ();

  /**
   * Solution operator.
   */
  virtual void operator() (const unsigned int level,
                           VectorType         &dst,
                           const VectorType   &src) const = 0;
};


/**
 * Base class used to declare the operations needed by a concrete class
 * implementing prolongation and restriction of vectors in the multigrid
 * context. This class is abstract and has no implementation of these
 * operations.
 *
 * There are several derived classes, reflecting the fact that vector types
 * and numbering of the fine-grid discretization and of the multi-level
 * implementation are independent.
 *
 * If you use multigrid for a single PDE or for your complete system of
 * equations, you will use MGTransferPrebuilt together with Multigrid. The
 * vector types used on the fine grid as well as for the multilevel operations
 * may be Vector or BlockVector. In both cases, MGTransferPrebuilt will
 * operate on all components of the solution.
 *
 * @note For the following, it is important to realize the difference between
 * a solution
 * @ref GlossComponent "component"
 * and a solution
 * @ref GlossBlock "block".
 * The distinction only applies if vector valued elements are used, but is
 * quite important then. This is reflected in the fact that it is not possible
 * right now to use transfer classes based on MGTransferComponentBase for
 * genuine vector valued elements, but descendants of MGTransferBlockBase
 * would have to be applied. In the following text, we will use the term
 * <em>block</em>, but remark that it might refer to components as well.
 *
 * @todo update the following documentation, since it does not reflect the
 * latest changes in structure.
 *
 * For mixed systems, it may be required to do multigrid only for a single
 * component or for some components. The classes MGTransferSelect and
 * MGTransferBlock handle these cases.
 *
 * MGTransferSelect is used if you use multigrid (on Vector objects) for a
 * single component, possibly grouped using <tt>mg_target_component</tt>.
 *
 * The class MGTransferBlock handles the case where your multigrid method
 * operates on BlockVector objects. These can contain all or a consecutive set
 * of the blocks of the complete system. Since most smoothers cannot operate
 * on block structures, it is not clear whether this case is really useful.
 * Therefore, a tested implementation of this case will be supplied when
 * needed.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2002, 2007
 */
template <typename VectorType>
class MGTransferBase : public Subscriptor
{
public:
  /**
   * Destructor. Does nothing here, but needs to be declared virtual anyway.
   */
  virtual ~MGTransferBase();

  /**
   * Prolongate a vector from level <tt>to_level-1</tt> to level
   * <tt>to_level</tt>. The previous content of <tt>dst</tt> is overwritten.
   *
   * @arg src is a vector with as many elements as there are degrees of
   * freedom on the coarser level involved.
   *
   * @arg dst has as many elements as there are degrees of freedom on the
   * finer level.
   */
  virtual void prolongate (const unsigned int to_level,
                           VectorType         &dst,
                           const VectorType   &src) const = 0;

  /**
   * Restrict a vector from level <tt>from_level</tt> to level
   * <tt>from_level-1</tt> and add this restriction to <tt>dst</tt>. If the
   * region covered by cells on level <tt>from_level</tt> is smaller than that
   * of level <tt>from_level-1</tt> (local refinement), then some degrees of
   * freedom in <tt>dst</tt> are active and will not be altered. For the other
   * degrees of freedom, the result of the restriction is added.
   *
   * @arg src is a vector with as many elements as there are degrees of
   * freedom on the finer level
   *
   * @arg dst has as many elements as there are degrees of freedom on the
   * coarser level.
   *
   */
  virtual void restrict_and_add (const unsigned int from_level,
                                 VectorType         &dst,
                                 const VectorType   &src) const = 0;
};



/**
 * Base class for multigrid smoothers. Does nothing but defining the interface
 * used by multigrid methods.
 *
 * @author Guido Kanschat, 2002
 */
template <typename VectorType>
class MGSmootherBase : public Subscriptor
{
public:
  /**
   * Virtual destructor.
   */
  virtual ~MGSmootherBase();
  /**
   * Release matrices.
   */
  virtual void clear() = 0;

  /**
   * Smoothing function. This is the function used in multigrid methods.
   */
  virtual void smooth (const unsigned int level,
                       VectorType         &u,
                       const VectorType   &rhs) const = 0;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
