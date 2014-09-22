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

#ifndef __deal2__mg_smoother_h
#define __deal2__mg_smoother_h


#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/base/mg_level_object.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/*
 * MGSmootherBase is defined in mg_base.h
 */

/*!@addtogroup mg */
/*@{*/

/**
 * A base class for smoother handling information on smoothing.
 * While not adding to the abstract interface in MGSmootherBase, this
 * class stores information on the number and type of smoothing steps,
 * which in turn can be used by a derived class.
 *
 * @author Guido Kanschat 2009
 */
template <class VECTOR>
class MGSmoother : public MGSmootherBase<VECTOR>
{
public:
  /**
   * Constructor. Sets smoothing
   * parameters and creates a
   * private GrowingVectorMemory
   * object to be used to retrieve vectors.
   */
  MGSmoother(const unsigned int steps = 1,
             const bool variable = false,
             const bool symmetric = false,
             const bool transpose = false);

  /**
   * @deprecated Since
   * GrowingVectorMemory now uses a
   * joint memory pool, it is
   * recommended to use the
   * constructor without the memory
   * object.
   *
   * Constructor. Sets memory and
   * smoothing parameters.
   */
  MGSmoother(VectorMemory<VECTOR> &mem,
             const unsigned int steps = 1,
             const bool variable = false,
             const bool symmetric = false,
             const bool transpose = false) DEAL_II_DEPRECATED;

  /**
   * Modify the number of smoothing
   * steps on finest level.
   */
  void set_steps (const unsigned int);

  /**
   * Switch on/off variable
   * smoothing.
   */
  void set_variable (const bool);

  /**
   * Switch on/off symmetric
   * smoothing.
   */
  void set_symmetric (const bool);

  /**
   * Switch on/off transposed
   * smoothing. The effect is
   * overriden by set_symmetric().
   */
  void set_transpose (const bool);

  /**
   * Set @p debug to a nonzero value
   * to get debug information
   * logged to @p deallog. Increase
   * to get more information
   */
  void set_debug (const unsigned int level);

private:
  /**
   * The memory object to be used
   * if none is given to the
   * constructor.
   */
  GrowingVectorMemory<VECTOR> my_memory;

protected:
  /**
   * Number of smoothing steps on
   * the finest level. If no
   * #variable smoothing is chosen,
   * this is the number of steps on
   * all levels.
   */
  unsigned int steps;

  /**
   * Variable smoothing: double the
   * number of smoothing steps
   * whenever going to the next
   * coarser level
   */
  bool variable;

  /**
   * Symmetric smoothing: in the
   * smoothing iteration, alternate
   * between the relaxation method
   * and its transpose.
   */
  bool symmetric;

  /**
   * Use the transpose of the
   * relaxation method instead of
   * the method itself. This has no
   * effect if #symmetric smoothing
   * is chosen.
   */
  bool transpose;

  /**
   * Output debugging information
   * to @p deallog if this is
   * nonzero.
   */
  unsigned int debug;

  /**
   * Memory for auxiliary vectors.
   */
  SmartPointer<VectorMemory<VECTOR>, MGSmoother<VECTOR> > mem;
};


/**
 * Smoother doing nothing. This class is not useful for many applications other
 * than for testing some multigrid procedures. Also some applications might
 * get convergence without smoothing and then this class brings you the
 * cheapest possible multigrid.
 *
 * @author Guido Kanschat, 1999, 2002
 */
template <class VECTOR>
class MGSmootherIdentity : public MGSmootherBase<VECTOR>
{
public:
  /**
   * Implementation of the
   * interface for @p Multigrid.
   * This function does nothing,
   * which by comparison with the
   * definition of this function
   * means that the the smoothing
   * operator equals the null
   * operator.
   */
  virtual void smooth (const unsigned int level,
                       VECTOR            &u,
                       const VECTOR      &rhs) const;
  virtual void clear ();
};


namespace mg
{
  /**
   * Smoother using relaxation classes.
   *
   * A relaxation class is an object that has two member functions,
   * @code
   * void  step(VECTOR& x, const VECTOR& d) const;
   * void Tstep(VECTOR& x, const VECTOR& d) const;
   * @endcode
   * performing one step of the smoothing scheme.
   *
   * This class performs smoothing on each level. The operation can be
   * controlled by several parameters. First, the relaxation parameter
   * @p omega is used in the underlying relaxation method. @p steps is
   * the number of relaxation steps on the finest level (on all levels
   * if @p variable is off). If @p variable is @p true, the number of
   * smoothing steps is doubled on each coarser level. This results in a
   * method having the complexity of the W-cycle, but saving grid
   * transfers. This is the method proposed by Bramble at al.
   *
   * The option @p symmetric switches on alternating between the
   * smoother and its transpose in each step as proposed by Bramble.
   *
   * @p transpose uses the transposed smoothing operation using
   * <tt>Tstep</tt> instead of the regular <tt>step</tt> of the
   * relaxation scheme.
   *
   * If you are using block matrices, the second @p initialize function
   * offers the possibility to extract a single block for smoothing. In
   * this case, the multigrid method must be used only with the vector
   * associated to that single block.
   *
   * @author Guido Kanschat,
   * @date 2003, 2009, 2010
   */
  template<class RELAX, class VECTOR>
  class SmootherRelaxation : public MGLevelObject<RELAX>, public MGSmoother<VECTOR>
  {
  public:
    /**
     * Constructor. Sets memory and
     * smoothing parameters.
     */
    SmootherRelaxation(const unsigned int steps = 1,
                       const bool variable = false,
                       const bool symmetric = false,
                       const bool transpose = false);

    /**
     * Initialize for matrices. This
     * function initializes the
     * smoothing operator with the
     * same smoother for each level.
     *
     * @p additional_data is an
     * object of type
     * @p RELAX::AdditionalData and
     * is handed to the
     * initialization function of the
     * relaxation method.
     */
    template <class MATRIX2>
    void initialize (const MGLevelObject<MATRIX2> &matrices,
                     const typename RELAX::AdditionalData &additional_data = typename RELAX::AdditionalData());

    /**
     * Initialize matrices and
     * additional data for each
     * level.
     *
     * If minimal or maximal level of
     * the two objects differ, the
     * greatest common range is
     * utilized. This way, smoothing
     * can be restricted to certain
     * levels even if the matrix was
     * generated for all levels.
     */
    template <class MATRIX2, class DATA>
    void initialize (const MGLevelObject<MATRIX2> &matrices,
                     const MGLevelObject<DATA> &additional_data);
    /**
     * Initialize for matrix
     * blocks. This function
     * initializes the smoothing
     * operator with the same
     * smoother for each level.
     *
     * @p additional_data is an
     * object of type
     * @p RELAX::AdditionalData and
     * is handed to the
     * initialization function of the
     * relaxation method.
     */
//     template <class MATRIX2>
//     void initialize (const MGLevelObject<MatrixBlock<MATRIX2> >& matrices,
//                   const typename RELAX::AdditionalData & additional_data = typename RELAX::AdditionalData());

    /**
     * Empty all vectors.
     */
    void clear ();

    /**
     * The actual smoothing method.
     */
    virtual void smooth (const unsigned int level,
                         VECTOR            &u,
                         const VECTOR      &rhs) const;

    /**
     * Memory used by this object.
     */
    std::size_t memory_consumption () const;
  };
}

/**
 * Smoother using relaxation classes.
 *
 * A relaxation class is an object that has two member functions,
 * @code
 * void  step(VECTOR& x, const VECTOR& d) const;
 * void Tstep(VECTOR& x, const VECTOR& d) const;
 * @endcode
 * performing one step of the smoothing scheme.
 *
 * This class performs smoothing on each level. The operation can be
 * controlled by several parameters. First, the relaxation parameter
 * @p omega is used in the underlying relaxation method. @p steps is
 * the number of relaxation steps on the finest level (on all levels
 * if @p variable is off). If @p variable is @p true, the number of
 * smoothing steps is doubled on each coarser level. This results in a
 * method having the complexity of the W-cycle, but saving grid
 * transfers. This is the method proposed by Bramble at al.
 *
 * The option @p symmetric switches on alternating between the
 * smoother and its transpose in each step as proposed by Bramble.
 *
 * @p transpose uses the transposed smoothing operation using
 * <tt>Tstep</tt> instead of the regular <tt>step</tt> of the
 * relaxation scheme.
 *
 * If you are using block matrices, the second @p initialize function
 * offers the possibility to extract a single block for smoothing. In
 * this case, the multigrid method must be used only with the vector
 * associated to that single block.
 *
 * The library contains instantiation for <tt>SparseMatrix<.></tt> and
 * <tt>Vector<.></tt>, where the template arguments are all combinations of
 * @p float and @p double. Additional instantiations may be created
 * by including the file mg_smoother.templates.h.
 *
 * @note If '--enable-mgcompatibility' was used on configuring
 * deal.II, this class behaves like MGSmootherPrecondition.
 *
 * @author Guido Kanschat, 2003
 */
template<class MATRIX, class RELAX, class VECTOR>
class MGSmootherRelaxation : public MGSmoother<VECTOR>
{
public:
  /**
   * Constructor. Sets smoothing parameters.
   */
  MGSmootherRelaxation(const unsigned int steps = 1,
                       const bool variable = false,
                       const bool symmetric = false,
                       const bool transpose = false);

  /**
   * Constructor. Sets memory and
   * smoothing parameters.
   *
   * @deprecated Use the constructor without the vector memory
   * object
   */
  MGSmootherRelaxation(VectorMemory<VECTOR> &mem,
                       const unsigned int steps = 1,
                       const bool variable = false,
                       const bool symmetric = false,
                       const bool transpose = false) DEAL_II_DEPRECATED;

  /**
   * Initialize for matrices. This
   * function stores pointers to
   * the level matrices and
   * initializes the smoothing
   * operator with the same smoother
   * for each level.
   *
   * @p additional_data is an
   * object of type
   * @p RELAX::AdditionalData and
   * is handed to the
   * initialization function of the
   * relaxation method.
   */
  template <class MATRIX2>
  void initialize (const MGLevelObject<MATRIX2> &matrices,
                   const typename RELAX::AdditionalData &additional_data = typename RELAX::AdditionalData());

  /**
   * Initialize for matrices. This
   * function stores pointers to
   * the level matrices and
   * initializes the smoothing
   * operator with the according
   * smoother for each level.
   *
   * @p additional_data is an
   * object of type
   * @p RELAX::AdditionalData and
   * is handed to the
   * initialization function of the
   * relaxation method.
   */
  template <class MATRIX2, class DATA>
  void initialize (const MGLevelObject<MATRIX2> &matrices,
                   const MGLevelObject<DATA> &additional_data);

  /**
   * Initialize for single blocks
   * of matrices. Of this block
   * matrix, the block indicated by
   * @p block_row and
   * @p block_col is selected on
   * each level.  This function
   * stores pointers to the level
   * matrices and initializes the
   * smoothing operator with the
   * same smoother for each
   * level.
   *
   * @p additional_data is an
   * object of type
   * @p RELAX::AdditionalData and
   * is handed to the
   * initialization function of the
   * relaxation method.
   */
  template <class MATRIX2, class DATA>
  void initialize (const MGLevelObject<MATRIX2> &matrices,
                   const DATA &additional_data,
                   const unsigned int block_row,
                   const unsigned int block_col);

  /**
   * Initialize for single blocks
   * of matrices. Of this block
   * matrix, the block indicated by
   * @p block_row and
   * @p block_col is selected on
   * each level.  This function
   * stores pointers to the level
   * matrices and initializes the
   * smoothing operator with the
   * according smoother for each
   * level.
   *
   * @p additional_data is an
   * object of type
   * @p RELAX::AdditionalData and
   * is handed to the
   * initialization function of the
   * relaxation method.
   */
  template <class MATRIX2, class DATA>
  void initialize (const MGLevelObject<MATRIX2> &matrices,
                   const MGLevelObject<DATA> &additional_data,
                   const unsigned int block_row,
                   const unsigned int block_col);

  /**
   * Empty all vectors.
   */
  void clear ();

  /**
   * The actual smoothing method.
   */
  virtual void smooth (const unsigned int level,
                       VECTOR            &u,
                       const VECTOR      &rhs) const;

  /**
   * Object containing relaxation
   * methods.
   */
  MGLevelObject<RELAX> smoothers;

  /**
   * Memory used by this object.
   */
  std::size_t memory_consumption () const;


private:
  /**
   * Pointer to the matrices.
   */
  MGLevelObject<PointerMatrix<MATRIX, VECTOR> > matrices;

};



/**
 * Smoother using preconditioner classes.
 *
 * This class performs smoothing on each level. The operation can be
 * controlled by several parameters. First, the relaxation parameter
 * @p omega is used in the underlying relaxation method. @p steps is
 * the number of relaxation steps on the finest level (on all levels
 * if @p variable is off). If @p variable is @p true, the number of
 * smoothing steps is doubled on each coarser level. This results in a
 * method having the complexity of the W-cycle, but saving grid
 * transfers. This is the method proposed by Bramble at al.
 *
 * The option @p symmetric switches on alternating between the
 * smoother and its transpose in each step as proposed by Bramble.
 *
 * @p transpose uses the transposed smoothing operation using
 * <tt>Tvmult</tt> instead of the regular <tt>vmult</tt> of the
 * relaxation scheme.
 *
 * If you are using block matrices, the second @p initialize function
 * offers the possibility to extract a single block for smoothing. In
 * this case, the multigrid method must be used only with the vector
 * associated to that single block.
 *
 * The library contains instantiation for <tt>SparseMatrix<.></tt> and
 * <tt>Vector<.></tt>, where the template arguments are all combinations of
 * @p float and @p double. Additional instantiations may be created
 * by including the file mg_smoother.templates.h.
 *
 * @author Guido Kanschat, 2009
 */
template<class MATRIX, class PRECONDITIONER, class VECTOR>
class MGSmootherPrecondition : public MGSmoother<VECTOR>
{
public:
  /**
   * Constructor. Sets smoothing parameters.
   */
  MGSmootherPrecondition(const unsigned int steps = 1,
                         const bool variable = false,
                         const bool symmetric = false,
                         const bool transpose = false);

  /**
   * Constructor. Sets memory and
   * smoothing parameters.
   *
   * @deprecated Use the constructor without the vector memory
   * object
   */
  MGSmootherPrecondition(VectorMemory<VECTOR> &mem,
                         const unsigned int steps = 1,
                         const bool variable = false,
                         const bool symmetric = false,
                         const bool transpose = false) DEAL_II_DEPRECATED;

  /**
   * Initialize for matrices. This
   * function stores pointers to
   * the level matrices and
   * initializes the smoothing
   * operator with the same smoother
   * for each level.
   *
   * @p additional_data is an
   * object of type
   * @p PRECONDITIONER::AdditionalData and
   * is handed to the
   * initialization function of the
   * relaxation method.
   */
  template <class MATRIX2>
  void initialize (const MGLevelObject<MATRIX2> &matrices,
                   const typename PRECONDITIONER::AdditionalData &additional_data = typename PRECONDITIONER::AdditionalData());

  /**
   * Initialize for matrices. This
   * function stores pointers to
   * the level matrices and
   * initializes the smoothing
   * operator with the according
   * smoother for each level.
   *
   * @p additional_data is an
   * object of type
   * @p PRECONDITIONER::AdditionalData and
   * is handed to the
   * initialization function of the
   * relaxation method.
   */
  template <class MATRIX2, class DATA>
  void initialize (const MGLevelObject<MATRIX2> &matrices,
                   const MGLevelObject<DATA> &additional_data);

  /**
   * Initialize for single blocks
   * of matrices. Of this block
   * matrix, the block indicated by
   * @p block_row and
   * @p block_col is selected on
   * each level.  This function
   * stores pointers to the level
   * matrices and initializes the
   * smoothing operator with the
   * same smoother for each
   * level.
   *
   * @p additional_data is an
   * object of type
   * @p PRECONDITIONER::AdditionalData and
   * is handed to the
   * initialization function of the
   * relaxation method.
   */
  template <class MATRIX2, class DATA>
  void initialize (const MGLevelObject<MATRIX2> &matrices,
                   const DATA &additional_data,
                   const unsigned int block_row,
                   const unsigned int block_col);

  /**
   * Initialize for single blocks
   * of matrices. Of this block
   * matrix, the block indicated by
   * @p block_row and
   * @p block_col is selected on
   * each level.  This function
   * stores pointers to the level
   * matrices and initializes the
   * smoothing operator with the
   * according smoother for each
   * level.
   *
   * @p additional_data is an
   * object of type
   * @p PRECONDITIONER::AdditionalData and
   * is handed to the
   * initialization function of the
   * relaxation method.
   */
  template <class MATRIX2, class DATA>
  void initialize (const MGLevelObject<MATRIX2> &matrices,
                   const MGLevelObject<DATA> &additional_data,
                   const unsigned int block_row,
                   const unsigned int block_col);

  /**
   * Empty all vectors.
   */
  void clear ();

  /**
   * The actual smoothing method.
   */
  virtual void smooth (const unsigned int level,
                       VECTOR            &u,
                       const VECTOR      &rhs) const;

  /**
   * Object containing relaxation
   * methods.
   */
  MGLevelObject<PRECONDITIONER> smoothers;

  /**
   * Memory used by this object.
   */
  std::size_t memory_consumption () const;


private:
  /**
   * Pointer to the matrices.
   */
  MGLevelObject<PointerMatrix<MATRIX, VECTOR> > matrices;

};

/*@}*/

/* ------------------------------- Inline functions -------------------------- */

#ifndef DOXYGEN

template <class VECTOR>
inline void
MGSmootherIdentity<VECTOR>::smooth (
  const unsigned int, VECTOR &,
  const VECTOR &) const
{}

template <class VECTOR>
inline void
MGSmootherIdentity<VECTOR>::clear ()
{}

//---------------------------------------------------------------------------

template <class VECTOR>
inline
MGSmoother<VECTOR>::MGSmoother(
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose)
  :
  steps(steps),
  variable(variable),
  symmetric(symmetric),
  transpose(transpose),
  debug(0),
  mem(&my_memory)
{}


template <class VECTOR>
inline
MGSmoother<VECTOR>::MGSmoother(
  VectorMemory<VECTOR> &mem,
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose)
  :
  steps(steps),
  variable(variable),
  symmetric(symmetric),
  transpose(transpose),
  debug(0),
  mem(&mem)
{}


template <class VECTOR>
inline void
MGSmoother<VECTOR>::set_steps (const unsigned int s)
{
  steps = s;
}


template <class VECTOR>
inline void
MGSmoother<VECTOR>::set_debug (const unsigned int s)
{
  debug = s;
}


template <class VECTOR>
inline void
MGSmoother<VECTOR>::set_variable (const bool flag)
{
  variable = flag;
}


template <class VECTOR>
inline void
MGSmoother<VECTOR>::set_symmetric (const bool flag)
{
  symmetric = flag;
}


template <class VECTOR>
inline void
MGSmoother<VECTOR>::set_transpose (const bool flag)
{
  transpose = flag;
}

//----------------------------------------------------------------------//

namespace mg
{
  template <class RELAX, class VECTOR>
  inline
  SmootherRelaxation<RELAX, VECTOR>::SmootherRelaxation(
    const unsigned int steps,
    const bool variable,
    const bool symmetric,
    const bool transpose)
    : MGSmoother<VECTOR>(steps, variable, symmetric, transpose)
  {}


  template <class RELAX, class VECTOR>
  inline void
  SmootherRelaxation<RELAX, VECTOR>::clear ()
  {
    MGLevelObject<RELAX>::clear();
  }


  template <class RELAX, class VECTOR>
  template <class MATRIX2>
  inline void
  SmootherRelaxation<RELAX, VECTOR>::initialize (
    const MGLevelObject<MATRIX2> &m,
    const typename RELAX::AdditionalData &data)
  {
    const unsigned int min = m.min_level();
    const unsigned int max = m.max_level();

    this->resize(min, max);

    for (unsigned int i=min; i<=max; ++i)
      (*this)[i].initialize(m[i], data);
  }


  template <class RELAX, class VECTOR>
  template <class MATRIX2, class DATA>
  inline void
  SmootherRelaxation<RELAX, VECTOR>::initialize (
    const MGLevelObject<MATRIX2> &m,
    const MGLevelObject<DATA> &data)
  {
    const unsigned int min = std::max(m.min_level(), data.min_level());
    const unsigned int max = std::min(m.max_level(), data.max_level());

    this->resize(min, max);

    for (unsigned int i=min; i<=max; ++i)
      (*this)[i].initialize(m[i], data[i]);
  }


  template <class RELAX, class VECTOR>
  inline void
  SmootherRelaxation<RELAX, VECTOR>::smooth(
    const unsigned int level,
    VECTOR &u,
    const VECTOR &rhs) const
  {
    unsigned int maxlevel = this->max_level();
    unsigned int steps2 = this->steps;

    if (this->variable)
      steps2 *= (1<<(maxlevel-level));

    bool T = this->transpose;
    if (this->symmetric && (steps2 % 2 == 0))
      T = false;
    if (this->debug > 0)
      deallog << 'S' << level << ' ';

    for (unsigned int i=0; i<steps2; ++i)
      {
        if (T)
          (*this)[level].Tstep(u, rhs);
        else
          (*this)[level].step(u, rhs);
        if (this->symmetric)
          T = !T;
      }
  }


  template <class RELAX, class VECTOR>
  inline
  std::size_t
  SmootherRelaxation<RELAX, VECTOR>::
  memory_consumption () const
  {
    return sizeof(*this)
           -sizeof(MGLevelObject<RELAX>)
           + MGLevelObject<RELAX>::memory_consumption()
           + this->mem->memory_consumption();
  }
}


//----------------------------------------------------------------------//

template <class MATRIX, class RELAX, class VECTOR>
inline
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::MGSmootherRelaxation(
  VectorMemory<VECTOR> &mem,
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose)
  :
  MGSmoother<VECTOR>(mem, steps, variable, symmetric, transpose)
{}



template <class MATRIX, class RELAX, class VECTOR>
inline
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::MGSmootherRelaxation(
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose)
  :
  MGSmoother<VECTOR>(steps, variable, symmetric, transpose)
{}



template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::clear ()
{
  smoothers.clear();

  unsigned int i=matrices.min_level(),
               max_level=matrices.max_level();
  for (; i<=max_level; ++i)
    matrices[i]=0;
}


template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2> &m,
  const typename RELAX::AdditionalData &data)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min; i<=max; ++i)
    {
      matrices[i] = &m[i];
      smoothers[i].initialize(m[i], data);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2> &m,
  const MGLevelObject<DATA> &data)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  Assert (data.min_level() == min,
          ExcDimensionMismatch(data.min_level(), min));
  Assert (data.max_level() == max,
          ExcDimensionMismatch(data.max_level(), max));

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min; i<=max; ++i)
    {
      matrices[i] = &m[i];
      smoothers[i].initialize(m[i], data[i]);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2> &m,
  const DATA &data,
  const unsigned int row,
  const unsigned int col)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min; i<=max; ++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2> &m,
  const MGLevelObject<DATA> &data,
  const unsigned int row,
  const unsigned int col)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  Assert (data.min_level() == min,
          ExcDimensionMismatch(data.min_level(), min));
  Assert (data.max_level() == max,
          ExcDimensionMismatch(data.max_level(), max));

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min; i<=max; ++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data[i]);
    }
}


template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::smooth(
  const unsigned int level,
  VECTOR &u,
  const VECTOR &rhs) const
{
  unsigned int maxlevel = smoothers.max_level();
  unsigned int steps2 = this->steps;

  if (this->variable)
    steps2 *= (1<<(maxlevel-level));

  bool T = this->transpose;
  if (this->symmetric && (steps2 % 2 == 0))
    T = false;
  if (this->debug > 0)
    deallog << 'S' << level << ' ';

  for (unsigned int i=0; i<steps2; ++i)
    {
      if (T)
        smoothers[level].Tstep(u, rhs);
      else
        smoothers[level].step(u, rhs);
      if (this->symmetric)
        T = !T;
    }
}



template <class MATRIX, class RELAX, class VECTOR>
inline
std::size_t
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::
memory_consumption () const
{
  return sizeof(*this)
         + matrices.memory_consumption()
         + smoothers.memory_consumption()
         + this->mem->memory_consumption();
}


//----------------------------------------------------------------------//

template <class MATRIX, class PRECONDITIONER, class VECTOR>
inline
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::MGSmootherPrecondition(
  VectorMemory<VECTOR> &mem,
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose)
  :
  MGSmoother<VECTOR>(mem, steps, variable, symmetric, transpose)
{}



template <class MATRIX, class PRECONDITIONER, class VECTOR>
inline
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::MGSmootherPrecondition(
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose)
  :
  MGSmoother<VECTOR>(steps, variable, symmetric, transpose)
{}



template <class MATRIX, class PRECONDITIONER, class VECTOR>
inline void
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::clear ()
{
  smoothers.clear();

  unsigned int i=matrices.min_level(),
               max_level=matrices.max_level();
  for (; i<=max_level; ++i)
    matrices[i]=0;
}



template <class MATRIX, class PRECONDITIONER, class VECTOR>
template <class MATRIX2>
inline void
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::initialize (
  const MGLevelObject<MATRIX2> &m,
  const typename PRECONDITIONER::AdditionalData &data)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min; i<=max; ++i)
    {
      matrices[i] = &m[i];
      smoothers[i].initialize(m[i], data);
    }
}



template <class MATRIX, class PRECONDITIONER, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::initialize (
  const MGLevelObject<MATRIX2> &m,
  const MGLevelObject<DATA> &data)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  Assert (data.min_level() == min,
          ExcDimensionMismatch(data.min_level(), min));
  Assert (data.max_level() == max,
          ExcDimensionMismatch(data.max_level(), max));

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min; i<=max; ++i)
    {
      matrices[i] = &m[i];
      smoothers[i].initialize(m[i], data[i]);
    }
}



template <class MATRIX, class PRECONDITIONER, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::initialize (
  const MGLevelObject<MATRIX2> &m,
  const DATA &data,
  const unsigned int row,
  const unsigned int col)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min; i<=max; ++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data);
    }
}



template <class MATRIX, class PRECONDITIONER, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::initialize (
  const MGLevelObject<MATRIX2> &m,
  const MGLevelObject<DATA> &data,
  const unsigned int row,
  const unsigned int col)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  Assert (data.min_level() == min,
          ExcDimensionMismatch(data.min_level(), min));
  Assert (data.max_level() == max,
          ExcDimensionMismatch(data.max_level(), max));

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min; i<=max; ++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data[i]);
    }
}



template <class MATRIX, class PRECONDITIONER, class VECTOR>
inline void
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::smooth(
  const unsigned int level,
  VECTOR &u,
  const VECTOR &rhs) const
{
  unsigned int maxlevel = matrices.max_level();
  unsigned int steps2 = this->steps;

  if (this->variable)
    steps2 *= (1<<(maxlevel-level));

  typename VectorMemory<VECTOR>::Pointer r(*this->mem);
  typename VectorMemory<VECTOR>::Pointer d(*this->mem);

  r->reinit(u,true);
  d->reinit(u,true);

  bool T = this->transpose;
  if (this->symmetric && (steps2 % 2 == 0))
    T = false;
  if (this->debug > 0)
    deallog << 'S' << level << ' ';

  for (unsigned int i=0; i<steps2; ++i)
    {
      if (T)
        {
          if (this->debug > 0)
            deallog << 'T';
          if (i == 0 && u.all_zero())
            *r = rhs;
          else
            {
              matrices[level].Tvmult(*r,u);
              r->sadd(-1.,1.,rhs);
            }
          if (this->debug > 2)
            deallog << ' ' << r->l2_norm() << ' ';
          smoothers[level].Tvmult(*d, *r);
          if (this->debug > 1)
            deallog << ' ' << d->l2_norm() << ' ';
        }
      else
        {
          if (this->debug > 0)
            deallog << 'N';
          if (i == 0 && u.all_zero())
            *r = rhs;
          else
            {
              matrices[level].vmult(*r,u);
              r->sadd(-1.,rhs);
            }
          if (this->debug > 2)
            deallog << ' ' << r->l2_norm() << ' ';
          smoothers[level].vmult(*d, *r);
          if (this->debug > 1)
            deallog << ' ' << d->l2_norm() << ' ';
        }
      u += *d;
      if (this->symmetric)
        T = !T;
    }
  if (this->debug > 0)
    deallog << std::endl;
}



template <class MATRIX, class PRECONDITIONER, class VECTOR>
inline
std::size_t
MGSmootherPrecondition<MATRIX, PRECONDITIONER, VECTOR>::
memory_consumption () const
{
  return sizeof(*this)
         + matrices.memory_consumption()
         + smoothers.memory_consumption()
         + this->mem->memory_consumption();
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
