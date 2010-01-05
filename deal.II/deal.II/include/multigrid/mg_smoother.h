//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_smoother_h
#define __deal2__mg_smoother_h


#include <base/config.h>
#include <base/smartpointer.h>
#include <lac/pointer_matrix.h>
#include <lac/vector_memory.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_level_object.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class MGDoFHandler;


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
				      * Constructor. Sets memory and
				      * smoothing parameters.
				      */
    MGSmoother(VectorMemory<VECTOR>& mem,
	       const unsigned int steps = 1,
	       const bool variable = false,
	       const bool symmetric = false,
	       const bool transpose = false);
    
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
				      * Set #debug to a nonzero value
				      * to get debug information
				      * logged to #deallog. Increase
				      * to get more information
				      */
    void set_debug (const unsigned int level);
    
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

				     /*
				      * Use the transpose of the
				      * relaxation method instead of
				      * the method itself. This has no
				      * effect if #symmetric smoothing
				      * is chosen.
				      */
    bool transpose;
    
				     /**
				      * Output debugging information
				      * to #deallog if this is
				      * nonzero.
				      */
    unsigned int debug;

				     /**
				      * Memory for auxiliary vectors.
				      */
    VectorMemory<VECTOR>& mem;    
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
			 VECTOR&            u,
			 const VECTOR&      rhs) const;
};


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
				      * Constructor. Sets memory and
				      * smoothing parameters.
				      */
    MGSmootherRelaxation(VectorMemory<VECTOR>& mem,
			 const unsigned int steps = 1,
			 const bool variable = false,
			 const bool symmetric = false,
			 const bool transpose = false);

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
    template <class MATRIX2, class DATA>
    void initialize (const MGLevelObject<MATRIX2>& matrices,
		     const DATA& additional_data);

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
    void initialize (const MGLevelObject<MATRIX2>& matrices,
		     const MGLevelObject<DATA>& additional_data);

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
    void initialize (const MGLevelObject<MATRIX2>& matrices,
		     const DATA& additional_data,
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
    void initialize (const MGLevelObject<MATRIX2>& matrices,
		     const MGLevelObject<DATA>& additional_data,
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
			 VECTOR&            u,
			 const VECTOR&      rhs) const;

 				     /**
				      * Object containing relaxation
				      * methods.
				      */
    MGLevelObject<RELAX> smoothers;

    				     /**
				      * Memory used by this object.
				      */
    unsigned int memory_consumption () const;
    

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
template<class MATRIX, class RELAX, class VECTOR>
class MGSmootherPrecondition : public MGSmoother<VECTOR>
{
  public:
				     /**
				      * Constructor. Sets memory and
				      * smoothing parameters.
				      */
    MGSmootherPrecondition(VectorMemory<VECTOR>& mem,
			   const unsigned int steps = 1,
			   const bool variable = false,
			   const bool symmetric = false,
			   const bool transpose = false);

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
    template <class MATRIX2, class DATA>
    void initialize (const MGLevelObject<MATRIX2>& matrices,
		     const DATA& additional_data);

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
    void initialize (const MGLevelObject<MATRIX2>& matrices,
		     const MGLevelObject<DATA>& additional_data);

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
    void initialize (const MGLevelObject<MATRIX2>& matrices,
		     const DATA& additional_data,
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
    void initialize (const MGLevelObject<MATRIX2>& matrices,
		     const MGLevelObject<DATA>& additional_data,
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
			 VECTOR&            u,
			 const VECTOR&      rhs) const;

 				     /**
				      * Object containing relaxation
				      * methods.
				      */
    MGLevelObject<RELAX> smoothers;

    				     /**
				      * Memory used by this object.
				      */
    unsigned int memory_consumption () const;
    

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
  const unsigned int, VECTOR&,
  const VECTOR&) const
{}

//---------------------------------------------------------------------------

template <class VECTOR>
inline
MGSmoother<VECTOR>::MGSmoother(
  VectorMemory<VECTOR>& mem,
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
		mem(mem)
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

template <class MATRIX, class RELAX, class VECTOR>
inline
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::MGSmootherRelaxation(
  VectorMemory<VECTOR>& mem,
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose)
		:
		MGSmoother<VECTOR>(mem, steps, variable, symmetric, transpose)
{}


template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::clear ()
{
  smoothers.clear();
  
  unsigned int i=matrices.get_minlevel(),
       max_level=matrices.get_maxlevel();
  for (; i<=max_level; ++i)
    matrices[i]=0;
}


template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2>& m,
  const DATA& data)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &m[i];
      smoothers[i].initialize(m[i], data);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2>& m,
  const MGLevelObject<DATA>& data)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &m[i];
      smoothers[i].initialize(m[i], data[i]);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2>& m,
  const DATA& data,
  const unsigned int row,
  const unsigned int col)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2>& m,
  const MGLevelObject<DATA>& data,
  const unsigned int row,
  const unsigned int col)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data[i]);
    }
}

#ifdef DEAL_II_MULTIGRID_COMPATIBILITY

template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::smooth(
  const unsigned int level,
  VECTOR& u,
  const VECTOR& rhs) const
{
  unsigned int maxlevel = matrices.get_maxlevel();
  unsigned int steps2 = this->steps;

  if (this->variable)
    steps2 *= (1<<(maxlevel-level));

  VECTOR* r = this->mem.alloc();
  VECTOR* d = this->mem.alloc();
  r->reinit(u);
  d->reinit(u);

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
	  matrices[level].vmult(*r,u);
	  r->sadd(-1.,1.,rhs);
	  if (this->debug > 2)
	    deallog << ' ' << r->l2_norm() << ' ';
	  smoothers[level].Tvmult(*d, *r);
	  if (this->debug > 1)
	    deallog << ' ' << d->l2_norm() << ' ';
	} else {
	  if (this->debug > 0)
	    deallog << 'N';
	  matrices[level].vmult(*r,u);
	  r->sadd(-1.,1.,rhs);
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
  
  this->mem.free(r);
  this->mem.free(d);
}

#else

template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::smooth(
  const unsigned int level,
  VECTOR& u,
  const VECTOR& rhs) const
{
  unsigned int maxlevel = matrices.get_maxlevel();
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

#endif


template <class MATRIX, class RELAX, class VECTOR>
inline unsigned int
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::
memory_consumption () const
{
  return sizeof(*this)
    + matrices.memory_consumption()
    + smoothers.memory_consumption()
    + this->mem.memory_consumption();
}


//----------------------------------------------------------------------//

template <class MATRIX, class RELAX, class VECTOR>
inline
MGSmootherPrecondition<MATRIX, RELAX, VECTOR>::MGSmootherPrecondition(
  VectorMemory<VECTOR>& mem,
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose)
		:
		MGSmoother<VECTOR>(mem, steps, variable, symmetric, transpose)
{}


template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherPrecondition<MATRIX, RELAX, VECTOR>::clear ()
{
  smoothers.clear();
  
  unsigned int i=matrices.get_minlevel(),
       max_level=matrices.get_maxlevel();
  for (; i<=max_level; ++i)
    matrices[i]=0;
}


template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherPrecondition<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2>& m,
  const DATA& data)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &m[i];
      smoothers[i].initialize(m[i], data);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherPrecondition<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2>& m,
  const MGLevelObject<DATA>& data)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &m[i];
      smoothers[i].initialize(m[i], data[i]);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherPrecondition<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2>& m,
  const DATA& data,
  const unsigned int row,
  const unsigned int col)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
template <class MATRIX2, class DATA>
inline void
MGSmootherPrecondition<MATRIX, RELAX, VECTOR>::initialize (
  const MGLevelObject<MATRIX2>& m,
  const MGLevelObject<DATA>& data,
  const unsigned int row,
  const unsigned int col)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data[i]);
    }
}

template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherPrecondition<MATRIX, RELAX, VECTOR>::smooth(
  const unsigned int level,
  VECTOR& u,
  const VECTOR& rhs) const
{
  unsigned int maxlevel = matrices.get_maxlevel();
  unsigned int steps2 = this->steps;

  if (this->variable)
    steps2 *= (1<<(maxlevel-level));

  VECTOR* r = this->mem.alloc();
  VECTOR* d = this->mem.alloc();
  r->reinit(u);
  d->reinit(u);

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
	  matrices[level].Tvmult(*r,u);
	  r->sadd(-1.,1.,rhs);
	  if (this->debug > 2)
	    deallog << ' ' << r->l2_norm() << ' ';
	  smoothers[level].Tvmult(*d, *r);
	  if (this->debug > 1)
	    deallog << ' ' << d->l2_norm() << ' ';
	} else {
	  if (this->debug > 0)
	    deallog << 'N';
	  matrices[level].vmult(*r,u);
	  r->sadd(-1.,1.,rhs);
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
  
  this->mem.free(r);
  this->mem.free(d);
}


template <class MATRIX, class RELAX, class VECTOR>
inline unsigned int
MGSmootherPrecondition<MATRIX, RELAX, VECTOR>::
memory_consumption () const
{
  return sizeof(*this)
    + matrices.memory_consumption()
    + smoothers.memory_consumption()
    + this->mem.memory_consumption();
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
