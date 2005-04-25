//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_block_smoother_h
#define __deal2__mg_block_smoother_h


#include <base/config.h>
#include <base/smartpointer.h>
#include <lac/pointer_matrix.h>
#include <lac/vector_memory.h>
#include <lac/block_vector.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_level_object.h>
#include <vector>

template <int dim> class MGDoFHandler;


/*
 * MGSmootherBase is defined in mg_base.h
 */

/*!@addtogroup mg */
/*@{*/

/**
 * General smoother class for block vectors. This class gives complete
 * freedom to the choice of a block smoother by being initialized with
 * a matrix and a smoother object. Therefore, the smoother object for
 * each level must be constructed by hand.
 *
 * @author Guido Kanschat, 2005
 */
template <class MATRIX, class RELAX, typename number>
class MGSmootherBlock
  : public MGSmootherBase<BlockVector<number> >
{
  public:
    				     /**
				      * Constructor. Sets memory and
				      * smoothing parameters.
				      */
    MGSmootherBlock(VectorMemory<BlockVector<number> >& mem,
		    const unsigned int steps = 1,
		    const bool variable = false,
		    const bool symmetric = false,
		    const bool transpose = false,
		    const bool reverse = false);

				     /**
				      * Initialize for matrices. The
				      * parameter <tt>matrices</tt> can be
				      * any object having functions
				      * <tt>get_minlevel()</tt> and
				      * <tt>get_maxlevel()</tt> as well as
				      * an <tt>operator[]</tt> returning a
				      * reference to @p MATRIX.
				      *
				      * The same convention is used
				      * for the parameter
				      * <tt>smoothers</tt>, such that
				      * <tt>operator[]</tt> returns
				      * the object doing the
				      * block-smoothing on a single
				      * level.
				      *
				      * This function stores pointers
				      * to the level matrices and
				      * smoothing operator for each
				      * level.
				      */
    template <class MGMATRIX, class MGRELAX>
    void initialize (const MGMATRIX& matrices,
		     const MGRELAX& smoothers);

				     /**
				      * Empty all vectors.
				      */
    void clear ();
    
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
				      * Switch on/off transposed. This
				      * is mutually exclusive with
				      * reverse().
				      */
    void set_transpose (const bool);
    
				     /**
				      * Switch on/off reversed. This
				      * is mutually exclusive with
				      * transpose().
				      */
    void set_reverse (const bool);

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
    virtual void smooth (const unsigned int         level,
			 BlockVector<number>&       u,
			 const BlockVector<number>& rhs) const;
 private:
				     /**
				      * Pointer to the matrices.
				      */
    MGLevelObject<PointerMatrix<MATRIX, BlockVector<number> > > matrices;

				     /**
				      * Pointer to the matrices.
				      */
    MGLevelObject<PointerMatrix<RELAX, BlockVector<number> > > smoothers;

				     /**
				      * Number of smoothing steps.
				      */
    unsigned int steps;

				     /**
				      * Variable smoothing?
				      */
    bool variable;

				     /**
				      * Symmetric smoothing?
				      */
    bool symmetric;

				     /*
				      * Transposed?
				      */
    bool transpose;

				     /**
				      * Reverse?
				      */
    bool reverse;

				     /**
				      * Memory for auxiliary vectors.
				      */
    VectorMemory<BlockVector<number> >& mem;

};


//---------------------------------------------------------------------------

template <class MATRIX, class RELAX, typename number>
inline
MGSmootherBlock<MATRIX, RELAX, number>::MGSmootherBlock(
  VectorMemory<BlockVector<number> >& mem,
  const unsigned int steps,
  const bool variable,
  const bool symmetric,
  const bool transpose,
  const bool reverse)
		:
		steps(steps),
		variable(variable),
		symmetric(symmetric),
		transpose(transpose),
		reverse(reverse),
		mem(mem)
{}


template <class MATRIX, class RELAX, typename number>
inline void
MGSmootherBlock<MATRIX, RELAX, number>::clear ()
{  
  unsigned int i=matrices.get_minlevel(),
       max_level=matrices.get_maxlevel();
  for (; i<=max_level; ++i)
    {
      smoothers[i] = 0;
      matrices[i] = 0;
    }
}


template <class MATRIX, class RELAX, typename number>
template <class MGMATRIX, class MGRELAX>
inline void
MGSmootherBlock<MATRIX, RELAX, number>::initialize (
  const MGMATRIX& m,
  const MGRELAX& s)
{
  const unsigned int min = m.get_minlevel();
  const unsigned int max = m.get_maxlevel();
  
  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i=min;i<=max;++i)
    {
      matrices[i] = &m[i];
      smoothers[i] = &s[i];
    }
}

template <class MATRIX, class RELAX, typename number>
inline void
MGSmootherBlock<MATRIX, RELAX, number>::
set_steps (const unsigned int s)
{
  steps = s;
}


template <class MATRIX, class RELAX, typename number>
inline void
MGSmootherBlock<MATRIX, RELAX, number>::
set_variable (const bool flag)
{
  variable = flag;
}


template <class MATRIX, class RELAX, typename number>
inline void
MGSmootherBlock<MATRIX, RELAX, number>::
set_symmetric (const bool flag)
{
  symmetric = flag;
}


template <class MATRIX, class RELAX, typename number>
inline void
MGSmootherBlock<MATRIX, RELAX, number>::
set_transpose (const bool flag)
{
  transpose = flag;
}


template <class MATRIX, class RELAX, typename number>
inline void
MGSmootherBlock<MATRIX, RELAX, number>::
set_reverse (const bool flag)
{
  reverse = flag;
}


template <class MATRIX, class RELAX, typename number>
inline void
MGSmootherBlock<MATRIX, RELAX, number>::smooth(
  const unsigned int level,
  BlockVector<number>& u,
  const BlockVector<number>& rhs) const
{
  unsigned int maxlevel = matrices.get_maxlevel();
  unsigned int steps2 = steps;

  if (variable)
    steps2 *= (1<<(maxlevel-level));

  BlockVector<number>* r = mem.alloc();
  BlockVector<number>* d = mem.alloc();
  r->reinit(u);
  d->reinit(u);

  bool T = transpose;
  if (symmetric && (steps2 % 2 == 0))
    T = false;
//  cerr << 'S' << level;
//  cerr << '(' << matrices[level]->m() << ',' << matrices[level]->n() << ')';
  
  for (unsigned int i=0; i<steps2; ++i)
    {
      if (T)
	{
//	  cerr << 'T';
	  matrices[level].vmult(*r,u);
	  r->sadd(-1.,1.,rhs);
	  smoothers[level].Tvmult(*d, *r);
	} else {
//	  cerr << 'N';
	  matrices[level].vmult(*r,u);
	  r->sadd(-1.,1.,rhs);
	  smoothers[level].vmult(*d, *r);
	}
//      cerr << '{' << r->l2_norm() << '}';
      u += *d;
      if (symmetric)
	T = !T;
 }

  mem.free(r);
  mem.free(d);
}



#endif
