//----------------------------  mg_smoother.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_smoother.h  ---------------------------
#ifndef __deal2__mg_smoother_h
#define __deal2__mg_smoother_h


#include <base/config.h>
#include <base/smartpointer.h>
#include <lac/pointer_matrix.h>
#include <multigrid/mg_base.h>
#include <vector>

template <int dim> class MGDoFHandler;


/**
 * Smoother doing nothing. This class is not useful for many applications other
 * than for testing some multigrid procedures. Also some applications might
 * get convergence without smoothing and then this class brings you the
 * cheapest possible multigrid.
 *
 * @author Guido Kanschat, 1999, 2002
 */
template <class VECTOR>
class MGSmootherIdentity : public MGSmoother<VECTOR>
{
  public:
				     /**
				      * Implementation of the
				      * interface for @p{Multigrid}.
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
 * Base class for multigrid smoothers. It implements some
 * functionality for setting the values of vectors at interior
 * boundaries (i.e. boundaries between differing levels of the
 * triangulation) to zero, by building a list of these degrees of
 * freedom's indices at construction time.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2002
 */
template <class VECTOR>
class MGSmootherContinuous : public MGSmoother<VECTOR>
{
  private:
				     /**
				      * Default constructor. Made private
				      * to prevent it being called, which
				      * is necessary since this could
				      * circumvent the set-up of the list
				      * if interior boundary dofs.
				      */
    MGSmootherContinuous ();
    
  public:

				     /**
				      * Constructor. This one collects
				      * the indices of the degrees of freedom
				      * on the interior boundaries between
				      * the different levels, which are
				      * needed by the function
				      * @p{set_zero_interior_boundaries}.
				      *
				      * Since this function is
				      * implemented a bit different in
				      * 1d (there are no faces of
				      * cells, just vertices), there
				      * are actually two sets of
				      * constructors, namely this one
				      * for 1d and the following one
				      * for all other
				      * dimensions. Really amusing
				      * about this text is, that there
				      * is no 1d implementation.
				      */
    MGSmootherContinuous (const MGDoFHandler<1> &mg_dof,
			  const unsigned int     steps);

				     /**
				      * Constructor. This one collects
				      * the indices of the degrees of
				      * freedom on the interior
				      * boundaries between the
				      * different levels, which are
				      * needed by the function
				      * @p{set_zero_interior_boundaries}.
				      *
				      * The parameter steps indicates
				      * the number of smoothing steps
				      * to be executed by @p{smooth}.
				      */
    template <int dim>
    MGSmootherContinuous (const MGDoFHandler<dim> &mg_dof,
			  const unsigned int       steps);    

				     /**
				      * Reset the values of the degrees of
				      * freedom on interior boundaries between
				      * different levels to zero in the given
				      * data vector @p{u}.
				      *
				      * Since the coarsest level (@p{level==0})
				      * has no interior boundaries, this
				      * function does nothing in this case.
				      */
    void set_zero_interior_boundary (const unsigned int level,
				     VECTOR&            u) const;
    
				     /**
				      * Modify the number of smoothing steps.
				      */
    void set_steps(unsigned int steps);
    
				     /**
				      * How many steps should be used?
				      */
    unsigned int get_steps() const;

  private:
				     /**
				      * Number of smoothing steps.
				      */
    unsigned int steps;
    
				     /**
				      * For each level, we store a list of
				      * degree of freedom indices which are
				      * located on interior boundaries between
				      * differing levels of the triangulation.
				      * Since the coarsest level has no
				      * interior boundary dofs, the first
				      * entry refers to the second level.
				      *
				      * These arrays are set by the constructor.
				      * The entries for each level are sorted ascendingly.
				      */
    std::vector<std::vector<unsigned int> > interior_boundary_dofs;
};


/**
 * Smoother using relaxation classes.
 *
 * This class performs smoothing on each level. The operation can be
 * controlled by several parameters. First, the relaxation parameter
 * @p{omega} is used in the underlying relaxation method. @p{steps} is
 * the number of relaxation steps on the finest level (on all levels
 * if @p{variable} is off). If @p{variable} is @p{true}, the number of
 * smoothing steps is doubled on each coarser level. This results in a
 * method having the complexity of the W-cycle, but saving grid
 * transfers. This is the method proposed by Bramble at al.
 *
 * The option @p{symmetric} switches on alternating between the
 * smoother and its transpose in each step as proposed by Bramble.
 *
 * @p{transpose} and @p{reverse} are similar in the effect that
 * instead of the smoother the transposed is used. Typically, this is
 * off for pre-smoothing and on for post-smoothing. While
 * @p{transpose} is the true transposed smoothing operation,
 * @p{reverse} just uses the transposed of the smoother, but the
 * non-transposed matrix-vector multiplication; this can be used to
 * invert the direction of the Gauss-Seidel method.
 *
 * If you are using block matrices, the second @p{initialize} function
 * offers the possibility to extract a single block for smoothing. In
 * this case, the multigrid method must be used only with the vector
 * associated to that single block.
 *
 * The library contains instantiation for @p{SparseMatrix<.>} and
 * @p{Vector<.>}, where the template arguments are all combinations of
 * @p{float} and @p{double}. Additional instantiations may be created
 * by including the file mg_smoother.templates.h.
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
			 const bool transpose = false,
			 const bool reverse = false);

				     /**
				      * Initialize for matrices. The
				      * parameter @p{matrices} can be
				      * any object having functions
				      * @p{get_minlevel()} and
				      * @p{get_maxlevel()} as well as
				      * an @p{operator[]} returning a
				      * reference to @p{MATRIX}. This
				      * function stores pointers to
				      * the level matrices and
				      * initializes the smoothing
				      * operator for each level.
				      *
				      * @p{additional_data} is an
				      * object of type
				      * @p{RELAX::AdditionalData} and
				      * is handed to the
				      * initialization function of the
				      * relaxation method.
				      */
    template <class MGMATRIXOBJECT, class DATA>
    void initialize (const MGMATRIXOBJECT& matrices,
		     const DATA& additional_data);

				     /**
				      * Initialize for single blocks
				      * of matrices. The parameter
				      * @p{matrices} can be any object
				      * having functions
				      * @p{get_minlevel()} and
				      * @p{get_maxlevel()} as well as
				      * an @p{operator[]} returning a
				      * reference to a block matrix
				      * where each block is of type
				      * @p{MATRIX}. Of this block
				      * matrix, the block indicated by
				      * @p{block_row} and
				      * @p{block_col} is selected on
				      * each level.  This function
				      * stores pointers to the level
				      * matrices and initializes the
				      * smoothing operator for each
				      * level.
				      *
				      * @p{additional_data} is an
				      * object of type
				      * @p{RELAX::AdditionalData} and
				      * is handed to the
				      * initialization function of the
				      * relaxation method.
				      */
    template <class MGMATRIXOBJECT, class DATA>
    void initialize (const MGMATRIXOBJECT& matrices,
		     const DATA& additional_data,
		     const unsigned int block_row,
		     const unsigned int block_col);

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
				      * @ref{reverse}.
				      */
    void set_transpose (const bool);
    
				     /**
				      * Switch on/off reversed. This
				      * is mutually exclusive with
				      * @ref{transpose}.
				      */
    void set_reverse (const bool);

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
    
 private:
				     /**
				      * Pointer to the matrices.
				      */
    MGLevelObject<PointerMatrix<MATRIX, VECTOR> > matrices;

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
    VectorMemory<VECTOR>& mem;
};


/* ------------------------------- Inline functions -------------------------- */

template <class VECTOR>
inline void
MGSmootherIdentity<VECTOR>::smooth (
  const unsigned int, VECTOR&,
  const VECTOR&) const
{}

/*----------------------------------------------------------------------*/

template <class VECTOR>
inline
void
MGSmootherContinuous<VECTOR>::set_steps(unsigned int i)
{
  steps = i;
}

template <class VECTOR>
inline
unsigned int
MGSmootherContinuous<VECTOR>::get_steps() const
{
  return steps;
}

//----------------------------------------------------------------------//

template <class MATRIX, class RELAX, class VECTOR>
inline
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::MGSmootherRelaxation(
  VectorMemory<VECTOR>& mem,
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
template <class MGMATRIXOBJECT, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGMATRIXOBJECT& m,
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
template <class MGMATRIXOBJECT, class DATA>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::initialize (
  const MGMATRIXOBJECT& m,
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
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::
set_steps (const unsigned int s)
{
  steps = s;
}


template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::
set_variable (const bool flag)
{
  variable = flag;
}


template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::
set_symmetric (const bool flag)
{
  symmetric = flag;
}


template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::
set_transpose (const bool flag)
{
  transpose = flag;
}


template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::
set_reverse (const bool flag)
{
  reverse = flag;
}


template <class MATRIX, class RELAX, class VECTOR>
inline void
MGSmootherRelaxation<MATRIX, RELAX, VECTOR>::smooth(
  const unsigned int level,
  VECTOR& u,
  const VECTOR& rhs) const
{
  unsigned int maxlevel = matrices.get_maxlevel();
  unsigned int steps2 = steps;

  if (variable)
    steps2 *= (1<<(maxlevel-level));

  Vector<double>* r = mem.alloc();
  Vector<double>* d = mem.alloc();
  r->reinit(u);
  d->reinit(u);

  bool T = transpose;
  if (symmetric && (steps2 % 2 == 0))
    T = false;
//  cerr << 'S' << level;
  
  for (unsigned int i=0; i<steps2; ++i)
    {
      if (T)
	{
					   // This is not really the
					   // transposed smoother, but
					   // just Gauss-Seidel with
					   // reverse numbering.  For
					   // a symmetric matrix, it
					   // is the transpose,
					   // though.
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
