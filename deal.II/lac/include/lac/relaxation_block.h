//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__relaxation_block_h
#define __deal2__relaxation_block_h

#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/vector.h>
#include <lac/precondition_block_base.h>
#include <lac/block_list.h>

#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

/**
 * Base class for the implementation of overlapping, multiplicative
 * Schwarz relaxation methods and smoothers.
 *
 * This class uses the infrastructure provided by
 * PreconditionBlockBase. It adds functions to initialize with a block
 * list and to do the relaxation step. The actual relaxation method
 * with the interface expected by SolverRelaxation and
 * MGSmootherRelaxation is in the derived classes.
 *
 * This class allows for more general relaxation methods than
 * PreconditionBlock, since the index sets may be arbitrary and
 * overlapping, while there only contiguous, disjoint sets of equal
 * size are allowed. As a drawback, this class cannot be used as a
 * preconditioner, since its implementation relies on a straight
 * forward implementation of the Gauss-Seidel process.
 *
 * @author Guido Kanschat
 * @date 2010
 */
template <class MATRIX, typename inverse_type=typename MATRIX::value_type>
class RelaxationBlock :
  protected PreconditionBlockBase<inverse_type>
{
  private:
				     /**
				      * Define number type of matrix.
				      */
    typedef typename MATRIX::value_type number;

				     /**
				      * Value type for inverse matrices.
				      */
    typedef inverse_type value_type;

  public:
				     /**
				      * Parameters for block
				      * relaxation methods.
				      */
    class AdditionalData
    {
      public:
					 /**
					  * Default constructor
					  */
	AdditionalData ();
	
					 /**
					  * Constructor.
					  */
	AdditionalData (const BlockList& block_list,
			const double relaxation = 1.,
			const bool invert_diagonal = true,
			const bool same_diagonal = false);
	
					 /**
					  * Pointer to the DoFHandler
					  * providing the cells.
					  */
	SmartPointer<const BlockList, typename RelaxationBlock<MATRIX, inverse_type>::AdditionalData> block_list;
	
					 /**
					  * Relaxation parameter.
					  */
	double relaxation;

					 /**
					  * Invert diagonal during initialization.
					  */
	bool invert_diagonal;

					 /**
					  * Assume all diagonal blocks
					  * are equal to save memory.
					  */
	bool same_diagonal;
    };
    
				     /**
				      * Initialize matrix and block
				      * size.  We store the matrix and
				      * the block size in the
				      * preconditioner object. In a
				      * second step, the inverses of
				      * the diagonal blocks may be
				      * computed.
				      *
				      * Additionally, a relaxation
				      * parameter for derived classes
				      * may be provided.
				      */
    void initialize (const MATRIX& A,
		     const AdditionalData parameters);

				     /**
				      * Deletes the inverse diagonal
				      * block matrices if existent,
				      * sets the blocksize to 0, hence
				      * leaves the class in the state
				      * that it had directly after
				      * calling the constructor.
				      */
    void clear();

				     /**
				      * Checks whether the object is empty.
				      */
    bool empty () const;

				     /**
				      * Read-only access to entries.
				      * This function is only possible
				      * if the inverse diagonal blocks
				      * are stored.
				      */
    value_type el(unsigned int i,
		  unsigned int j) const;

    				     /**
				      * Stores the inverse of the
				      * diagonal blocks in
				      * @p inverse. This costs some
				      * additional memory - for DG
				      * methods about 1/3 (for double
				      * inverses) or 1/6 (for float
				      * inverses) of that used for the
				      * matrix - but it makes the
				      * preconditioning much faster.
				      *
				      * It is not allowed to call this
				      * function twice (will produce
				      * an error) before a call of
				      * <tt>clear(...)</tt>  because at the
				      * second time there already
				      * exist the inverse matrices.
				      *
				      * After this function is called,
				      * the lock on the matrix given
				      * through the @p use_matrix
				      * function is released, i.e. you
				      * may overwrite of delete it.
				      * You may want to do this in
				      * case you use this matrix to
				      * precondition another matrix.
				      */
    void invert_diagblocks();

  protected:
				     /**
				      * Perform one block relaxation
				      * step.
				      *
				      * Depending on the arguments @p
				      * dst and @p pref, this performs
				      * an SOR step (both reference
				      * the same vector) of a Jacobi
				      * step (both are different
				      * vectors). For the Jacobi step,
				      * the calling function must copy
				      * @p dst to @p pref after this.
				      */
    template <typename number2>
    void do_step (
      Vector<number2>       &dst,
      const Vector<number2> &prev,
      const Vector<number2> &src,
      const bool backward) const;
				     /**
				      * Pointer to the matrix. Make
				      * sure that the matrix exists as
				      * long as this class needs it,
				      * i.e. until calling
				      * @p invert_diagblocks, or (if
				      * the inverse matrices should
				      * not be stored) until the last
				      * call of the preconditoining
				      * @p vmult function of the
				      * derived classes.
				      */
    SmartPointer<const MATRIX,RelaxationBlock<MATRIX,inverse_type> > A;

				     /**
				      * Control information.
				      */
    AdditionalData additional_data;
};


/**
 * Block Gauss-Seidel method with possibly overlapping blocks.
 *
 * This class implements the step() and Tstep() functions expected by
 * SolverRelaxation and MGSmootherRelaxation. They perform a
 * multiplicative Schwarz method on the blocks provided in the
 * BlockList of AdditionalData. Differing from PreconditionBlockSOR,
 * these blocks may be of varying size, non-contiguous, and
 * overlapping. On the other hand, this class does not implement the
 * preconditioner interface expected by Solver objects.
 *
 * @author Guido Kanschat
 * @data 2010
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class RelaxationBlockSOR : public virtual Subscriptor,
			     protected RelaxationBlock<MATRIX, inverse_type>
{
  public:
				     /**
				      * Default constructor.
				      */
//    RelaxationBlockSOR();
    
				     /**
				      * Define number type of matrix.
				      */
    typedef typename MATRIX::value_type number;

				     /**
				      * Make type publicly available.
				      */
    using RelaxationBlock<MATRIX,inverse_type>::AdditionalData;

				     /**
				      * Make initialization function
				      * publicly available.
				      */
    using RelaxationBlock<MATRIX, inverse_type>::initialize;

				     /**
				      * Make function of base class public again.
				      */
    using RelaxationBlock<MATRIX, inverse_type>::clear;

				     /**
				      * Make function of base class public again.
				      */
    using RelaxationBlock<MATRIX, inverse_type>::empty;
				     /**
				      * Perform one step of the SOR
				      * iteration.
				      */
    template <typename number2>
    void step (Vector<number2>& dst, const Vector<number2>& rhs) const;

				     /**
				      * Perform one step of the
				      * transposed SOR iteration.
				      */
    template <typename number2>
    void Tstep (Vector<number2>& dst, const Vector<number2>& rhs) const;

  protected:
				     /**
				      * Constructor to be used by
				      * RelaxationBlockSSOR.
				      */
//    RelaxationBlockSOR(bool store);
};


/**
 * Symmetric block Gauss-Seidel method with possibly overlapping blocks.
 *
 * This class implements the step() and Tstep() functions expected by
 * SolverRelaxation and MGSmootherRelaxation. They perform a
 * multiplicative Schwarz method on the blocks provided in the
 * BlockList of AdditionalData in symmetric fashion. Differing from
 * PreconditionBlockSSOR, these blocks may be of varying size,
 * non-contiguous, and overlapping. On the other hand, this class does
 * not implement the preconditioner interface expected by Solver
 * objects.
 *
 * @author Guido Kanschat
 * @data 2010
 */
template<class MATRIX, typename inverse_type = typename MATRIX::value_type>
class RelaxationBlockSSOR : public virtual Subscriptor,
			    protected RelaxationBlock<MATRIX, inverse_type>
{
  public:    
				     /**
				      * Define number type of matrix.
				      */
    typedef typename MATRIX::value_type number;

				     /**
				      * Make type publicly available.
				      */
    using RelaxationBlock<MATRIX,inverse_type>::AdditionalData;

				     /**
				      * Make initialization function
				      * publicly available.
				      */
    using RelaxationBlock<MATRIX, inverse_type>::initialize;

				     /**
				      * Make function of base class public again.
				      */
    using RelaxationBlock<MATRIX, inverse_type>::clear;

				     /**
				      * Make function of base class public again.
				      */
    using RelaxationBlock<MATRIX, inverse_type>::empty;
				     /**
				      * Perform one step of the SOR
				      * iteration.
				      */
    template <typename number2>
    void step (Vector<number2>& dst, const Vector<number2>& rhs) const;

				     /**
				      * Perform one step of the
				      * transposed SOR iteration.
				      */
    template <typename number2>
    void Tstep (Vector<number2>& dst, const Vector<number2>& rhs) const;
};


DEAL_II_NAMESPACE_CLOSE

#endif
