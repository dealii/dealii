//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_base_h
#define __deal2__mg_base_h

/*
 * This file contains MGLevelObject and some abstract base classes
 * used by Multigrid.
 */

#include <base/config.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/vector.h>


/*!@addtogroup mg */
/*@{*/


/**
 * Multilevel matrix base. This class sets up the interface needed by
 * multilevel algorithms. It has no relation to the actual matrix type
 * and takes the vector class as only template argument.
 *
 * Usually, the derived class MGMatrix, operating on an
 * MGLevelObject of matrices will be sufficient for applications.
 *
 * @author Guido Kanschat, 2002
 */
template <class VECTOR>
class MGMatrixBase : public Subscriptor
{
  public:
				   /*
				    * Virtual destructor.
				    */
  virtual ~MGMatrixBase();

				   /**
				    * Matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void vmult(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const = 0;

				   /**
				    * Adding matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void vmult_add(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const = 0;

				   /**
				    * Transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const = 0;

				   /**
				    * Adding transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult_add(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const = 0;
};


/**
 * Base class for coarse grid solvers.  This defines the virtual
 * parenthesis operator, being the interface used by multigrid
 * methods. Any implementation will be done by derived classes.
 *
 * @author Guido Kanschat, 2002
 */
template <class VECTOR>
class MGCoarseGridBase : public Subscriptor
{
  public:
				     /**
				      * Virtual destructor.
				      */
    virtual ~MGCoarseGridBase ();

				     /**
				      * Solver method implemented by
				      * derived classes.
				      */
    virtual void operator() (const unsigned int   level,
			     VECTOR       &dst,
			     const VECTOR &src) const = 0;
};


/**
 * Base class used to declare the operations needed by a concrete class
 * implementing prolongation and restriction of vectors in the multigrid
 * context. This class is abstract and has no implementation of
 * these operations.
 *
 * There are several derived classes, reflecting the fact that vector
 * types and numbering of the fine-grid discretization and of the
 * multi-level implementation are independent.
 *
 * If you use multigrid for your complete system of equations, you
 * will use MGTransferPrebuilt together with Multigrid. The vector
 * types used on the fine grid as well as for the multilevel
 * operations may be Vector or BlockVector. In both cases,
 * MGTransferPrebuilt will operate on all components of the solution.
 *
 * For mixed systems, it may be required to do multigrid only for a
 * single component or for some components. The classes
 * MGTransferSelect and MGTransferBlock handle these cases.
 *
 * MGTransferSelect is used if you use mutligrid (on Vector objects)
 * for a single component, possibly grouped using
 * <tt>mg_target_component</tt>.
 *
 * The class MGTransferBlock handles the case where your multigrid
 * method operates on BlockVector objects. These can contain all or a
 * consecutive set of the blocks of the complete system. Since most
 * smoothers cannot operate on block structures, it is not clear
 * whether this case is really useful. Therefore, a tested
 * implementation of this case will be supplied when needed.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2002
 */
 template <class VECTOR>
 class MGTransferBase : public Subscriptor
 {
   public:
 				     /**
 				      * Destructor. Does nothing here, but
 				      * needs to be declared virtual anyway.
 				      */
     virtual ~MGTransferBase();

 				     /**
 				      * Prolongate a vector from level
 				      * <tt>to_level-1</tt> to level
 				      * <tt>to_level</tt>. The previous
 				      * content of <tt>dst</tt> is
 				      * overwritten.
 				      *
 				      * @arg src is a vector with as
 				      * many elements as there are
 				      * degrees of freedom on the
 				      * coarser level involved.
				      *
 				      * @arg dst has as many elements
 				      * as there are degrees of
 				      * freedom on the finer level.
 				      */
     virtual void prolongate (const unsigned int to_level,
 			     VECTOR&            dst,
 			     const VECTOR&      src) const = 0;

 				     /**
 				      * Restrict a vector from level
 				      * <tt>from_level</tt> to level
 				      * <tt>from_level-1</tt> and add
 				      * this restriction to
 				      * <tt>dst</tt>. If the region
 				      * covered by cells on level
 				      * <tt>from_level</tt> is smaller
 				      * than that of level
 				      * <tt>from_level-1</tt> (local
 				      * refinement), then some degrees
 				      * of freedom in <tt>dst</tt> are
 				      * active and will not be
 				      * altered. For the other degress
 				      * of freedom, the result of the
 				      * restriction is added.
 				      *
 				      * @arg src is a vector with as
 				      * many elements as there are
 				      * degrees of freedom on the
 				      * finer level
				      *
				      * @arg dst has as many elements as there
 				      * are degrees of freedom on the coarser
 				      * level.
				      *
 				      */
     virtual void restrict_and_add (const unsigned int from_level,
 				   VECTOR&            dst,
 				   const VECTOR&      src) const = 0;
 };



/**
 * Base class for multigrid smoothers. Does nothing but defining the
 * interface used by multigrid methods.
 *
 * @author Guido Kanschat, 2002
 */
template <class VECTOR>
class MGSmootherBase : public Subscriptor
{
  public:
				   /**
				    * Virtual destructor.
				    */
  virtual ~MGSmootherBase();

				   /**
				    * Smoothing function. This is the
				    * function used in multigrid
				    * methods.
				    */
  virtual void smooth (const unsigned int level,
		       VECTOR&            u,
		       const VECTOR&      rhs) const = 0;  
};

/*@}*/


#endif
