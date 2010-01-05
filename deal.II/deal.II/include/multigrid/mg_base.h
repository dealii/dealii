//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2010 by the deal.II authors
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


DEAL_II_NAMESPACE_OPEN


template <int, int> class MGDoFHandler;
template <typename> class MGLevelObject;


namespace internal
{
  namespace mg
  {
				     /**
				      * Ajust vectors on all levels to
				      * correct size.  Here, we just
				      * count the numbers of degrees
				      * of freedom on each level and
				      * @p reinit each level vector
				      * to this length.
				      */
    template <int dim, typename number, int spacedim>
    static void
    reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
		   MGLevelObject<dealii::Vector<number> > &vector);

				     /**
				      * Ajust vectors on all levels to
				      * correct size.  Here, we just
				      * count the numbers of degrees
				      * of freedom on each level and
				      * @p reinit each level vector
				      * to this length.
				      */
    template <int dim, typename number, int spacedim>
    static void
    reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
		   MGLevelObject<BlockVector<number> > &vector);


				     /**
				      * Adjust vectors on all levels
				      * to correct size. The degrees
				      * of freedom on each level are
				      * counted by block and only the
				      * block selected is used.
				      */
    template <int dim, typename number, int spacedim>
    static void
    reinit_vector_by_blocks (
      const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
      MGLevelObject<dealii::Vector<number> > &v,
      const unsigned int selected,
      std::vector<std::vector<unsigned int> >& cached_sizes);

				     /**
				      * Adjust block vectors on all
				      * levels to correct size. The
				      * degrees of freedom on each
				      * level are counted by block.
				      */
    template <int dim, typename number, int spacedim>
    static void
    reinit_vector_by_blocks (
      const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
      MGLevelObject<BlockVector<number> > &v,
      const std::vector<bool> &selected,
      std::vector<std::vector<unsigned int> >& cached_sizes);

				     /**
				      * Adjust block-vectors on all
				      * levels to correct size.  Count
				      * the numbers of degrees of
				      * freedom on each level
				      * component-wise. Then, assign
				      * each block of @p vector the
				      * corresponding size.
				      *
				      * The boolean field @p selected
				      * allows restricting this
				      * operation to certain
				      * components. In this case, @p
				      * vector will only have as many
				      * blocks as there are true
				      * values in @p selected (no
				      * blocks of length zero are
				      * padded in). If this argument
				      * is omitted, all blocks will be
				      * considered.
				      *
				      * Degrees of freedom must be
				      * sorted by component in order
				      * to obtain reasonable results
				      * from this function.
				      *
				      * The argument
				      * @p target_component allows to
				      * re-sort and group components
				      * as in
				      * DoFRenumbering::component_wise.
				      *
				      *
				      */
    template <int dim, typename number, int spacedim>
    static void
    reinit_vector_by_components (const dealii::MGDoFHandler<dim,spacedim>& mg_dof,
				 MGLevelObject<BlockVector<number> >& v,
				 const std::vector<bool>& selected,
				 const std::vector<unsigned int>& target_component,
				 std::vector<std::vector<unsigned int> >& cached_sizes);
				     /**
				      * Adjust vectors on all levels
				      * to correct size.  Count the
				      * numbers of degrees of freedom
				      * on each level component-wise
				      * in a single component. Then,
				      * assign @p vector the
				      * corresponding size.
				      *
				      * The boolean field @p selected
				      * may be nonzero in a single
				      * component, indicating the
				      * block of a block vector the
				      * argument @p v corresponds to.
				      *
				      * Degrees of freedom must be
				      * sorted by component in order
				      * to obtain reasonable results
				      * from this function.
				      *
				      * The argument
				      * @p target_component allows to
				      * re-sort and group components
				      * as in
				      * DoFRenumbering::component_wise.
				      */
    template <int dim, typename number, int spacedim>
    static void
    reinit_vector_by_components (
      const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
      MGLevelObject<dealii::Vector<number> > &v,
      const std::vector<bool> &selected,
      const std::vector<unsigned int> &target_component,
      std::vector<std::vector<unsigned int> >& cached_sizes);

  }
}



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
  virtual void vmult (const unsigned int level,
		      VECTOR& dst,
		      const VECTOR& src) const = 0;

				   /**
				    * Adding matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void vmult_add (const unsigned int level,
			  VECTOR& dst,
			  const VECTOR& src) const = 0;

				   /**
				    * Transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult (const unsigned int level,
		       VECTOR& dst,
		       const VECTOR& src) const = 0;

				   /**
				    * Adding transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult_add (const unsigned int level,
			   VECTOR& dst,
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
 * If you use multigrid for a single PDF or for your complete system
 * of equations, you
 * will use MGTransferPrebuilt together with Multigrid. The vector
 * types used on the fine grid as well as for the multilevel
 * operations may be Vector or BlockVector. In both cases,
 * MGTransferPrebuilt will operate on all components of the solution.
 *
 * @note For the following, it is important to realize the difference
 * between a solution @ref GlossComponent "component" and a solution
 * @ref GlossBlock "block". The distinction only applies if vector valued
 * elements are used, but is quite important then. This is reflected
 * in the fact that it is not possible right now to use transfer
 * classes based on MGTransferComponentBase for genuine vector valued
 * elements, but descendants of MGTransferBlockBase would have to be
 * applied. In the following text, we will use the term
 * <em>block</em>, but remark that it might refer to components as well.
 *
 * @todo update the following documentation, since it does not reflect
 * the latest changes in structure.
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
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2002, 2007
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
				      * Release matrices.
				      */
    virtual void clear() = 0;

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

DEAL_II_NAMESPACE_CLOSE

#endif
