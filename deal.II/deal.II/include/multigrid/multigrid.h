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
#ifndef __deal2__multigrid_h
#define __deal2__multigrid_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <multigrid/mg_base.h>
#include <multigrid/mg_level_object.h>
#include <multigrid/mg_dof_handler.h>

#include <vector>

/*!@addtogroup mg */
/*@{*/

/**
 * Implementation of the multigrid method.
 *
 * @warning multigrid on locally refined meshes only works with
 * <b>discontinuous finite elements</b> right now. It is not clear,
 * whether the paradigm of local smoothing we use is applicable to
 * continuous elements with hanging nodes; in fact, most people you
 * meet on conferences seem to deny this.
 *
 * The function which starts a multigrid cycle on the finest level is
 * cycle(). Depending on the cycle type chosen with the constructor
 * (see enum Cycle), this function triggers one of the cycles
 * level_v_step() or level_step(), where the latter one can do
 * different types of cycles.
 *
 * Using this class, it is expected that the right hand side has been
 * converted from a vector living on the locally finest level to a
 * multilevel vector. This is a nontrivial operation, usually
 * initiated automatically by the class PreconditionMG and performed
 * by the classes derived from MGTransferBase.
 *
 * @note The interface of this class is still very clumsy. In
 * particular, you will have to set up quite a few auxiliary objects
 * before you can use it. Unfortunately, it seems that this can be
 * avoided only be restricting the flexibility of this class in an
 * unacceptable way.
 *
 * @author Guido Kanschat, 1999 - 2005
 */
template <class VECTOR>
class Multigrid : public Subscriptor
{
  public:
				     /**
				      * List of implemented cycle types.
				      */
    enum Cycle
    {
					   /// The V-cycle
	  v_cycle,
					   /// The W-cycle
	  w_cycle,
					   /// The F-cycle
	  f_cycle
    };
    
    typedef VECTOR vector_type;
    typedef const VECTOR const_vector_type;
  
				     /**
				      * Constructor. The
				      * MGDoFHandler is used to
				      * determine the highest possible
				      * level. <tt>transfer</tt> is an
				      * object performing prolongation
				      * and restriction.
				      *
				      * This function already
				      * initializes the vectors which
				      * will be used later in the
				      * course of the
				      * computations. You should
				      * therefore create objects of
				      * this type as late as possible.
				      */
    template <int dim>
    Multigrid(const MGDoFHandler<dim>& mg_dof_handler,
	      const MGMatrixBase<VECTOR>& matrix,
	      const MGCoarseGridBase<VECTOR>& coarse,
	      const MGTransferBase<VECTOR>& transfer,
	      const MGSmootherBase<VECTOR>& pre_smooth,
	      const MGSmootherBase<VECTOR>& post_smooth,
	      Cycle cycle = v_cycle);
				     /**
				      * Experimental constructor for
				      * cases in which no MGDoFHandler
				      * is available.
				      *
				      * @warning Not intended for general use.
				      */
    Multigrid(const unsigned int minlevel,
	      const unsigned int maxlevel,
	      const MGMatrixBase<VECTOR>& matrix,
	      const MGCoarseGridBase<VECTOR>& coarse,
	      const MGTransferBase<VECTOR>& transfer,
	      const MGSmootherBase<VECTOR>& pre_smooth,
	      const MGSmootherBase<VECTOR>& post_smooth,
	      Cycle cycle = v_cycle);

				     /**
				      * Reinit this class according to
				      * #minlevel and #maxlevel.
				      */
    void reinit (const unsigned int minlevel,
		 const unsigned int maxlevel);

				     /**
				      * Execute one multigrid
				      * cycle. The type of cycle is
				      * selected by the constructor
				      * argument cycle. See the enum
				      * Cycle for available types.
				      */
    void cycle ();
    
				     /**
				      * Execute one step of the
				      * V-cycle algorithm.  This
				      * function assumes, that the
				      * multilevel vector #defect is
				      * filled with the residual of an
				      * outer defect correction
				      * scheme. This is usually taken
				      * care of by
				      * PreconditionMG). After
				      * vcycle(), the result is in the
				      * multilevel vector
				      * #solution. See
				      * <tt>copy_*_mg</tt> in class
				      * MGTools if you want to use
				      * these vectors yourself.
				      *
				      * The actual work for this
				      * function is done in
				      * level_v_step().
				      */
    void vcycle ();

				     /**
				      * @deprecated This function is
				      * purely experimental and will
				      * probably never be implemented
				      * in a way that it can be
				      * released.
				      *
				      * Perform a multigrid cycle with
				      * a vector which is already a
				      * level vector. Use of this
				      * function assumes that there is
				      * NO local refinement and that
				      * both vectors are on the finest
				      * level of this Multigrid
				      * object.
				      */
    void vmult(VECTOR& dst, const VECTOR& src) const;
    
				     /**
				      * @deprecated Even worse than
				      * vmult(), this function is not
				      * even implemented, but just
				      * declared such that certain
				      * objects relying on it can be
				      * constructed.
				      */
    void Tvmult(VECTOR& dst, const VECTOR& src) const;
    
				     /**
				      * Set additional matrices to
				      * correct residual computation
				      * at refinement edges. These
				      * matrices originate from
				      * discontinuous Galerkin methods
				      * (see FE_DGQ etc.), where they
				      * correspond tu the edge fluxes
				      * at the refinement edge between
				      * two levels.
				      */
    void set_edge_matrices (const MGMatrixBase<VECTOR>& edge_down,
			    const MGMatrixBase<VECTOR>& edge_up);

				     /**
				      * Return the finest level for
				      * multigrid.
				      */
    unsigned int get_maxlevel() const;
    
				     /**
				      * Return the coarsest level for
				      * multigrid.
				      */
    unsigned int get_minlevel() const;
    
				     /**
				      * Set the highest level for
				      * which the multilevel method is
				      * performed. By default, this is
				      * the finest level of the
				      * Triangulation; therefore, this
				      * function will only accept
				      * arguments smaller than the
				      * current #maxlevel and not
				      * smaller than the current
				      * #minlevel.
				      */
    void set_maxlevel(unsigned int);

				     /**
				      * Set the coarse level for which
				      * the multilevel method is
				      * performed. By default, this is
				      * zero. Accepted are
				      * non-negative values not larger than
				      * than the current #maxlevel.
				      *
				      * If <tt>relative</tt> ist
				      * <tt>true</tt>, then this
				      * function determins the number
				      * of levels used, that is, it
				      * sets #minlevel to
				      * #maxlevel-<tt>level</tt>.
				      *
				      * @note The mesh on the coarsest
				      * level must cover the whole
				      * domain. There may not be
				      * hanging nodes on #minlevel.
				      *
				      * @note If #minlevel is set to a
				      * nonzero value, do not forget
				      * to adjust your coarse grid
				      * solver!
				      */
    void set_minlevel(unsigned int level,
		      bool relative = false);

				     /**
				      * Chance #cycle_type used in cycle().
				      */
    void set_cycle(Cycle);
    
				     /**
				      * Set the debug level. Higher
				      * values will create more
				      * debugging output during the
				      * multigrid cycles.
				      */
    void set_debug(unsigned int);
    
  private:
    
				     /**
				      * The V-cycle multigrid method.
				      * <tt>level</tt> is the level the
				      * function starts on. It
				      * will usually be called for the
				      * highest level from outside,
				      * but will then call itself
				      * recursively for <tt>level-1</tt>,
				      * unless we are on #minlevel
				      * where the coarse grid solver
				      * solves the problem exactly.
				      */
    void level_v_step (const unsigned int level);

				     /**
				      * The actual W-cycle or F-cycle
				      * multigrid method.
				      * <tt>level</tt> is the level
				      * the function starts on. It
				      * will usually be called for the
				      * highest level from outside,
				      * but will then call itself
				      * recursively for
				      * <tt>level-1</tt>, unless we
				      * are on #minlevel where the
				      * coarse grid solver solves the
				      * problem exactly.
				      */
    void level_step (const unsigned int level, Cycle cycle);

				     /**
				      * Cycle type performed by the method cycle().
				      */
    Cycle cycle_type;
    
				     /**
				      * Level for coarse grid solution.
				      */
    unsigned int minlevel;

				     /**
				      * Highest level of cells.
				      */
    unsigned int maxlevel;

  public:
				     /**
				      * Input vector for the
				      * cycle. Contains the defect of
				      * the outer method projected to
				      * the multilevel vectors.
				      */
    MGLevelObject<VECTOR> defect;

				     /**
				      * The solution update after the
				      * multigrid step.
				      */
    MGLevelObject<VECTOR> solution;
    
  private:
				     /**
				      * Auxiliary vector.
				      */
    MGLevelObject<VECTOR> t;    

				     /**
				      * Auxiliary vector for W- and
				      * F-cycles. Left uninitialized
				      * in V-cycle.
				      */
    MGLevelObject<VECTOR> defect2;


				     /**
				      * The matrix for each level.
				      */
    SmartPointer<const MGMatrixBase<VECTOR> > matrix;
    
				     /**
				      * The matrix for each level.
				      */
    SmartPointer<const MGCoarseGridBase<VECTOR> > coarse;
    
				     /**
				      * Object for grid tranfer.
				      */
    SmartPointer<const MGTransferBase<VECTOR> > transfer;

				     /**
				      * The pre-smoothing object.
				      */
    SmartPointer<const MGSmootherBase<VECTOR> > pre_smooth;

				     /**
				      * The post-smoothing object.
				      */
    SmartPointer<const MGSmootherBase<VECTOR> > post_smooth;
    
				     /**
				      * Edge matrix from fine to coarse.
				      */
    SmartPointer<const MGMatrixBase<VECTOR> > edge_down;

				     /**
				      * Transpose edge matrix from coarse to fine.
				      */
    SmartPointer<const MGMatrixBase<VECTOR> > edge_up;

				     /**
				      * Level for debug
				      * output. Defaults to zero and
				      * can be set by set_debug().
				      */
    unsigned int debug;
    
    template<int dim, class VECTOR2, class TRANSFER> friend class PreconditionMG;
};


/**
 * Multi-level preconditioner.
 * Here, we collect all information needed for multi-level preconditioning
 * and provide the standard interface for LAC iterative methods.
 *
 * Furthermore, it needs functions <tt>void copy_to_mg(const VECTOR&)</tt>
 * to store @p src in the right hand side of the multi-level method and
 * <tt>void copy_from_mg(VECTOR&)</tt> to store the result of the v-cycle in @p dst.
 *
 * @author Guido Kanschat, 1999, 2000, 2001, 2002
 */
template<int dim, class VECTOR, class TRANSFER>
class PreconditionMG : public Subscriptor
{
  public:
				     /**
				      * Constructor.
				      * Arguments are the multigrid object,
				      * pre-smoother, post-smoother and
				      * coarse grid solver.
				      */
    PreconditionMG(const MGDoFHandler<dim>&     mg_dof,
		   Multigrid<VECTOR>&           mg,
		   const TRANSFER& transfer);

				     /**
				      * Dummy function needed by other classes.
				      */
    bool empty () const;
    
				     /**
				      * Preconditioning operator.
				      * Calls the @p vcycle function
				      * of the @p MG object passed to
				      * the constructor.
				      *
				      * This is the operator used by
				      * LAC iterative solvers.
				      */
    template<class VECTOR2>
    void vmult (VECTOR2       &dst,
		const VECTOR2 &src) const;
    
				     /**
				      * Preconditioning operator.
				      * Calls the @p vcycle function
				      * of the @p MG object passed to
				      * the constructor.
				      */
    template<class VECTOR2>
    void vmult_add (VECTOR2       &dst,
		    const VECTOR2 &src) const;
    
				     /**
				      * Tranposed preconditioning operator.
				      *
				      * Not implemented, but the
				      * definition may be needed.
				      */
    template<class VECTOR2>
    void Tvmult (VECTOR2       &dst,
		const VECTOR2 &src) const;
    
				     /**
				      * Tranposed preconditioning operator.
				      *
				      * Not implemented, but the
				      * definition may be needed.
				      */
    template<class VECTOR2>
    void Tvmult_add (VECTOR2       &dst,
		     const VECTOR2 &src) const;
    
  private:
				     /**
				      * Associated @p MGDoFHandler.
				      */
    SmartPointer<const MGDoFHandler<dim> > mg_dof_handler;

				     /**
				      * The multigrid object.
				      */
    SmartPointer<Multigrid<VECTOR> > multigrid;
    
				     /**
				      * Object for grid tranfer.
				      */
    SmartPointer<const TRANSFER> transfer;  
};

/*@}*/

/// @if NoDoc
/* --------------------------- inline functions --------------------- */


template <class VECTOR>
template <int dim>
Multigrid<VECTOR>::Multigrid (const MGDoFHandler<dim>& mg_dof_handler,
			      const MGMatrixBase<VECTOR>& matrix,
			      const MGCoarseGridBase<VECTOR>& coarse,
			      const MGTransferBase<VECTOR>& transfer,
			      const MGSmootherBase<VECTOR>& pre_smooth,
			      const MGSmootherBase<VECTOR>& post_smooth,
			      Multigrid<VECTOR>::Cycle cycle)
		:
		cycle_type(cycle),
		minlevel(0),
		maxlevel(mg_dof_handler.get_tria().n_levels()-1),
		defect(minlevel,maxlevel),
		solution(minlevel,maxlevel),
		t(minlevel,maxlevel),
		defect2(minlevel,maxlevel),
		matrix(&matrix, typeid(*this).name()),
		coarse(&coarse, typeid(*this).name()),
		transfer(&transfer, typeid(*this).name()),
		pre_smooth(&pre_smooth, typeid(*this).name()),
		post_smooth(&post_smooth, typeid(*this).name()),
		edge_down(0, typeid(*this).name()),
		edge_up(0, typeid(*this).name()),
		debug(0)
{}



template <class VECTOR>
inline
unsigned int
Multigrid<VECTOR>::get_maxlevel () const
{
  return maxlevel;
}



template <class VECTOR>
inline
unsigned int
Multigrid<VECTOR>::get_minlevel () const
{
  return minlevel;
}


/* --------------------------- inline functions --------------------- */


template<int dim, class VECTOR, class TRANSFER>
PreconditionMG<dim, VECTOR, TRANSFER>
::PreconditionMG(const MGDoFHandler<dim>&      mg_dof_handler,
		 Multigrid<VECTOR>& mg,
		 const TRANSFER&              transfer)
		:
  mg_dof_handler(&mg_dof_handler),
  multigrid(&mg),
  transfer(&transfer)
{}

template<int dim, class VECTOR, class TRANSFER>
inline bool
PreconditionMG<dim, VECTOR, TRANSFER>::empty () const
{
  return false;
}

template<int dim, class VECTOR, class TRANSFER>
template<class VECTOR2>
void
PreconditionMG<dim, VECTOR, TRANSFER>::vmult (
  VECTOR2& dst,
  const VECTOR2& src) const
{
  transfer->copy_to_mg(*mg_dof_handler,
		       multigrid->defect,
		       src);
  multigrid->cycle();
  transfer->copy_from_mg(*mg_dof_handler,
			 dst,
			 multigrid->solution);
}


template<int dim, class VECTOR, class TRANSFER>
template<class VECTOR2>
void
PreconditionMG<dim, VECTOR, TRANSFER>::vmult_add (
  VECTOR2& dst,
  const VECTOR2& src) const
{
  transfer->copy_to_mg(*mg_dof_handler,
		       multigrid->defect,
		       src);
  multigrid->cycle();
  transfer->copy_from_mg_add(*mg_dof_handler,
			     dst,
			     multigrid->solution);
}


template<int dim, class VECTOR, class TRANSFER>
template<class VECTOR2>
void
PreconditionMG<dim, VECTOR, TRANSFER>::Tvmult (
  VECTOR2&,
  const VECTOR2&) const
{
  Assert(false, ExcNotImplemented());
}


template<int dim, class VECTOR, class TRANSFER>
template<class VECTOR2>
void
PreconditionMG<dim, VECTOR, TRANSFER>::Tvmult_add (
  VECTOR2&,
  const VECTOR2&) const
{
  Assert(false, ExcNotImplemented());
}

/// @endif

#endif
