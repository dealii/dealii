//----------------------------  multigrid.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  multigrid.h  ---------------------------
#ifndef __deal2__multigrid_h
#define __deal2__multigrid_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <multigrid/mg_base.h>
#include <lac/forward_declarations.h>
#include <grid/forward_declarations.h>

#include <vector>


/**
 * Implementation of the multigrid method.
 *
 * The function actually performing a multi-level cycle,
 * @p{level_mgstep}, as well as the function @p{vcycle}, calling it,
 * require several helper classes handed over as template parameters.
 * These classes have to meet the following requirements:
 *
 * @p{MATRIX} is a matrix as specified for LAC-solvers, that is, it
 * has a function
 * \begin{verbatim}
 *   void MATRIX::vmult(VECTOR& x, const VECTOR& b)
 * \end{verbatim}
 * performing the matrix vector product @p{x=Ab}.
 *
 * @p{SMOOTH} is a class with a function
 * \begin{verbatim}
 *   void SMOOTH::smooth (unsigned int level, VECTOR& x, const VECTOR& b)
 * \end{verbatim}
 * modifying the vector @p{x} such that the residual @{b-Ax}
 * on level @p{l} is smoothened. Refer to @ref{MGSmootherRelaxation}
 * for an example.
 *
 * @p{COARSE} has an
 * \begin{verbatim}
 *   void COARSE::operator() (VECTOR& x, const VECTOR& b)
 * \end{verbatim}
 * returning the solution to the system @p{Ax=b} on the coarsest level
 * in @p{x}.
 * 
 * @author Guido Kanschat, 1999, 2001
 */
template <class VECTOR>
class Multigrid : public Subscriptor
{
  public:
  typedef VECTOR vector_type;
  typedef const VECTOR const_vector_type;
  
				     /**
				      * Constructor. The
				      * @p{MGDoFHandler} is used to
				      * determine the highest possible
				      * level. @p{transfer} is an
				      * object performing prolongation
				      * and restriction.
				      *
				      * The V-cycle will start on
				      * level @p{maxlevel} and goes
				      * down to level @p{minlevel},
				      * where the coarse grid solver
				      * will be used.
				      *
				      * This function already
				      * initializes the vectors which
				      * will be used later on in the
				      * course of the
				      * computations. You should
				      * therefore create objects of
				      * this type as late as possible.
				      */
  template <int dim>
    Multigrid(const MGDoFHandler<dim>& mg_dof_handler,
	      const unsigned int            minlevel = 0,
	      const unsigned int            maxlevel = 1000000);
    

				     /**
				      * Execute one step of the
				      * v-cycle algorithm.  This
				      * function assumes, that the
				      * vector @p{defect} is properly
				      * filled with the residual in
				      * the outer defect correction
				      * scheme (usually performed by
				      * @ref{PreconditionMG}). After
				      * execution of @p{vcycle()}, the
				      * result is in the vector
				      * @p{solution}. See
				      * @p{copy_*_mg} in class
				      * @p{MGTools} if you want to use
				      * these vectors yourself.
				      *
				      * The actual work for this
				      * function is done in
				      * @p{level_mgstep}.
				      */
    template <class MATRIX, class TRANSFER, class SMOOTHER, class COARSE>
    void vcycle(const MGLevelObject<MATRIX>& matrix,
		const TRANSFER&              transfer,
		const SMOOTHER&              pre_smooth,
		const SMOOTHER&              post_smooth,
		const COARSE&                cgs);

				     /**
				      * Print a level vector using
				      * @ref{DoFHandler}.
				      */
/*      virtual void print_vector (const unsigned int level, */
/*  			       const VECTOR& v, */
/*  			       const char* name) const; */
    

  private:
    
				     /**
				      * The actual v-cycle multigrid method.
				      * @p{level} is the level the
				      * function should work on. It
				      * will usually be called for the
				      * highest level from outside,
				      * but will then call itself
				      * recursively for @p{level-1},
				      * unless we are on @p{minlevel}
				      * where instead of recursing
				      * deeper, the coarse grid solver
				      * is used to solve the matrix of
				      * this level.
				      */
    template <class MATRIX, class TRANSFER, class SMOOTHER, class COARSE>
    void level_mgstep (const unsigned int           level,
		       const MGLevelObject<MATRIX>& matrix,
		       const TRANSFER&              transfer,
		       const SMOOTHER&              pre_smooth,
		       const SMOOTHER&              post_smooth,
		       const COARSE&                cgs);
				     /**
				      * Level for coarse grid solution.
				      */
    unsigned int minlevel;

				     /**
				      * Highest level of cells.
				      */
    unsigned int maxlevel;
    
				     /**
				      * Auxiliary vector. Contains the
				      * defect to be corrected before
				      * the multigrid step.
				      */
    MGLevelObject<VECTOR> defect;

				     /**
				      * Auxiliary vector. Contains the
				      * updated solution after the
				      * multigrid step.
				      */
    MGLevelObject<VECTOR> solution;
    
				     /**
				      * Auxiliary vector.
				      */
    MGLevelObject<VECTOR> t;    


  private:
    

				     /**
				      * Exception.
				      */
    DeclException2(ExcSwitchedLevels, int, int,
		   << "minlevel and maxlevel switched, should be: "
		   << arg1 << "<=" << arg2);
    template<int dim, class MATRIX, class VECTOR2, class TRANSFER, class SMOOTHER, class COARSE> friend class PreconditionMG;
};


/**
 * Multi-level preconditioner.
 * Here, we collect all information needed for multi-level preconditioning
 * and provide the standard interface for LAC iterative methods.
 *
 * The template parameter class @p{MG} is required to inherit @p{MGBase}.
 * Furthermore, it needs functions @p{void copy_to_mg(const VECTOR&)}
 * to store @p{src} in the right hand side of the multi-level method and
 * @p{void copy_from_mg(VECTOR&)} to store the result of the v-cycle in @p{dst}.
 *
 * @author Guido Kanschat, 1999
 */
template<int dim, class MATRIX, class VECTOR, class TRANSFER, class SMOOTHER, class COARSE>
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
		   const MGLevelObject<MATRIX>& matrix,
		   const TRANSFER&              transfer,
		   const SMOOTHER&              pre,
		   const SMOOTHER&              post,
		   const COARSE&                coarse);
    
				     /**
				      * Preconditioning operator.
				      * Calls the @p{vcycle} function
				      * of the @p{MG} object passed to
				      * the constructor.
				      *
				      * This is the operator used by
				      * LAC iterative solvers.
				      */
    void vmult (VECTOR       &dst,
		const VECTOR &src) const;
    
				     /**
				      * Preconditioning operator.
				      * Calls the @p{vcycle} function
				      * of the @p{MG} object passed to
				      * the constructor.
				      */
  void vmult_add (VECTOR       &dst,
		  const VECTOR &src) const;
    
				     /**
				      * Tranposed preconditioning operator.
				      *
				      * Not implemented, but the
				      * definition may be needed.
				      */
    void Tvmult (VECTOR       &dst,
		const VECTOR &src) const;
    
				     /**
				      * Tranposed preconditioning operator.
				      *
				      * Not implemented, but the
				      * definition may be needed.
				      */
    void Tvmult_add (VECTOR       &dst,
		     const VECTOR &src) const;
    
  private:
				     /**
				      * Associated @p{MGDoFHandler}.
				      */
    SmartPointer<const MGDoFHandler<dim> > mg_dof_handler;

				     /**
				      * The multigrid object.
				      */
    SmartPointer<Multigrid<VECTOR> > multigrid;
    
				     /**
				      * Matrices for each level. The
				      * matrices are prepared by the
				      * constructor of @p{Multigrid} and
				      * can be accessed for
				      * assembling.
				      */
    SmartPointer<const MGLevelObject<MATRIX> > level_matrices;

				   /**
				    * Object for grid tranfer.
				    */
    SmartPointer<const TRANSFER> transfer;
  
				     /**
				      * The pre-smoothing object.
				      */
    SmartPointer<const SMOOTHER> pre;

				     /**
				      * The post-smoothing object.
				      */
    SmartPointer<const SMOOTHER> post;
    
				     /**
				      * The coarse grid solver.
				      */
    SmartPointer<const COARSE> coarse;
};




/* --------------------------- inline functions --------------------- */


template <class VECTOR>
template <int dim>
Multigrid<VECTOR>::Multigrid (const MGDoFHandler<dim>& mg_dof_handler,
			      const unsigned int min_level,
			      const unsigned int max_level)
		:
		minlevel(min_level),
				     maxlevel(std::min(mg_dof_handler.get_tria().n_levels()-1,
						       max_level)),
				     defect(minlevel,maxlevel),
				     solution(minlevel,maxlevel),
				     t(minlevel,maxlevel)
{
};


template <class VECTOR>
template <class MATRIX, class TRANSFER, class SMOOTHER, class COARSE>
void
Multigrid<VECTOR>::level_mgstep(const unsigned int        level,
				const MGLevelObject<MATRIX>& matrix,
				const TRANSFER&              transfer,
				const SMOOTHER     &pre_smooth,
				const SMOOTHER     &post_smooth,
				const COARSE &coarse_grid_solver)
{
#ifdef MG_DEBUG
  char *name = new char[100];
  sprintf(name, "MG%d-defect",level);
  print_vector(level, defect[level], name);
#endif

  solution[level] = 0.;
  
  if (level == minlevel)
    {
      coarse_grid_solver(level, solution[level], defect[level]);
#ifdef MG_DEBUG
      sprintf(name, "MG%d-solution",level);
      print_vector(level, solution[level], name);
#endif
      return;
    }

			   // smoothing of the residual by modifying s
  pre_smooth.smooth(level, solution[level], defect[level]);

#ifdef MG_DEBUG
  sprintf(name, "MG%d-pre",level);
  print_vector(level, solution[level], name);
#endif
  
				   // t = -A*solution[level]
  matrix[level].vmult(t[level], solution[level]);
  t[level].scale(-1);
  
				   // make t rhs of lower level
				   // The non-refined parts of the
				   // coarse-level defect already contain
				   // the global defect, the refined parts
				   // its restriction.
  for (unsigned int l = level;l>minlevel;--l)
    {
      t[l-1] = 0.;
      transfer.restrict_and_add (l, t[l-1], t[l]);
      defect[l-1] += t[l-1];
    }

				   // add additional DG contribution
//  edge_vmult(level, defect[level-1], defect[level]);
  
				   // do recursion
  level_mgstep(level-1, matrix, transfer, pre_smooth, post_smooth, coarse_grid_solver);

				   // reset size of the auxiliary
				   // vector, since it has been
				   // resized in the recursive call to
				   // level_mgstep directly above
  t[level] = 0.;

				   // do coarse grid correction

  transfer.prolongate(level, t[level], solution[level-1]);

#ifdef MG_DEBUG
  sprintf(name, "MG%d-cgc",level);
  print_vector(level, t, name);
#endif

  solution[level] += t[level];
  
				   // post-smoothing
  post_smooth.smooth(level, solution[level], defect[level]);

#ifdef MG_DEBUG
  sprintf(name, "MG%d-post",level);
  print_vector(level, solution[level], name);

  delete[] name;
#endif
}


template <class VECTOR>
template <class MATRIX, class TRANSFER, class SMOOTHER, class COARSE>
void
Multigrid<VECTOR>::vcycle(const MGLevelObject<MATRIX>& matrix,
			  const TRANSFER&              transfer,
			  const SMOOTHER     &pre_smooth,
			  const SMOOTHER     &post_smooth,
			  const COARSE &coarse_grid_solver)
{
				   // The defect vector has been
				   // initialized by copy_to_mg. Now
				   // adjust the other vectors.
  solution.resize(minlevel, maxlevel);
  t.resize(minlevel, maxlevel);
  
  for (unsigned int level=minlevel; level<=maxlevel;++level)
    {
      solution[level].reinit(defect[level]);
      t[level].reinit(defect[level]);
    }

  level_mgstep (maxlevel, matrix, transfer, pre_smooth, post_smooth, coarse_grid_solver);
//  abort ();
}


/* ------------------------------------------------------------------- */


template<int dim, class MATRIX, class VECTOR, class TRANSFER, class SMOOTHER, class COARSE>
PreconditionMG<dim, MATRIX, VECTOR, TRANSFER, SMOOTHER, COARSE>
::PreconditionMG(const MGDoFHandler<dim>&      mg_dof_handler,
		 Multigrid<VECTOR>& mg,
		 const MGLevelObject<MATRIX>& matrix,
		 const TRANSFER&              transfer,
		 const SMOOTHER&              pre,
		 const SMOOTHER&              post,
		 const COARSE&                coarse)
		:
  mg_dof_handler(&mg_dof_handler),
  multigrid(&mg),
  level_matrices(&matrix),
  transfer(&transfer),
  pre(&pre),
  post(&post),
  coarse(&coarse)
{}


template<int dim, class MATRIX, class VECTOR, class TRANSFER, class SMOOTHER, class COARSE>
void
PreconditionMG<dim, MATRIX, VECTOR, TRANSFER, SMOOTHER, COARSE>::vmult (
  VECTOR& dst,
  const VECTOR& src) const
{
  transfer->copy_to_mg(*mg_dof_handler,
		       multigrid->defect,
		       src);
  
  multigrid->vcycle(*level_matrices,
		    *transfer,
		    *pre, *post, *coarse);

  transfer->copy_from_mg(*mg_dof_handler,
			 dst,
			 multigrid->solution);
}


template<int dim, class MATRIX, class VECTOR, class TRANSFER, class SMOOTHER, class COARSE>
void
PreconditionMG<dim, MATRIX, VECTOR, TRANSFER, SMOOTHER, COARSE>::vmult_add (
  VECTOR& dst,
  const VECTOR& src) const
{
  transfer->copy_to_mg(*mg_dof_handler,
		       multigrid->defect,
		       src);
  
  multigrid->vcycle(*level_matrices,
		    *transfer,
		    *pre, *post, *coarse);

  transfer->copy_from_mg_add(*mg_dof_handler,
			     dst,
			     multigrid->solution);
}


template<int dim, class MATRIX, class VECTOR, class TRANSFER, class SMOOTHER, class COARSE>
void
PreconditionMG<dim, MATRIX, VECTOR, TRANSFER, SMOOTHER, COARSE>::Tvmult (
  VECTOR&,
  const VECTOR&) const
{
  Assert(false, ExcNotImplemented());
}


template<int dim, class MATRIX, class VECTOR, class TRANSFER, class SMOOTHER, class COARSE>
void
PreconditionMG<dim, MATRIX, VECTOR, TRANSFER, SMOOTHER, COARSE>::Tvmult_add (
  VECTOR&,
  const VECTOR&) const
{
  Assert(false, ExcNotImplemented());
}


#endif
