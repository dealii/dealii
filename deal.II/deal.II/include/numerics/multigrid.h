/*----------------------------   multigrid.h     ---------------------------*/
/*      $Id$                 */
#ifndef __multigrid_H
#define __multigrid_H
/*----------------------------   multigrid.h     ---------------------------*/

class dVector;

/**
 * Vector with data for each level.
 *
 * This class provides auxiliary vectors for the multigrid method. It
 * seems, that it is used only locally in class #MultiGrid#, so do not
 * bother about it.
 *
 * In case your smoothing method needs auxiliary vectors, you should
 * use a memory class for each level separately. If memory is short,
 * one object for all levels could be a better choice.
 *
 * @author Guido Kanschat, 1999
 */
class MGVector
{
  MGVector(const MGVector&);
public:
				   /**
				    * Constructor using the number of
				    * levels.
				    */
  MGVector(unsigned l);
				   /**
				    * Access to data on a specific
				    level.
				   */
  dVector& operator[](unsigned l);

				   /**
				    * Access to data on a specific
				    level.
				   */
  const dVector& operator[](unsigned l) const;
  
};



/**
 * Basic matrix class for multigrid preconditioning.
 *
 * This matrix may be used in the iterative methods of LAC, where the
 * functions #vmult# and #precondition# and possibly their transposed
 * versions are needed.
 *
 * The functionality of the multigrid method is restricted to defect
 * correction. It is <B>not</B> iterative and the start solution is
 * always zero. Since by this <I>u<SUP>E</SUP><SUB>l</SUB></I> and
 * <I>u<SUP>A</SUP><SUB>l</SUB></I> (see report on multigrid) are
 * always zero, restriction is simplified a lot and maybe even the
 * seam condition on grids is oblivious. Still, I am not sure that
 * these restrictions on the class might cause numerical
 * inefficiencies.
 *
 * The function #precondition# is the actual multigrid method and
 * makes use of several operations to be implemented in derived
 * classes. It takes a defect in #src# and the result is the multigrid
 * preconditioned defect in #dst#.
 *
 * @author Guido Kanschat, 1999
 * */
class MultiGrid
{
  MultiGrid(const MultiGrid&);
  const MultiGrid& operator=(const MultiGrid&);

				   /** Auxiliary vectors.
				    */
  mutable MGVector d,s;

				   /**
				    * Highest level of cells.
				    */
  unsigned maxlevel;

				   /**
				    * Level for coarse grid solution.
				    */
  unsigned minlevel;
  
				   /** Tranfer from dVector to
				    * MGVector.
				    *
				    * This function copies data from a
				    * dVector, that is, data on the
				    * locally finest level, into the
				    * corresponding levels of an
				    * MGVector.
				    */
  void copy_to_mg(MGVector& dst, const dVector& src) const;

				   /**
				    * Transfer from MGVector to
				    * dVector.
				    *
				    * Copies data from active portions
				    * of an MGVector into the
				    * respective positions of a
				    * dVector. All other entries of
				    * #src# are zero.
				    */
  void copy_from_mg(dVector& dst, const MGVector& src) const;

				   /**
				    * The actual v-cycle multigrid method.
				    * This function is called on the
				    * highest level and recursively
				    * invokes itself down to the
				    * coarsest. There, it calls
				    * #coarse_grid_solution# and
				    * proceeds back up.
				    */
  void level_mgstep(unsigned level) const;
  
protected:
				   /**
				    * Number of pre-smoothing steps.
				    */
  unsigned n_pre_smooth;
  
				   /**
				    * Number of post-smoothing steps
				    */
  unsigned n_post_smooth;
				   /**
				    * The (pre-)smoothing algorithm.
				    * This function is required to
				    * perform #steps# iterations to
				    * smoothen the residual $Ax-b$.
				    *
				    * When overloading this function
				    * in derived classes, make sure
				    * that smoothing only applies to
				    * interior degrees of freedom of
				    * the actual level.
				    *
				    */
  virtual void smooth(unsigned level,
		      dVector& x,
		      const dVector& b,
		      unsigned steps) const;

				   /**
				    * The post-smoothing algorithm.
				    * Defaults to #smooth#.
				    */
  virtual void post_smooth(unsigned level,
			   dVector& dst,
			   const dVector& src,
			   unsigned steps) const;

				   /**
				    * Apply operator on all
				    * cells of a level.
				    *
				    */
  virtual void level_vmult(unsigned level,
			   dVector& dst,
			   const dVector& src) const;
				   /**
				    * Apply operator on non-refined
				    * cells.
				    *
				    * The sum over all levels
				    * of the results of this function
				    * is the multiplication with the
				    * normal fine-grid matrix.
				    */
  virtual void level_active_vmult(unsigned level,
				  dVector& dst,
				  const dVector& src) const;

				   /**
				    * Restriction from #level#.
				    *
				    * This function <B>adds</B> the
				    * restriction of #src# to #dst#,
				    * where #src# is a vector on #level#
				    * and #dst# is on #level#-1.
				    */
  virtual void restriction(unsigned level,
			   dVector& dst,
			   const dVector& src) const;
			   

				   /**
				    * Prolongation to #level#.
				    *
				    * <B>Adds</B> the prolongation of
				    * #src# to #dst#. Here, #dst# is on
				    * #level# and #src# on #level#-1.
				    */
  virtual void prolongation(unsigned level,
			    dVector& dst,
			    const dVector& src) const;

				   /**
				    * Solve exactly on coarsest grid.
				    */
  virtual void coarse_grid_solution(unsigned l,
				    dVector& dst,
				    const dVector& src) const;
  

  

  
public:
				   /**
				    * Constructor, subject to change.
				    */
  MultiGrid();

				   /** Matrix-vector multiplication.
				    *
				    * The global, non-multigrid
				    * matrix-vector multiplication
				    * used to compute the residual in
				    * the outer iteration.
				    *
				    */
  void vmult(dVector& dst, const dVector& src) const;
				   /** Multigrid preconditioning.
				    *
				    * This is the function, where the
				    * multigrid method is implemented.
				    *
				    */
  void precondition(dVector& dst, const dVector& src) const;
};




/*----------------------------   multigrid.h     ---------------------------*/
/* end of #ifndef __multigrid_H */
#endif
/*----------------------------   multigrid.h     ---------------------------*/
