/*----------------------------   multigrid.h     ---------------------------*/
/*      $Id$                 */
#ifndef __multigrid_H
#define __multigrid_H
/*----------------------------   multigrid.h     ---------------------------*/

class dVector;

/**
 * Vector with data for each level.
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
 * The function #precondition# is the actual multigrid method and
 * makes use of several operations to be implemented in derived
 * classes...
 *
 */
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
				    * dVector.
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
  unsigned n_presmooth;
  
				   /**
				    * Number of post-smoothing steps
				    */
  unsigned n_postsmooth;
				   /**
				    * The (pre-)smoothing algorithm.
				    */
  virtual void smooth(unsigned level,
		      dVector& dst,
		      const dVector& src) const;

				   /**
				    * The post-smoothing algorithm.
				    * Defaults to #smooth#.
				    */
  virtual void post_smooth(unsigned level,
		      dVector& dst,
		      const dVector& src) const;

				   /**
				    * Apply operator on non-refined
				    * cells.
				    */
  virtual void level_active_vmult(unsigned level,
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
				    */
  void precondition(dVector& dst, const dVector& src) const;
};




/*----------------------------   multigrid.h     ---------------------------*/
/* end of #ifndef __multigrid_H */
#endif
/*----------------------------   multigrid.h     ---------------------------*/
