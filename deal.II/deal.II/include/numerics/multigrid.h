/*----------------------------   multigrid.h     ---------------------------*/
/*      $Id$                 */
#ifndef __multigrid_H
#define __multigrid_H
/*----------------------------   multigrid.h     ---------------------------*/

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
class MGMatrix
{
  MGMatrix(const MGMatrix&);
  operator=(const MGMatrix&);
public:
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
