/*----------------------------   sparse_vanka.h     ---------------------------*/
/*      $Id$                 */
#ifndef __sparse_vanka_H
#define __sparse_vanka_H
/* Copyright Guido Kanschat, 1999 */
/*----------------------------   sparse_vanka.h     ---------------------------*/



#include <base/smartpointer.h>
#include <lac/forward-declarations.h>

#include <vector>


/**
 * Point-wise Vanka preconditioning.
 * This class does Vanka preconditioning  on a point-wise base.
 * Vanka preconditioners are used for saddle point problems like Stoke's
 * problem or problems arising in optimization where Lagrange multiplier
 * occur and let Netwon's matrix have a zero block. With these matrices the
 * application of Jacobi or Gauss-Seidel methods is impossible, because
 * some diagonal elements are zero in the rows of the Lagrange multiplier.
 * The approach of Vanka is to solve a small (usually indefinite) system
 * of equations for each Langrange multiplie variable (we will also call
 * the pressure in Stoke's equation a Langrange multiplier since it
 * can be interpreted as such).
 *
 * Objects of this class are constructed by passing a vector of indices
 * of the degrees of freedom of the Lagrange multiplier. In the actual
 * preconditioning method, these rows are traversed in the order in which
 * the appear in the matrix. Since this is a Gauﬂ-Seidel like procedure,
 * remember to have a good ordering in advance (for transport dominated
 * problems, Cuthill-McKee algorithms are a good means for this, if points
 * on the inflow boundary are chosen as starting points for the renumbering).
 *
 * For each selected degree of freedom, a local system of equations is built
 * by the degree of freedom itself and all other values coupling immediately,
 * i.e. the set of degrees of freedom considered for the local system of
 * equations for degree of freedom #i# is #i# itself and all #j# such that
 * the element #(i,j)# is a nonzero entry in the sparse matrix under 
 * consideration. The elements #(j,i)# are not considered. We now pick all
 * matrix entries from rows and columns out of the set of degrees of freedom
 * just described out of the global matrix and put it into a local matrix,
 * which is subsequently inverted. This system may be of different size for
 * each degree of freedom, depending for example on the local neighborhood of
 * the respective node on a computational grid.
 *
 * The right hand side is built up in the same way, i.e. by copying all entries
 * that coupled with the one under present consideration, but it is augmented
 * by all degrees of freedom coupling with the degrees from the set described
 * above (i.e. the DoFs coupling second order to the present one).
 *
 * This local system is solved and the values are updated into the
 * destination vector.
 *
 * Remark: the Vanka method is a non-symmetric preconditioning method.
 *
 *
 * \subsection{Example of Use}
 * This little example is taken from a program doing parameter optimization.
 * The Lagrange multiplier is the third component of the finite element
 * used. The system is solved by the GMRES method.
 * \begin{verbatim}
 *                        // tag the Lagrange multiplier variable
 *    vector<bool> signature(3);
 *    signature[0] = signature[1] = false;
 *    signature[2] = true;
 *
 *                        // tag all dofs belonging to the
 *                        // Lagrange multiplier
 *    vector<bool> selected_dofs (dof.n_dofs(), false);
 *    DoFTools::extract_dofs(dof, signature, p_select);
 *                        // create the Vanka object
 *    SparseVanka<double> vanka (global_matrix, selected_dofs);
 *
 *                        // create the solver
 *    SolverGMRES<PreconditionedSparseMatrix<double>,
 *                Vector<double> >    gmres(control,memory,504);
 *    
 *			  // solve
 *    gmres.solve (global_matrix, solution, right_hand_side,
 *	           vanka);
 * \end{verbatim}
 *
 * @author Guido Kanschat, 1999
 */
template<typename number>
class SparseVanka
{
  public:
				     /**
				      * Constructor. Gets the matrix
				      * for preconditioning and a bit
				      * vector with entries #true# for
				      * all rows to be updated. A
				      * reference to this vector will
				      * be stored, so it must persist
				      * longer than the Vanka
				      * object. The same is true for
				      * the matrix.
				      */
    SparseVanka(const SparseMatrix<number> &M,
		const vector<bool>         &selected);
    
				     /**
				      * Destructor.
				      * Delete all allocated matrices.
				      */
    ~SparseVanka();
    
				     /**
				      * Do the preconditioning. This
				      * function contains a dispatch
				      * mechanism to use the
				      * multiplicative version by
				      * default and the additive version
				      * if requested by #set_additive#.
				      */
    template<typename number2>
    void operator() (Vector<number2>       &dst,
		     const Vector<number2> &src) const;

				     /**
				      * Application of the Vanka operator.
				      * This function takes the residual
				      * in #src# and returns the resulting
				      * update vector in #dst#.
				      */
    template<typename number2>
    void forward (Vector<number2>       &dst,
		  const Vector<number2> &src) const;
    
				     /**
				      * Application of the transpose
				      * Vanka operator.
				      * This function takes the residual
				      * in #src# and returns the resulting
				      * update vector in #dst#.
				      */
    template<typename number2>
    void backward (Vector<number2>       &dst,
		   const Vector<number2> &src) const;
    
				     /**
				      * Minimize memory consumption.
				      * Activating this option reduces
				      * memory needs of the Vanka object
				      * to nearly zero. You pay for this
				      * by a high increase of computing
				      * time, since all local matrices
				      * are built up and inverted every
				      * time the Vanka operator is applied.
				      */
    void conserve_memory();

				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotSquare);
    
  private:
				     /**
				      * Pointer to the matrix.
				      */
    SmartPointer<const SparseMatrix<number> > matrix;
    
				     /**
				      * Indices of Lagrange
				      * multipliers.
				      */
    const vector<bool> &selected;
    
				     /**
				      * Conserve memory flag.
				      */
    bool conserve_mem;
    
				     /**
				      * Array of inverse matrices,
				      * one for each degree of freedom.
				      * Only those elements will be used
				      * that are tagged in #selected#.
				      */
    mutable vector<SmartPointer<FullMatrix<float> > > inverses;
};



/* ---------------------------- Inline functions -----------------------*/

template<typename number>
template<typename number2>
inline
void
SparseVanka<number>::operator() (Vector<number2>& dst,
				 const Vector<number2>& src) const
{
  dst = 0.;
  forward(dst, src);
}




/*----------------------------   sparse_vanka.h     ---------------------------*/
/* end of #ifndef __sparse_vanka_H */
#endif
/*----------------------------   sparse_vanka.h     ---------------------------*/
