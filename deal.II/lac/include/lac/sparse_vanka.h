/*----------------------------   sparse_vanka.h     ---------------------------*/
/*      $Id$                 */
#ifndef __sparse_vanka_H
#define __sparse_vanka_H
/* Copyright Guido Kanschat, 1999 */
/*----------------------------   sparse_vanka.h     ---------------------------*/



#include <base/smartpointer.h>
#include <lac/forward_declarations.h>

#include <vector>
#include <map>


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
 * The right hand side is built up in the same way, i.e. by copying
 * all entries that coupled with the one under present consideration,
 * but it is augmented by all degrees of freedom coupling with the
 * degrees from the set described above (i.e. the DoFs coupling second
 * order to the present one). The reason for this is, that the local
 * problems to be solved shall have Dirichlet boundary conditions on
 * the second order coupling DoFs, so we have to take them into
 * account but eliminate them before actually solving; this
 * elimination is done by the modification of the right hand side, and
 * in the end these degrees of freedom do not occur in the matrix and
 * solution vector any more at all.
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
 *
 * \subsubsection{Implementor's remark}
 * At present, the local matrices are built up such that the degree of freedom
 * associated with the local Lagrange multiplier is the first one. Thus, usually
 * the upper left entry in the local matrix is zero. It is not clear to me (W.B.)
 * whether this might pose some problems in the inversion of the local matrices.
 * Maybe someone would like to check this.
 *
 * @author Guido Kanschat, documentation and extensions by Wolfgang Bangerth; 1999, 2000
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
				      *
				      * The matrix #M# which is passed
				      * here may or may not be the
				      * same matrix for which this
				      * object shall act as
				      * preconditioner. In particular,
				      * it is conceivable that the
				      * preconditioner is build up for
				      * one matrix once, but is used
				      * for subsequent steps in a
				      * nonlinear process as well,
				      * where the matrix changes in
				      * each step slightly.
				      *
				      * If #conserve_mem# is #false#,
				      * then the inverses of the local
				      * systems are computed now, if
				      * the flag is #true#, then they
				      * are computed every time the
				      * preconditioner is
				      * applied. This saves some
				      * memory, but makes
				      * preconditioning very
				      * slow. Note also, that if the
				      * flag is #false#, the the
				      * contents of the matrix #M# at
				      * the time of calling this
				      * constructor are used, while if
				      * the flag is #true#, then the
				      * values in #M# at the time of
				      * preconditioning are used. This
				      * may lead to different results,
				      * obviously, of #M# changes.
				      *
				      * The parameter #n_threads#
				      * determines how many threads
				      * shall be used in parallel when
				      * building the inverses of the
				      * diagonal blocks. This
				      * parameter is ignored if not in
				      * multithreaded mode.
				      */
    SparseVanka(const SparseMatrix<number> &M,
		const vector<bool>         &selected,
		const bool                  conserve_memory = false,
		const unsigned int          n_threads       = 1);
    
				     /**
				      * Destructor.
				      * Delete all allocated matrices.
				      */
    ~SparseVanka();

				     /**
				      * Do the preconditioning.
				      * This function takes the residual
				      * in #src# and returns the resulting
				      * update vector in #dst#.
				      */
    template<typename number2>
    void operator() (Vector<number2>       &dst,
		     const Vector<number2> &src) const;

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
    const bool conserve_mem;

				     /**
				      * Number of threads to be used
				      * when building the
				      * inverses. Only relevant in
				      * multithread mode.
				      */
    const unsigned int n_threads;
    
				     /**
				      * Array of inverse matrices,
				      * one for each degree of freedom.
				      * Only those elements will be used
				      * that are tagged in #selected#.
				      */
    mutable vector<SmartPointer<FullMatrix<float> > > inverses;

				     /**
				      * Compute the inverses of all
				      * selected diagonal elements.
				      */
    void compute_inverses ();

				     /**
				      * Compute the inverses at
				      * positions in the range
				      * #[begin,end)#. In
				      * non-multithreaded mode,
				      * #compute_inverses()# calls
				      * this function with the whole
				      * range, but in multithreaded
				      * mode, several copies of this
				      * function are spawned.
				      */
    void compute_inverses (const unsigned int begin,
			   const unsigned int end);
    
				     /**
				      * Compute the inverse of the
				      * block located at position
				      * #row#. Since the map is used
				      * quite often, it is generated
				      * only once in the caller of
				      * this function and passed to
				      * this function which first
				      * clears it. Reusing the map
				      * makes the process
				      * significantly faster than in
				      * the case where this function
				      * re-creates it each time.
				      */
    void compute_inverse (const unsigned int               row,
			  map<unsigned int, unsigned int> &local_index);
};






/*----------------------------   sparse_vanka.h     ---------------------------*/
/* end of #ifndef __sparse_vanka_H */
#endif
/*----------------------------   sparse_vanka.h     ---------------------------*/
