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
 * @author Guido Kanschat, Wolfgang Bangerth; 1999, 2000
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
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidRange,
		    unsigned int, unsigned int,
		    << "The bounds [" << arg1 << ',' << arg2
		    << ") do not form a valid range.");
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidVectorSize,
		    unsigned int, unsigned int,
		    << "The dimensions of vectors and matrices, "
		    << arg1 << " and " << arg2 << " do not match.");

  protected:
				     /**
				      * Apply the inverses in the
				      * range #[begin,end)# to the
				      * #src# vector and move the
				      * result into #dst#. Actually,
				      * only values of #src# from
				      * within the range are taken
				      * (all others are set to zero),
				      * and only values inside the
				      * range are written to #dst#, so
				      * the application of this
				      * function only does what is
				      * announced in the general
				      * documentation if the given
				      * range is the whole interval.
				      *
				      * The reason for providing the
				      * interval anyway is that in
				      * derived classes we may want to
				      * apply the preconditioner to
				      * blocks of the matrix only, in
				      * order to parallelize the
				      * application. Then, it is
				      * important to only write to
				      * some slices of #dst# and only
				      * takes values from similar
				      * slices of #src#, in order to
				      * eliminate the dependencies of
				      * threads of each other.
				      *
				      * The #operator()# of this class
				      * of course calls this function
				      * with the whole interval
				      * #[begin,end)=[0,matrix.m())#.
				      */
    template<typename number2>
    void apply_preconditioner (Vector<number2>       &dst,
			       const Vector<number2> &src,
			       const unsigned int     begin,
			       const unsigned int     end) const;    
    
  private:
				     /**
				      * Pointer to the matrix.
				      */
    SmartPointer<const SparseMatrix<number> > matrix;
    
				     /**
				      * Conserve memory flag.
				      */
    const bool conserve_mem;

				     /**
				      * Indices of those degrees of
				      * freedom that we shall work on.
				      */
    const vector<bool> &selected;

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



/**
 * Block version of the sparse Vanka preconditioner. This class
 * divides the matrix into blocks and works on the diagonal blocks
 * only, which of course reduces the efficiency as preconditioner, but
 * is perfectly parallelizable. The constructor takes a parameter into
 * how many diagonal blocks the matrix shall be subdivided and then
 * lets the underlying class do the work.
 *
 * Division of the matrix is done in a way such that the blocks are
 * not necessarily of equal size, but such that the number of selected
 * degrees of freedom for which a local system is to be solved is
 * equal between blocks. The reason for this strategy to subdivision
 * is load-balancing for multithreading, but it is necessary to note
 * that this almost renders the capability as precondition useless if
 * the degrees of freedom are numbered by component, i.e. all Lagrange
 * multipliers en bloc.
 *
 * This class is probably useless if you don't have a multiprocessor
 * system, since then the amount of work per preconditioning step is
 * the same as for the #SparseVanka# class, but preconditioning
 * properties are worse. On the other hand, if you have a
 * multiprocessor system, the worse preconditioning quality (leading
 * to more iterations of the linear solver) usually is well balanced
 * by the increased speed of application due to the parallelization,
 * leading to an overall decrease in elapsed wall-time for solving
 * your linear system. It should be noted that the quality as
 * preconditioner reduces with growing number of blocks, so there may
 * be an optimal value (in terms of wall-time per linear solve) for
 * the number of blocks.
 *
 * To facilitate writing portable code, if the number of blocks into
 * which the matrix is to be subdivided, is set to one, then this
 * class acts just like the #SparseVanka# class. You may therefore
 * want to set the number of blocks equal to the number of processors
 * you have.
 *
 * Note that the parallelization is done if #deal.II# was configured
 * for multithread use and that the number of threads which is spawned
 * equals the number of blocks. This is reasonable since you will not
 * want to set the number of blocks unnecessarily large, since, as
 * mentioned, this reduces the preconditioning properties.
 *
 *
 * \subsection{Typical results}
 *
 * As a prototypical test case, we use a nonlinear problem from
 * optimization, which leads to a series of saddle point problems,
 * each of which is solved using GMRES with Vanka as
 * preconditioner. The equation had approx. 850 degrees of
 * freedom. With the non-blocked version #SparseVanka# (or
 * #SparseBlockVanka# with #n_blocks==1#), the following numbers of
 * iterations is needed to solver the linear system in each nonlinear
 * step: \begin{verbatim} 101 68 64 53 35 21 \end{verbatim} With four
 * blocks, we need the following numbers of iterations
 * \begin{verbatim} 124 88 83 66 44 28 \end{verbatim} As can be seen,
 * more iterations are needed. However, in terms of computing time,
 * the first version needs 72 seconds wall time (and 79 seconds CPU
 * time, which is more than wall time since some other parts of the
 * program were parallelized as well), while the second version needed
 * 53 second wall time (and 110 seconds CPU time) on a four processor
 * machine. The total time is in both cases dominated by the linear
 * solvers. In this case, it is therefore worth while using the
 * blocked version of the preconditioner if wall time is more
 * important than CPU time.
 *
 * @author Wolfgang Bangerth, 2000
 */
template<typename number>
class SparseBlockVanka : public SparseVanka<number>
{
  public:
				     /**
				      * Constructor. Pass all
				      * arguments except for
				      * #n_blocks# to the base class.
				      */
    SparseBlockVanka (const SparseMatrix<number> &M,
		      const vector<bool>         &selected,
		      const bool                  conserve_memory = false,
		      const unsigned int          n_threads       = 1,
		      const unsigned int          n_blocks        = 1);

				     /**
				      * Apply the preconditioner.
				      */
    template<typename number2>
    void operator() (Vector<number2>       &dst,
		     const Vector<number2> &src) const;
    
  private:
				     /**
				      * Store the number of blocks.
				      */
    const unsigned int n_blocks;

				     /**
				      * In this field, we precompute
				      * the first and the one after
				      * the last index of each
				      * block. This computation is
				      * done in the constructor, to
				      * avoid recomputing each time
				      * the preconditioner is called.
				      */
    vector<pair<unsigned int, unsigned int> > intervals;
};




/*----------------------------   sparse_vanka.h     ---------------------------*/
/* end of #ifndef __sparse_vanka_H */
#endif
/*----------------------------   sparse_vanka.h     ---------------------------*/
