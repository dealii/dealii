/*-------------------- precondition.h --------------------*/
//$Id$
// Guido Kanschat, University of Heidelberg, 1999

#ifndef __lac_precondition_h
#define __lac_precondition_h

/**
 * No preconditioning.
 * This class halps you, if you want to use a linear solver without
 * preconditioning. Since this is a strange idea, the documentation
 * here stays quite short.
 *
 * @author Guido Kanschat, 1999
 */
template<class VECTOR>
class PreconditionIdentity
{
  public:
  				     /**
				      * Execute preconditioning.
				      */
    void operator() (VECTOR&, const VECTOR&) const;
};


/**
 * Preconditioner using a matrix-builtin function.
 * This class forms a preconditioner suitable for the LAC solver
 * classes. Since many preconditioning methods are based on matrix
 * entries, these have to be implemented as member functions of the
 * underlying matrix implementation. This class now is intended to
 * allow easy access to these member functions from LAC solver
 * classes.
 *
 * It seems that all builtin preconditioners have a relaxation
 * parameter, so please use #PreconditionRelaxation# for these.
 * @author Guido Kanschat, 1999
 */ 

template<class MATRIX, class VECTOR>
class PreconditionUseMatrix
{
  public:
				     /**
				      * Type of the preconditioning
				      * function of the matrix.
				      */
    typedef void ( MATRIX::* function_ptr)(VECTOR&, const VECTOR&) const;
    
				     /**
				      * Constructor.
				      * This constructor stores a
				      * reference to the matrix object
				      * for later use and selects a
				      * preconditioning method, which
				      * must be a member function of
				      * that matrix.
				      */
    PreconditionUseMatrix(const MATRIX & M, function_ptr method);
    
				     /**
				      * Execute preconditioning.
				      */
    void operator() (VECTOR&, const VECTOR&) const;

  private:
				     /**
				      * Pointer to the matrix in use.
				      */
    const MATRIX& matrix;
				     /**
				      * Pointer to the preconditioning
				      * function.
				      */
    const function_ptr precondition;
};

/**
 * Preconditioner for builtin relaxation methods.
 * Application of this preconditioner involves
 * use of the #precondition_...# methods of #SparseMatrix#.
 *
 * Construction of objects of this class is quite strange, so read the
 * documentation of the constructor for an example.
 * @author Guido Kanschat, 1999
 */
template<class MATRIX, class VECTOR>
class PreconditionRelaxation
{
  public:
				     /**
				      * Type of the preconditioning
				      * function of the matrix.
				      */
    typedef void ( MATRIX::* function_ptr)(VECTOR&, const VECTOR&,
					   typename MATRIX::value_type) const;
    
				     /**
				      * Constructor.
				      * This constructor stores a
				      * reference to the matrix object
				      * for later use and selects a
				      * preconditioning method, which
				      * must be a member function of
				      * that matrix. An example of a
				      * typical usage is this:
	* \begin{verbatim}
	* SparseMatrix<float> A;
	* Vector<double> u;
	* PreconditioningRelaxation<SparseMatrix<float>,Vector<double> >
	*    precondition_ssor(A, &SparseMatrix<float>::template
	*       precondition_SSOR<double>, 1.2);
	* \end{verbatim}
	*/
    PreconditionRelaxation(const MATRIX & M, function_ptr method,
			   double omega = 1.);
    
				     /**
				      * Execute preconditioning.
				      */
    void operator() (VECTOR&, const VECTOR&) const;

  private:
				     /**
				      * Pointer to the matrix in use.
				      */
    const MATRIX& matrix;
				     /**
				      * Pointer to the preconditioning
				      * function.
				      */
    const function_ptr precondition;
				     /**
				      * Relaxation parameter.
				      */
    double omega;
};


template<class VECTOR>
void
PreconditionIdentity<VECTOR>::operator() (VECTOR& dst, const VECTOR& src) const
{
  dst = src;
}

template<class MATRIX, class VECTOR>
PreconditionUseMatrix<MATRIX, VECTOR>::PreconditionUseMatrix(const MATRIX& M,
							     function_ptr method)
		:
		matrix(M), precondition(method)
{}


template<class MATRIX, class VECTOR>
void
PreconditionUseMatrix<MATRIX, VECTOR>::operator() (VECTOR& dst,
						   const VECTOR& src) const
{
  (matrix.*precondition)(dst, src);
}

template<class MATRIX, class VECTOR>
PreconditionRelaxation<MATRIX, VECTOR>::PreconditionRelaxation(const MATRIX& M,
							     function_ptr method,
							       double omega)
		:
		matrix(M), precondition(method), omega(omega)
{}


template<class MATRIX, class VECTOR>
void
PreconditionRelaxation<MATRIX, VECTOR>::operator() (VECTOR& dst,
						   const VECTOR& src) const
{
  (matrix.*precondition)(dst, src, omega);
}

#endif
