/*----------------------------   matrices.h     ---------------------------*/
/*      $Id$                 */
#ifndef __matrices_H
#define __matrices_H
/*----------------------------   matrices.h     ---------------------------*/



#include <base/exceptions.h>

template <int dim> class Triangulation;
template <int dim> class DoFHandler;
template <int dim> class FiniteElement;
template <int dim> class FEValues;
template <int dim> class Quadrature;
template <int dim> class Function;
template <int dim> class Boundary;
template <int dim> class Equation;

class dVector;
class dFMatrix;
class dSMatrix;



/**
 * Provide a class which assembles certain standard matrices for a given
 * triangulation, using a given finite element and a quadrature formula.
 * All functions are static, so it is not necessary to create an object
 * of this type, though you may do so.
 *
 *
 * \subsection{Conventions for all functions}
 *
 * All functions take a sparse matrix object to hold the matrix to be
 * created. The functions assume that the matrix is initialized with a
 * sparsity pattern (#dSMatrixStruct#) corresponding to the given degree
 * of freedom handler, i.e. the sparsity structure is already as needed.
 * You can do this by calling the #DoFHandler<dim>::make_sparsity_pattern#
 * function.
 *
 * Furthermore it is assumed that no relevant data is in the matrix. All
 * entries will be overwritten. Entries which are not needed by the matrix
 * (and were thus added 'by hand' after #make_sparsity_pattern# was called)
 * are not touched and in special are not set to zero, so you have to care
 * yourself about that if you really need these entries.
 *
 *
 * \subsection{Supported matrices}
 *
 * At present there are functions to create the following matrices:
 * \begin{itemize}
 * \item #create_mass_matrix#: create the matrix with entries
 *   $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$. Uses the
 *   #MassMatrix# class. 
 *
 * \item #create_laplace_matrix#: there are two versions of this; the
 *   one which takes the #Function<dim># object creates
 *   $a_{ij} = \int_\Omega a(x) \nabla\phi_i(x) \nabla\phi_j(x) dx$,
 *   $a$ being the given function, while the other one assumes that
 *   $a=1$ which enables some optimzations. In fact the two versions
 *   are in one function, the coefficient being given as a defaulted
 *   argument, which is a pointer to a function and defaults to zero.
 *   This function uses the #LaplaceMatrix# class.
 * \end{itemize}
 */
template <int dim>
class MatrixCreator {
  public:
				     /**
				      * Assemble the mass matrix. See
				      * the general doc of this class
				      * for more information.
				      */
    static void create_mass_matrix (const DoFHandler<dim>    &dof,
				    const FiniteElement<dim> &fe,
				    const Quadrature<dim>    &q,
				    const Boundary<dim>      &boundary,
				    dSMatrix &matrix);

				     /**
				      * Assemble the laplace matrix with a
				      * variable weight function. If no 
				      * coefficient is given, it is assumed
				      * to be zero.
				      * 
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_laplace_matrix (const DoFHandler<dim>    &dof,
				       const FiniteElement<dim> &fe,
				       const Quadrature<dim>    &q,
				       const Boundary<dim>      &boundary,
				       dSMatrix &matrix,
				       const Function<dim>      *a = 0);
};





/**
 * Equation class to be passed to the #Assembler# if you want to make up the
 * mass matrix for your problem. The mass matrix is the matrix with
 * $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$.
 *
 * You may pass a coefficient function to the constructor. If you do so, the
 * assemble routines compute the matrix
 * $m_{ij} = \int_\Omega a(x) \phi_i(x) \phi_j(x) dx$
 * instead. The coefficient will in many cases be a strictly positive function.
 *
 * The class also has functions to create a right hand side
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$. The function $f(x)$ has to be
 * given to the constructor; if none is given, an error is issued if you
 * try to create a right hand side vector. The function to create right
 * hand side vectors is the same for all the matrix class in this file,
 * since it does not depend on the operator.
 *
 * The defaults for both right hand side and coefficient function is a
 * #NULL# pointer. If you need a coefficient but no right hand side object,
 * simply pass a #NULL# pointer to the constructor for its first argument.
 */
template <int dim>
class MassMatrix :  public Equation<dim> {
  public:
				     /**
				      * Constructor. Pass a function object if
				      * you want to create a right hand side
				      * vector, pass a function pointer (default
				      * is a NULL pointer). It is your duty to
				      * guarantee that the function object for
				      * the right hand side lives at least as
				      * long as this object does.
				      *
				      * You may also pass a function describing
				      * the weight to the integral (see the
				      * general docs for more information). The
				      * same applies for this object as said
				      * above.
				      */
    MassMatrix (const Function<dim> * const rhs = 0,
		const Function<dim> * const a = 0);

				     /**
				      * Assemble the cell matrix and right hand
				      * side vector for this cell. You need to
				      * give a right hand side object to the
				      * constructor to use this function. If
				      * a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &cell) const;

				     /**
				      * Construct the cell matrix for this cell.
				      * If a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &cell) const;

				     /**
				      * Only construct the right hand side
				      * vector for this cell. You need to give
				      * a right hand side function to the
				      * constructor in order to call this
				      * function.
				      */
    virtual void assemble (dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &cell) const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoRHSSelected);
    
  protected:
				     /**
				      * Pointer to a function describing the
				      * right hand side of the problem. Should
				      * be zero if not given to the constructor
				      * and should then not be used.
				      */
    const Function<dim> * const right_hand_side;

				     /**
				      * Pointer to a function describing the
				      * coefficient to the integral for the
				      * matrix entries. Should
				      * be zero if not given to the constructor
				      * and should then not be used.
				      */
    const Function<dim> * const coefficient;
};





/**
 * Equation class to be passed to the #Assembler# if you want to make up the
 * laplace matrix for your problem. The laplace matrix is the matrix with
 * $a_{ij} = \int_\Omega \nabla\phi_i(x) \cdot \nabla\phi_j(x) dx$.
 *
 * You may pass a coefficient function to the constructor. If you do so, the
 * assemble routines compute the matrix
 * $m_{ij} = \int_\Omega a(x) \nabla\phi_i(x) \cdot \nabla\phi_j(x) dx$
 * instead. The coefficient will in many cases be a strictly positive function.
 *
 * The class also has functions to create a right hand side
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$. The function $f(x)$ has to be
 * given to the constructor; if none is given, an error is issued if you
 * try to create a right hand side vector. The function to create right
 * hand side vectors is the same for all the matrix class in this file,
 * since it does not depend on the operator.
 *
 * The defaults for both right hand side and coefficient function is a
 * #NULL# pointer. If you need a coefficient but no right hand side object,
 * simply pass a #NULL# pointer to the constructor for its first argument.
 */
template <int dim>
class LaplaceMatrix :  public Equation<dim> {
  public:
				     /**
				      * Constructor. Pass a function object if
				      * you want to create a right hand side
				      * vector, pass a function pointer (default
				      * is a NULL pointer). It is your duty to
				      * guarantee that the function object for
				      * the right hand side lives at least as
				      * long as this object does.
				      *
				      * You may also pass a function describing
				      * the weight to the integral (see the
				      * general docs for more information). The
				      * same applies for this object as said
				      * above.
				      */
    LaplaceMatrix (const Function<dim> * const rhs = 0,
		   const Function<dim> * const a = 0);

				     /**
				      * Assemble the cell matrix and right hand
				      * side vector for this cell. You need to
				      * give a right hand side object to the
				      * constructor to use this function. If
				      * a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &cell) const;

				     /**
				      * Construct the cell matrix for this cell.
				      * If a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &cell) const;

				     /**
				      * Only construct the right hand side
				      * vector for this cell. You need to give
				      * a right hand side function to the
				      * constructor in order to call this
				      * function.
				      */
    virtual void assemble (dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &cell) const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcNoRHSSelected);
    
  protected:
				     /**
				      * Pointer to a function describing the
				      * right hand side of the problem. Should
				      * be zero if not given to the constructor
				      * and should then not be used.
				      */
    const Function<dim> * const right_hand_side;

    				     /**
				      * Pointer to a function describing the
				      * coefficient to the integral for the
				      * matrix entries. Should
				      * be zero if not given to the constructor
				      * and should then not be used.
				      */
    const Function<dim> * const coefficient;
};





/*----------------------------   matrices.h     ---------------------------*/
/* end of #ifndef __matrices_H */
#endif
/*----------------------------   matrices.h     ---------------------------*/
