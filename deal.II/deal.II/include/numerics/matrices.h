/*----------------------------   matrices.h     ---------------------------*/
/*      $Id$                 */
#ifndef __matrices_H
#define __matrices_H
/*----------------------------   matrices.h     ---------------------------*/



#include <base/exceptions.h>
#include <map>
#include <set>

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
 *   $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$. Here, the $\phi_i$
 *   are the basis functions of the finite element space given.
 *   This function uses the #MassMatrix# class. 
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
 *
 * All created matrices are `raw': they are not condensed, i.e. hanging
 * nodes are not eliminated. The reason is that you may want to add
 * several matrices and could then condense afterwards only once,
 * instead of for every matrix. To actually do computations with these
 * matrices, you have to condense the matrix using the
 * #ConstraintMatrix::condense# function; you also have to condense the
 * right hand side accordingly and distribute the solution afterwards.
 *
 * If you want to use boundary conditions, you have to use a function
 * like #ProblemBase<>::apply_dirichlet_bc# to matrix and right hand
 * side.
 *
 *
 * \subsection{Matrices on the boundary}
 *
 * The #create_boundary_mass_matrix# creates the matrix with entries
 * $m_{ij} = \int_{\Gamma} \phi_i \phi_j dx$, where $\Gamma$ is the union
 * of boundary parts with indicators contained in a set passed to the
 * function (i.e. if you want to set up the mass matrix for the parts of
 * the boundary with indicators zero and 2, you pass the function a set
 * of #unsigned char#s as parameter #boundary_parts# containing the elements
 * zero and 2). The $\phi_i$ are the basis functions which have at least
 * part of their support om $\Gamma$. The mapping between row and column
 * indices in the mass matrix and the right hand side and the global degree
 * of freedom numbers of the respective basis functions on the whole domain
 * is returned as a vector of numbers which has the same size as the dimension
 * of matrix and right hand side.
 *
 * Since in most cases we are not interested in the pure mass matrix on the
 * boundary, but rather need it to compute the projection of a function to
 * the boundary, no function is provided to only create the matrix.
 *
 * This function needs to get passed a matrix object to hold the resulting sparse
 * matrix. This object is supposed to be initialized with a suitable sparsity
 * pattern, which can be created using the
 * #DoFHandler<>::make_boundary_sparsity_pattern# function.
 *
 *
 * \subsection{Right hand sides}
 *
 * In many cases, you will not only want to build the matrix, but also
 * a right hand side, which will give a vector with
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$. For this purpose, each function
 * exists in two versions, one only building the matrix and one also
 * building the right hand side vector.
 *
 * Creation of the right hand side
 * is the same for all operators and therefore for all of the functions
 * below. It would be most orthogonal to write one single function which
 * builds up the right hand side and not provide many functions doing
 * the same thing. However, this may result in a heavy performance
 * penalty, since then many values of a certain finite element have to
 * be computed twice, so it is more economical to  implement it more than
 * once. If you only want to create a right hand side as above, there is
 * a function in the #VectorCreator# class. The use of the latter may be
 * useful if you want to create many right hand side vectors.
 */
template <int dim>
class MatrixCreator {
  public:
				     /**
				      *	Declare a data type which denotes a
				      *	mapping between a boundary indicator
				      *	and the function denoting the boundary
				      *	values on this part of the boundary.
				      *	Only one boundary function may be given
				      *	for each boundary indicator, which is
				      *	guaranteed by the #map# data type.
				      *	
				      *	See the general documentation of this
				      *	class for more detail.
				      */
    typedef map<unsigned char,const Function<dim>*> FunctionMap;

				     /**
				      * Assemble the mass matrix. If no 
				      * coefficient is given, it is assumed
				      * to be constant one.
				      * 
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_mass_matrix (const DoFHandler<dim>    &dof,
				    const FiniteElement<dim> &fe,
				    const Quadrature<dim>    &q,
				    const Boundary<dim>      &boundary,
				    dSMatrix                 &matrix,
				    const Function<dim>      *a = 0);

    				     /**
				      * Assemble the mass matrix and a right
				      * hand side vector. If no 
				      * coefficient is given, it is assumed
				      * to be constant one.
				      * 
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_mass_matrix (const DoFHandler<dim>    &dof,
				    const FiniteElement<dim> &fe,
				    const Quadrature<dim>    &q,
				    const Boundary<dim>      &boundary,
				    dSMatrix                 &matrix,
				    const Function<dim>      &rhs,
				    dVector                  &rhs_vector,
				    const Function<dim>      *a = 0);

				     /**
				      * Assemble the mass matrix and a right
				      * hand side vector along the boundary.
				      * If no 
				      * coefficient is given, it is assumed
				      * to be constant one.
				      *
				      * The matrix is assumed to already be
				      * initialized with a suiting sparsity
				      * pattern (the #DoFHandler# provides an
				      * appropriate function).
				      *
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_boundary_mass_matrix (const DoFHandler<dim>    &dof,
					     const FiniteElement<dim> &fe,
					     const Quadrature<dim-1>  &q,
					     const Boundary<dim>      &boundary,
					     dSMatrix                 &matrix,
					     const FunctionMap        &rhs,
					     dVector                  &rhs_vector,
					     vector<int>              &vec_to_dof_mapping,
					     const Function<dim>      *a = 0);

				     /**
				      * Assemble the Laplace matrix. If no 
				      * coefficient is given, it is assumed
				      * to be constant one.
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

				     /**
				      * Assemble the Laplace matrix and a right
				      * hand side vector. If no 
				      * coefficient is given, it is assumed
				      * to be constant one.
				      * 
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_laplace_matrix (const DoFHandler<dim>    &dof,
				       const FiniteElement<dim> &fe,
				       const Quadrature<dim>    &q,
				       const Boundary<dim>      &boundary,
				       dSMatrix &matrix,
				       const Function<dim>      &rhs,
				       dVector                  &rhs_vector,
				       const Function<dim>      *a = 0);
};





/**
 * Provide a collection of functions operating on matrices. These include
 * the application of boundary conditions to a linear system of equations
 * and others.
 *
 *
 * \subsection{Boundary conditions}
 *
 * The #apply_boundar_values# function inserts boundary conditions of
 * into a system of equations.  To actually do this you have to specify
 * a list of degree of freedom indices along with the value this degree of
 * freedom shall assume. To see how to get such a list, see the discussion
 * of the #VectorTools::interpolate_boundary_values# function.
 *
 * The inclusion into the assemblage process is as follows: when the matrix and
 * vectors are set up, a list of nodes subject to dirichlet bc is made and
 * matrix and vectors are changed accordingly. This is done by deleting all
 * entries in the matrix in the line of this degree of freedom, setting the
 * main diagonal entry to one and the right hand side element to the
 * boundary value at this node. This forces this node's value to be as specified.
 * To decouple the remaining linear system of equations and to make the system
 * symmetric again (at least if it was before), one Gauss elimination
 * step is performed with this line, by adding this (now almost empty) line to
 * all other lines which couple with the given degree of freedom and thus
 * eliminating all coupling between this degree of freedom and others. Now
 * also the column consists only of zeroes, apart from the main diagonal entry.
 *
 * It seems as if we had to make clear not to overwrite the lines of other
 * boundary nodes when doing the Gauss elimination step. However, since we
 * reset the right hand side when passing such a node, it is not a problem
 * to change the right hand side values of other boundary nodes not yet
 * processed. It would be a problem to change those entries of nodes already
 * processed, but since the matrix entry of the present column on the row
 * of an already processed node is zero, the Gauss step does not change
 * the right hand side. We need therefore not take special care of other
 * boundary nodes.
 * 
 * To make solving faster, we preset the solution vector with the right boundary
 * values. Since boundary nodes can never be hanging nodes, and since all other
 * entries of the solution vector are zero, we need not condense the solution
 * vector if the condensation process is done in-place. If done by copying
 * matrix and vectors to smaller ones, it would also be necessary to condense
 * the solution vector to preserve the preset boundary values.
 * 
 * It it not clear whether the deletion of coupling between the boundary degree
 * of freedom and other dofs really forces the corresponding entry in the
 * solution vector to have the right value when using iterative solvers,
 * since their search directions may contain components in the direction
 * of the boundary node. For this reason, we perform a very simple line
 * balancing by not setting the main diagonal entry to unity, but rather
 * to the value it had before deleting this line, or to the first nonzero
 * main diagonal entry if it is zero from a previous Gauss elimination
 * step. Of course we have to change
 * the right hand side appropriately. This is not a very good
 * strategy, but it at least should give the main diagonal entry a value
 * in the right order of dimension, which makes the solving process a bit
 * more stable. A refined algorithm would set the entry to the mean of the
 * other diagonal entries, but this seems to be too expensive.
 *
 * Because of the mentioned question, whether or not a preset solution value
 * which does not couple with other degrees of freedom remains its value or
 * not during solving iteratively, it may or may not be necessary to set
 * the correct value after solving again. This question is an open one as of
 * now and may be answered by future experience.
 *
 * 
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class MatrixTools : public MatrixCreator<dim> {
  public:
				     /**
				      * Apply dirichlet boundary conditions
				      * to the system matrix and vectors
				      * as described in the general
				      * documentation.
				      */
    static void apply_boundary_values (const map<int,double> &boundary_values,
				       dSMatrix              &matrix,
				       dVector               &solution,
				       dVector               &right_hand_side);

				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
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
			   const typename Triangulation<dim>::cell_iterator &) const;

				     /**
				      * Construct the cell matrix for this cell.
				      * If a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &) const;

				     /**
				      * Only construct the right hand side
				      * vector for this cell. You need to give
				      * a right hand side function to the
				      * constructor in order to call this
				      * function.
				      */
    virtual void assemble (dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &) const;
    
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
			   const typename Triangulation<dim>::cell_iterator &) const;

				     /**
				      * Construct the cell matrix for this cell.
				      * If a coefficient was given to the
				      * constructor, it is used.
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &) const;

				     /**
				      * Only construct the right hand side
				      * vector for this cell. You need to give
				      * a right hand side function to the
				      * constructor in order to call this
				      * function.
				      */
    virtual void assemble (dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const typename Triangulation<dim>::cell_iterator &) const;

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
