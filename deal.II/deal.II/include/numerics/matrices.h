//----------------------------  matrices.h  ---------------------------
//    Version: $Name$
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  matrices.h  ---------------------------
#ifndef __deal2__matrices_h
#define __deal2__matrices_h


#include <base/exceptions.h>
#include <map>


// forward declarations
template <int dim> class Quadrature;

template<typename number> class Vector;
template<typename number> class FullMatrix;
template<typename number> class SparseMatrix;

template <typename number> class BlockSparseMatrix;
template <typename Number> class BlockVector;

template <int dim> class DoFHandler;
template <int dim> class MGDoFHandler;
template <int dim> class FEValues;



/**
 * Provide a class which assembles certain standard matrices for a given
 * triangulation, using a given finite element and a quadrature formula.
 * All functions are static, so it is not necessary to create an object
 * of this type, though you may do so.
 *
 *
 * @sect3{Conventions for all functions}
 *
 * All functions take a sparse matrix object to hold the matrix to be
 * created. The functions assume that the matrix is initialized with a
 * sparsity pattern (@ref{SparsityPattern}) corresponding to the given degree
 * of freedom handler, i.e. the sparsity structure is already as needed.
 * You can do this by calling the @ref{DoFHandler}@p{<dim>::make_sparsity_pattern}
 * function.
 *
 * Furthermore it is assumed that no relevant data is in the matrix. All
 * entries will be overwritten. Entries which are not needed by the matrix
 * (and were thus added 'by hand' after @p{make_sparsity_pattern} was called)
 * are not touched and in special are not set to zero, so you have to care
 * yourself about that if you really need these entries.
 *
 *
 * @sect3{Supported matrices}
 *
 * At present there are functions to create the following matrices:
 * @begin{itemize}
 * @item @p{create_mass_matrix}: create the matrix with entries
 *   $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$. Here, the $\phi_i$
 *   are the basis functions of the finite element space given.
 *   This function uses the @ref{MassMatrix} class.
 *
 *   Two ways to create this matrix are offered. The first one uses
 *   numerical quadrature and the @ref{MassMatrix} class. In this case,
 *   a coefficient may be given to evaluate
 *   $m_{ij} = \int_\Omega a(x) \phi_i(x) \phi_j(x) dx$ instead.
 *   This way of setting up the mass matrix is quite general, but has
 *   some drawbacks, see the documentation of the @ref{MassMatrix} class.
 *
 *   The other way uses exact integration, as offered by the finite
 *   element class used. This way you can avoid quadrature errors and
 *   the assemblage is much faster. However, no coefficient can be
 *   given.
 *
 *   Note that the effect of the two ways of setting up the mass
 *   matrix is not the same if you use finite elements which are
 *   composed of several subelements. In this case, using the
 *   quadrature free way (without coefficient) results in a matrix
 *   which does not couple the subelements, as described in the
 *   @ref{FESystem}@p{::get_local_mass_matrix} documentation, while the
 *   way using quadrature sets up the full matrix, i.e. with the
 *   cross coupling of shape functions belonging to different subelements.
 *
 *   If the finite element for which the mass matrix is to be built
 *   has more than one component, the resulting matrix will not couple
 *   the different components. It will furthermore accept a single
 *   coefficient through the @ref{Function} parameter for all
 *   components. If you want different coefficients for the different
 *   parameters, you have to assemble the matrix yourself, sorry; the
 *   implementation of the function will serve as a good starting
 *   point, though. (You may also modify the implementation to accept
 *   vector-valued functions and send this implementation to us -- we
 *   will then include this implementation into the library.)
 *
 * @item @p{create_laplace_matrix}: there are two versions of this; the
 *   one which takes the @ref{Function} object creates
 *   $a_{ij} = \int_\Omega a(x) \nabla\phi_i(x) \nabla\phi_j(x) dx$,
 *   $a$ being the given function, while the other one assumes that
 *   $a=1$ which enables some optimizations. In fact the two versions
 *   are in one function, the coefficient being given as a defaulted
 *   argument, which is a pointer to a function and defaults to zero.
 *   This function uses the @ref{LaplaceMatrix} class.
 *
 *   If the finite element in use presently has more than only one
 *   component, this function may not be overly useful. It assembles a
 *   Laplace matrix block for each component (with the same
 *   coefficient for each component). These blocks are not coupled.  *
 *   @end{itemize}
 *
 * All created matrices are `raw': they are not condensed, i.e. hanging
 * nodes are not eliminated. The reason is that you may want to add
 * several matrices and could then condense afterwards only once,
 * instead of for every matrix. To actually do computations with these
 * matrices, you have to condense the matrix using the
 * @ref{ConstraintMatrix}@p{::condense} function; you also have to condense the
 * right hand side accordingly and distribute the solution afterwards.
 *
 * In all cases, the elements of the matrix to be assembled are simply
 * summed up from the contributions of each cell. Therefore you may want
 * to clear the matrix before assemblage.
 *
 * If you want to use boundary conditions, you have to use a function
 * like @p{ProblemBase<>::apply_dirichlet_bc} to matrix and right hand
 * side.
 *
 *
 * @sect3{Matrices on the boundary}
 *
 * The @p{create_boundary_mass_matrix} creates the matrix with entries
 * $m_{ij} = \int_{\Gamma} \phi_i \phi_j dx$, where $\Gamma$ is the union
 * of boundary parts with indicators contained in a set passed to the
 * function (i.e. if you want to set up the mass matrix for the parts of
 * the boundary with indicators zero and 2, you pass the function a set
 * of @p{unsigned char}s as parameter @p{boundary_parts} containing the elements
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
 * @ref{DoFHandler}@p{::make_boundary_sparsity_pattern} function.
 *
 * The object describing the exact form of the boundary is obtained from the
 * triangulation object.
 *
 * @sect3{Right hand sides}
 *
 * In many cases, you will not only want to build the matrix, but also
 * a right hand side, which will give a vector with
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$. For this purpose, each function
 * exists in two versions, one only building the matrix and one also
 * building the right hand side vector. (The @p{create_mass_matrix} function
 * which does not use quadrature does not offer a version to evaluate a right
 * hand side also, since this needs quadrature anyway. Take look at the
 * @ref{VectorTools} class to find a function to set up a right hand side vector
 * only.)
 *
 * Creation of the right hand side
 * is the same for all operators and therefore for all of the functions
 * below. It would be most orthogonal to write one single function which
 * builds up the right hand side and not provide many functions doing
 * the same thing. However, this may result in a heavy performance
 * penalty, since then many values of a certain finite element have to
 * be computed twice, so it is more economical to  implement it more than
 * once. If you only want to create a right hand side as above, there is
 * a function in the @p{VectorCreator} class. The use of the latter may be
 * useful if you want to create many right hand side vectors.
 *
 *
 * All functions in this collection use the finite elemen given to the
 * @ref{DoFHandler} object the last time that the degrees of freedom were
 * distributed on the triangulation.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class MatrixCreator
{
  public:
				     /**
				      *	Declare a data type which denotes a
				      *	mapping between a boundary indicator
				      *	and the function denoting the boundary
				      *	values on this part of the boundary.
				      *	Only one boundary function may be given
				      *	for each boundary indicator, which is
				      *	guaranteed by the @p{map} data type.
				      *	
				      *	See the general documentation of this
				      *	class for more detail.
				      */
    typedef typename std::map<unsigned char,const Function<dim>*> FunctionMap;

				     /**
				      * Assemble the mass matrix. If no 
				      * coefficient is given, it is assumed
				      * to be unity.
				      * 
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_mass_matrix (const DoFHandler<dim>    &dof,
				    const Quadrature<dim>    &q,
				    SparseMatrix<double>     &matrix,
				    const Function<dim>      *a = 0);

    				     /**
				      * Assemble the mass matrix and a right
				      * hand side vector. If no 
				      * coefficient is given, it is assumed
				      * to be unity.
				      *
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_mass_matrix (const DoFHandler<dim>    &dof,
				    const Quadrature<dim>    &q,
				    SparseMatrix<double>     &matrix,
				    const Function<dim>      &rhs,
				    Vector<double>           &rhs_vector,
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
				      * pattern (the @ref{DoFHandler} provides an
				      * appropriate function).
				      *
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_boundary_mass_matrix (const DoFHandler<dim>    &dof,
					     const Quadrature<dim-1>  &q,
					     SparseMatrix<double>     &matrix,
					     const FunctionMap        &rhs,
					     Vector<double>           &rhs_vector,
					     std::vector<unsigned int>&vec_to_dof_mapping,
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
				       const Quadrature<dim>    &q,
				       SparseMatrix<double>     &matrix,
				       const Function<dim>      *a = 0);

				     /**
				      * Generate Laplace matrix for a given level.
				      * 
				      * See the general doc of this class
				      * for more information.
				      */
    static void create_level_laplace_matrix (unsigned int             level,
					     const MGDoFHandler<dim>& dof,
					     const Quadrature<dim>&   q,
					     SparseMatrix<float>&     matrix,
					     const Function<dim>*     a = 0);

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
				       const Quadrature<dim>    &q,
				       SparseMatrix<double>     &matrix,
				       const Function<dim>      &rhs,
				       Vector<double>           &rhs_vector,
				       const Function<dim>      *a = 0);

				     /**
				      * Exception
				      */
    DeclException0 (ExcComponentMismatch);
};



/**
 * Provide a collection of functions operating on matrices. These include
 * the application of boundary conditions to a linear system of equations
 * and others.
 *
 *
 * @sect3{Boundary conditions}
 *
 * The @p{apply_boundary_values} function inserts boundary conditions of
 * into a system of equations.  To actually do this you have to specify
 * a list of degree of freedom indices along with the values these degrees of
 * freedom shall assume. To see how to get such a list, see the discussion
 * of the @ref{VectorTools}@p{::interpolate_boundary_values} function.
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
 * Finding which rows contain an entry in the column for which we are
 * presently performing a Gauss elimination step is either difficult
 * or very simple, depending on the circumstances. If the sparsity
 * pattern is symmetric (whether the matrix is symmetric is irrelevant
 * here), then we can infer the rows which have a nonzero entry in the
 * present column by looking at which columns in the present row are
 * nonempty. In this case, we only need to look into a fixed number of
 * rows and need not search all rows. On the other hand, if the
 * sparsity pattern is nonsymmetric, then we need to use an iterative
 * solver which can handle nonsymmetric matrices in any case, so there
 * may be no need to do the Gauss elimination anyway. In fact, this is
 * the way the function works: it takes a parameter
 * (@p{elininate_columns}) that specifies whether the sparsity pattern
 * is symmetric; if so, then the column is eliminated and the right
 * hand side is also modified accordingly. If not, then only the row
 * is deleted and the column is not touched at all, and all right hand
 * side values apart from the one corresponding to the present row
 * remain unchanged.
 *
 * If the sparsity pattern for your matrix is non-symmetric, you must
 * set the value of this parameter to @p{false} in any case, since then
 * we can't eliminate the column without searching all rows, which
 * would be too expensive (if @p{N} be the number of rows, and @p{m} the
 * number of nonzero elements per row, then eliminating one column is
 * an @p{O(N*log(m))} operation, since searching in each row takes
 * @p{log(m)} operations). If your sparsity pattern is symmetric, but
 * your matrix is not, then you might specify @p{false} as well. If your
 * sparsity pattern and matrix are both symmetric, you might want to
 * specify @p{true} (the complexity of eliminating one row is then
 * @p{O(m*log(m))}, since we only have to search @p{m} rows for the
 * respective element of the column). Given the fact that @p{m} is
 * roughly constant, irrespective of the discretization, and that the
 * number of boundary nodes is @p{sqrt(N)} in 2d, the algorithm for
 * symmetric sparsity patterns is @p{O(sqrt(N)*m*log(m))}, while it
 * would be @p{O(N*sqrt(N)*log(m))} for the general case; the latter
 * is too expensive to be performed.
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
 * To make solving faster, we preset the solution vector with the
 * right boundary values. It it not clear whether the deletion of
 * coupling between the boundary degree of freedom and other dofs
 * really forces the corresponding entry in the solution vector to
 * have the right value when using iterative solvers, since their
 * search directions may contain components in the direction of the
 * boundary node. For this reason, we perform a very simple line
 * balancing by not setting the main diagonal entry to unity, but
 * rather to the value it had before deleting this line, or to the
 * first nonzero main diagonal entry if it is zero for some reason.
 * Of course we have to change the right hand side appropriately. This
 * is not a very good strategy, but it at least should give the main
 * diagonal entry a value in the right order of dimension, which makes
 * the solvution process a bit more stable. A refined algorithm would
 * set the entry to the mean of the other diagonal entries, but this
 * seems to be too expensive.
 *
 * 
 * @author Wolfgang Bangerth, 1998, 2000
 */
template <int dim>
class MatrixTools : public MatrixCreator<dim>
{
  public:
				     /**
				      * Apply dirichlet boundary conditions
				      * to the system matrix and vectors
				      * as described in the general
				      * documentation.
				      */
    template <typename number>
    static void
    apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
			   SparseMatrix<number>  &matrix,
			   Vector<number>        &solution,
			   Vector<number>        &right_hand_side,
			   const bool             eliminate_columns = true);

				     /**
				      * Apply dirichlet boundary
				      * conditions to the system
				      * matrix and vectors as
				      * described in the general
				      * documentation. This function
				      * works for block sparse
				      * matrices and block vectors
				      */
    static void
    apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
			   BlockSparseMatrix<double>           &matrix,
			   BlockVector<double>                 &solution,
			   BlockVector<double>                 &right_hand_side,
			   const bool           eliminate_columns = true);
    
				     /**
				      * Exception
				      */
    DeclException2 (ExcDimensionsDontMatch,
		    int, int,
		    << "The dimensions " << arg1 << " and " << arg2
		    << " don't match.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotBlockSquare);
				     /**
				      * Exception
				      */
    DeclException0 (ExcBlocksDontMatch);
};



#endif
