//----------------------------  matrices.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  matrices.h  ---------------------------
#ifndef __deal2__matrices_h
#define __deal2__matrices_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/thread_management.h>
#include <dofs/function_map.h>
#include <map>


// forward declarations
template <int dim> class Quadrature;

template<typename number> class Vector;
template<typename number> class FullMatrix;
template<typename number> class SparseMatrix;

template <typename number> class BlockSparseMatrix;
template <typename Number> class BlockVector;

template <int dim> class Mapping;
template <int dim> class DoFHandler;
template <int dim> class MGDoFHandler;
template <int dim> class FEValues;

#ifdef DEAL_II_USE_PETSC
namespace PETScWrappers
{
  class SparseMatrix;
  class Vector;
}
#endif


/**
 * Provide a class which assembles certain standard matrices for a
 * given triangulation, using a given finite element, a given mapping
 * and a quadrature formula.  All functions are static, so it is not
 * necessary to create an object of this type, though you may do so.
 *
 *
 * @sect3{Conventions for all functions}
 *
 * There exist two versions of each function. One with a Mapping
 * argument and one without. If a code uses a mapping different from
 * MappingQ1 the functions <em>with</em> mapping argument should
 * be used. Code that uses only MappingQ1 may also use the
 * functions <em>without</em> Mapping argument. Each of these
 * latter functions create a MappingQ1 object and just call the
 * respective functions with that object as mapping argument. The
 * functions without Mapping argument still exist to ensure
 * backward compatibility. Nevertheless it is advised to change the
 * user's codes to store a specific Mapping object and to use
 * the functions that take this @p Mapping object as argument. This
 * gives the possibility to easily extend the user codes to work also
 * on mappings of higher degree, this just by exchanging
 * MappingQ1 by, for example, a MappingQ or another
 * Mapping object of interest.
 *
 * All functions take a sparse matrix object to hold the matrix to be
 * created. The functions assume that the matrix is initialized with a
 * sparsity pattern (SparsityPattern) corresponding to the given degree
 * of freedom handler, i.e. the sparsity structure is already as needed.
 * You can do this by calling the DoFTools::make_sparsity_pattern()
 * function.
 *
 * Furthermore it is assumed that no relevant data is in the matrix. All
 * entries will be overwritten. Entries which are not needed by the matrix
 * are not touched and in special are not set to zero.
 * In all cases, the elements of the matrix to be assembled are simply
 * summed up from the contributions of each cell. Therefore you may want
 * to clear the matrix before assemblage.
 *
 * All created matrices are `raw': they are not condensed,
 * i.e. hanging nodes are not eliminated. The reason is that you may
 * want to add several matrices and could then condense afterwards
 * only once, instead of for every matrix. To actually do computations
 * with these matrices, you have to condense the matrix using the
 * ConstraintMatrix@p ::condense function; you also have to
 * condense the right hand side accordingly and distribute the
 * solution afterwards.
 *
 * If you want to use boundary conditions with the matrices generated
 * by the functions of this class, you have to use a function like
 * <tt>ProblemBase<>::apply_dirichlet_bc</tt> to matrix and right hand
 * side.
 *
 *
 * @sect3{Supported matrices}
 *
 * At present there are functions to create the following matrices:
 * <ul>
 * <li> @p create_mass_matrix: create the matrix with entries
 *   $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$ by numerical
 *   quadrature. Here, the $\phi_i$ are the basis functions of the
 *   finite element space given.
 *
 *   A coefficient may be given to evaluate
 *   $m_{ij} = \int_\Omega a(x) \phi_i(x) \phi_j(x) dx$ instead.
 *
 * <li> @p create_laplace_matrix: create the matrix with entries
 *   $m_{ij} = \int_\Omega \nabla\phi_i(x) \nabla\phi_j(x) dx$ by
 *   numerical quadrature.
 *
 *   Again, a coefficient may be given to evaluate
 *   $m_{ij} = \int_\Omega a(x) \nabla\phi_i(x) \phi_j(x) dx$ instead.
 * </ul>
 *
 * Make sure that the order of the Quadrature formula given to these
 * functions is sufficiently high to compute the matrices with the
 * required accuracy. For the choice of this quadrature rule you need
 * to take into account the polynomial degree of the FiniteElement
 * basis functions, the roughness of the coefficient @p a, as well as
 * the degree of the given @p Mapping.
 *
 * Note, that for system elements the mass matrix and the laplace
 * matrix is implemented such that each components couples only with
 * itself. I.e. there is no coupling of shape functions belonging to
 * different components.
 *
 * If the finite element for which the mass matrix or the laplace
 * matrix is to be built has more than one component, this function
 * accepts a single coefficient as well as a vector valued coefficient
 * function. For the latter case make sure that the number of
 * components coincides with the number of components of the system
 * finite element.
 *
 *
 * @sect3{Matrices on the boundary}
 *
 * The @p create_boundary_mass_matrix creates the matrix with entries
 * $m_{ij} = \int_{\Gamma} \phi_i \phi_j dx$, where $\Gamma$ is the
 * union of boundary parts with indicators contained in a
 * FunctioMap@p ::FunctionMap passed to the function (i.e. if
 * you want to set up the mass matrix for the parts of the boundary
 * with indicators zero and 2, you pass the function a map of
 * <tt>unsigned char</tt>s as parameter @p boundary_functions containing
 * the keys zero and 2). The $\phi_i$ are the basis functions which
 * have at least part of their support on $\Gamma$. The mapping
 * @p dof_to_boundary_mapping required by this function maps global
 * DoF numbers to a numbering of the degrees of freedom located on the
 * boundary, and can be obtained using the function
 * @p DoFTools::map_dof_to_boundary_indices.
 *
 * Since in most cases we are not interested in the pure mass matrix on the
 * boundary, but rather need it to compute the projection of a function to
 * the boundary, no function is provided to only create the matrix.
 *
 * This function needs to get passed a matrix object to hold the resulting sparse
 * matrix. This object is supposed to be initialized with a suitable sparsity
 * pattern, which can be created using the
 * DoFHandler@p ::make_boundary_sparsity_pattern function.
 *
 * The object describing the exact form of the boundary is obtained from the
 * triangulation object.
 *
 *
 * @sect3{Right hand sides}
 *
 * In many cases, you will not only want to build the matrix, but also
 * a right hand side, which will give a vector with
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$. For this purpose, each function
 * exists in two versions, one only building the matrix and one also
 * building the right hand side vector. If you want to create a right
 * hand side vector without creating a matrix, you can use the
 * VectorTools::create_right_hand_side() function. The use of the
 * latter may be useful if you want to create many right hand side
 * vectors.
 *
 * Creation of the right hand side is the same for all operators and
 * therefore for all of the functions below. It would be most
 * orthogonal to write one single function which builds up the right
 * hand side and not provide many functions doing the same
 * thing. However, this may result in a heavy performance penalty,
 * since then many values of a certain finite element have to be
 * computed twice, so it is more economical to implement it more than
 * once.
 *
 * All functions in this collection use the finite element given to
 * the DoFHandler object the last time that the degrees of
 * freedom were distributed on the triangulation.
 *
 * @author Wolfgang Bangerth, 1998, Ralf Hartmann, 2001
 */
class MatrixCreator
{
  public:
				     /**
				      * Assemble the mass matrix. If no 
				      * coefficient is given, it is assumed
				      * to be unity.
				      *
				      * If the library is configured
				      * to use multithreading, this
				      * function works in parallel.
				      *
				      * See the general doc of this class
				      * for more information.
				      */
    template <int dim, typename number>
    static void create_mass_matrix (const Mapping<dim>       &mapping,
				    const DoFHandler<dim>    &dof,
				    const Quadrature<dim>    &q,
				    SparseMatrix<number>     &matrix,
				    const Function<dim> * const a = 0);

				     /**
				      * Calls the create_mass_matrix()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim, typename number>
    static void create_mass_matrix (const DoFHandler<dim>    &dof,
				    const Quadrature<dim>    &q,
				    SparseMatrix<number>     &matrix,
				    const Function<dim> * const a = 0);

    				     /**
				      * Assemble the mass matrix and a
				      * right hand side vector. If no
				      * coefficient is given, it is
				      * assumed to be unity.
				      *
				      * If the library is configured
				      * to use multithreading, this
				      * function works in parallel.
				      *
				      * See the general doc of this
				      * class for more information.
				      */
    template <int dim, typename number>
    static void create_mass_matrix (const Mapping<dim>       &mapping,
				    const DoFHandler<dim>    &dof,
				    const Quadrature<dim>    &q,
				    SparseMatrix<number>     &matrix,
				    const Function<dim>      &rhs,
				    Vector<double>           &rhs_vector,
				    const Function<dim> * const a = 0);

				     /**
				      * Calls the create_mass_matrix()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim, typename number>
    static void create_mass_matrix (const DoFHandler<dim>    &dof,
				    const Quadrature<dim>    &q,
				    SparseMatrix<number>     &matrix,
				    const Function<dim>      &rhs,
				    Vector<double>           &rhs_vector,
				    const Function<dim> * const a = 0);
    
				     /**
				      * Assemble the mass matrix and a
				      * right hand side vector along
				      * the boundary.  If no
				      * coefficient is given, it is
				      * assumed to be constant one.
				      *
				      * The matrix is assumed to
				      * already be initialized with a
				      * suiting sparsity pattern (the
				      * DoFHandler provides an
				      * appropriate function).
				      *
				      * If the library is configured
				      * to use multithreading, this
				      * function works in parallel.
				      *
				      * See the general doc of this
				      * class for more information.
				      */
    template <int dim>
    static
    void create_boundary_mass_matrix (const Mapping<dim>       &mapping,
				      const DoFHandler<dim>    &dof,
				      const Quadrature<dim-1>  &q,
				      SparseMatrix<double>     &matrix,
				      const typename FunctionMap<dim>::type &boundary_functions,
				      Vector<double>           &rhs_vector,
				      std::vector<unsigned int>&dof_to_boundary_mapping,
				      const Function<dim> * const a = 0);

				     /**
				      * Same function, but for 1d.
				      */
    static
    void create_boundary_mass_matrix (const Mapping<1>       &mapping,
				      const DoFHandler<1>    &dof,
				      const Quadrature<0>    &q,
				      SparseMatrix<double>   &matrix,
				      const FunctionMap<1>::type &boundary_functions,
				      Vector<double>         &rhs_vector,
				      std::vector<unsigned int>&dof_to_boundary_mapping,
				      const Function<1> * const a = 0);


				     /**
				      * Calls the
				      * create_boundary_mass_matrix()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim>
    static
    void create_boundary_mass_matrix (const DoFHandler<dim>    &dof,
				      const Quadrature<dim-1>  &q,
				      SparseMatrix<double>     &matrix,
				      const typename FunctionMap<dim>::type        &boundary_functions,
				      Vector<double>           &rhs_vector,
				      std::vector<unsigned int>&dof_to_boundary_mapping,
				      const Function<dim> * const a = 0);

				     /**
				      * Assemble the Laplace
				      * matrix. If no coefficient is
				      * given, it is assumed to be
				      * constant one.
				      * 
				      * If the library is configured
				      * to use multithreading, this
				      * function works in parallel.
				      *
				      * See the general doc of this
				      * class for more information.
				      */
    template <int dim>
    static void create_laplace_matrix (const Mapping<dim>       &mapping,
				       const DoFHandler<dim>    &dof,
				       const Quadrature<dim>    &q,
				       SparseMatrix<double>     &matrix,
				       const Function<dim> * const a = 0);
    
				     /**
				      * Calls the
				      * create_laplace_matrix()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim>
    static void create_laplace_matrix (const DoFHandler<dim>    &dof,
				       const Quadrature<dim>    &q,
				       SparseMatrix<double>     &matrix,
				       const Function<dim> * const a = 0);

				     /**
				      * Assemble the Laplace matrix
				      * and a right hand side
				      * vector. If no coefficient is
				      * given, it is assumed to be
				      * constant one.
				      * 
				      * If the library is configured
				      * to use multithreading, this
				      * function works in parallel.
				      *
				      * See the general doc of this
				      * class for more information.
				      */
    template <int dim>
    static void create_laplace_matrix (const Mapping<dim>       &mapping,
				       const DoFHandler<dim>    &dof,
				       const Quadrature<dim>    &q,
				       SparseMatrix<double>     &matrix,
				       const Function<dim>      &rhs,
				       Vector<double>           &rhs_vector,
				       const Function<dim> * const a = 0);

				     /**
				      * Calls the
				      * create_laplace_matrix()
				      * function, see above, with
				      * <tt>mapping=MappingQ1@<dim@>()</tt>.
				      */
    template <int dim>
    static void create_laplace_matrix (const DoFHandler<dim>    &dof,
				       const Quadrature<dim>    &q,
				       SparseMatrix<double>     &matrix,
				       const Function<dim>      &rhs,
				       Vector<double>           &rhs_vector,
				       const Function<dim> * const a = 0);

				     /**
				      * Exception
				      */
    DeclException0 (ExcComponentMismatch);    
  private:
				     /**
				      * Convenience abbreviation for
				      * pairs of DoF handler cell
				      * iterators. This type works
				      * just like a
				      * <tt>std::pair<iterator,iterator></tt>
				      * but is templatized on the
				      * space dimension.
				      */
    template <int dim>
    struct IteratorRange 
    {
					 /**
					  * Typedef for the iterator type.
					  */
	typedef typename DoFHandler<dim>::active_cell_iterator active_cell_iterator;

					 /**
					  * Abbreviation for a pair of
					  * iterators.
					  */
	typedef std::pair<active_cell_iterator,active_cell_iterator> iterator_pair;
	
					 /**
					  * Constructor. Initialize
					  * the two values by the
					  * given values.
					  */
	IteratorRange (const active_cell_iterator &first,
		       const active_cell_iterator &second);

					 /**
					  * Constructor taking a pair
					  * of values for
					  * initialization.
					  */
	IteratorRange (const iterator_pair &ip);
	
					 /**
					  * Pair of iterators denoting
					  * a half-open range.
					  */
	active_cell_iterator first, second;
    };
    
    
				     /**
				      * Version of the same function
				      * (without suffix @p _1) with
				      * the same argument list that
				      * operates only on an interval
				      * of iterators. Used for
				      * parallelization. The mutex is
				      * used to synchronise access to
				      * the matrix.
				      */
    template <int dim, typename number>
    static
    void create_mass_matrix_1 (const Mapping<dim>       &mapping,
			       const DoFHandler<dim>    &dof,
			       const Quadrature<dim>    &q,
			       SparseMatrix<number>     &matrix,
			       const Function<dim> * const a,
			       const IteratorRange<dim>  range,
			       Threads::ThreadMutex     &mutex);

				     /**
				      * Version of the same function
				      * (without suffix @p _2) with
				      * the same argument list that
				      * operates only on an interval
				      * of iterators. Used for
				      * parallelization. The mutex is
				      * used to synchronise access to
				      * the matrix.
				      */
    template <int dim, typename number>
    static
    void create_mass_matrix_2 (const Mapping<dim>       &mapping,
			       const DoFHandler<dim>    &dof,
			       const Quadrature<dim>    &q,
			       SparseMatrix<number>     &matrix,
			       const Function<dim>      &rhs,
			       Vector<double>           &rhs_vector,
			       const Function<dim> * const a,
			       const IteratorRange<dim>  range,
			       Threads::ThreadMutex     &mutex);

				     /**
				      * Version of the same function
				      * (without suffix @p _1) with
				      * the same argument list that
				      * operates only on an interval
				      * of iterators. Used for
				      * parallelization. The mutex is
				      * used to synchronise access to
				      * the matrix.
				      */
    template <int dim>
    static
    void create_laplace_matrix_1 (const Mapping<dim>       &mapping,
				  const DoFHandler<dim>    &dof,
				  const Quadrature<dim>    &q,
				  SparseMatrix<double>     &matrix,
				  const Function<dim> * const a,
				  const IteratorRange<dim>  range,
				  Threads::ThreadMutex     &mutex);

				     /**
				      * Version of the same function
				      * (without suffix @p _2) with
				      * the same argument list that
				      * operates only on an interval
				      * of iterators. Used for
				      * parallelization. The mutex is
				      * used to synchronise access to
				      * the matrix.
				      */
    template <int dim>
    static
    void create_laplace_matrix_2 (const Mapping<dim>       &mapping,
				  const DoFHandler<dim>    &dof,
				  const Quadrature<dim>    &q,
				  SparseMatrix<double>     &matrix,
				  const Function<dim>      &rhs,
				  Vector<double>           &rhs_vector,
				  const Function<dim> * const a,
				  const IteratorRange<dim>  range,
				  Threads::ThreadMutex     &mutex);

				     /**
				      * Version of the same function
				      * (without suffix @p _1) with
				      * the same argument list that
				      * operates only on an interval
				      * of iterators. Used for
				      * parallelization. The mutex is
				      * used to synchronise access to
				      * the matrix.
				      */
    template <int dim>
    static
    void create_boundary_mass_matrix_1 (const Mapping<dim>       &mapping,
					const DoFHandler<dim>    &dof,
					const Quadrature<dim-1>  &q,
					SparseMatrix<double>     &matrix,
					const typename FunctionMap<dim>::type        &boundary_functions,
					Vector<double>           &rhs_vector,
					std::vector<unsigned int>&dof_to_boundary_mapping,
					const Function<dim> * const a,
					const IteratorRange<dim>  range,
					Threads::ThreadMutex     &mutex);
};



/**
 * Provide a collection of functions operating on matrices. These include
 * the application of boundary conditions to a linear system of equations
 * and others.
 *
 *
 * @sect3{Boundary conditions}
 *
 * The apply_boundary_values() function inserts boundary conditions
 * into a system of equations.  To actually do this you have to
 * specify a list of degree of freedom indices along with the values
 * these degrees of freedom shall assume. To see how to get such a
 * list, see the discussion of the
 * VectorTools::interpolate_boundary_values function.
 *
 * There are two ways to incorporate fixed degrees of freedom such as boundary
 * nodes into a linear system, as discussed below.
 * 
 *
 * @sect3{Global elimination}
 *
 * In the first method, we first assemble the global linear system without
 * respect for fixed degrees of freedom, and in a second step eliminate them
 * again from the linear system. The inclusion into the assembly process is as
 * follows: when the matrix and vectors are set up, a list of nodes subject to
 * dirichlet bc is made and matrix and vectors are modified accordingly. This
 * is done by deleting all entries in the matrix in the line of this degree of
 * freedom, setting the main diagonal entry to a suitable positive value and
 * the right hand side element to a value so that the solution of the linear
 * system will have the boundary value at this node. To decouple the remaining
 * linear system of equations and to make the system symmetric again (at least
 * if it was before), one Gauss elimination step is performed with this line,
 * by adding this (now almost empty) line to all other lines which couple with
 * the given degree of freedom and thus eliminating all coupling between this
 * degree of freedom and others. Now the respective column also consists only
 * of zeroes, apart from the main diagonal entry. Alternatively, the functions
 * in this class take a boolean parameter that allows to omit this last step,
 * if symmetry of the resulting linear system is not required. Note that
 * usually even CG can cope with a non-symmetric linear system with this
 * particular structure.
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
 * (@p elininate_columns) that specifies whether the sparsity pattern
 * is symmetric; if so, then the column is eliminated and the right
 * hand side is also modified accordingly. If not, then only the row
 * is deleted and the column is not touched at all, and all right hand
 * side values apart from the one corresponding to the present row
 * remain unchanged.
 *
 * If the sparsity pattern for your matrix is non-symmetric, you must
 * set the value of this parameter to @p false in any case, since then
 * we can't eliminate the column without searching all rows, which
 * would be too expensive (if @p N be the number of rows, and @p m the
 * number of nonzero elements per row, then eliminating one column is
 * an <tt>O(N*log(m))</tt> operation, since searching in each row takes
 * <tt>log(m)</tt> operations). If your sparsity pattern is symmetric, but
 * your matrix is not, then you might specify @p false as well. If your
 * sparsity pattern and matrix are both symmetric, you might want to
 * specify @p true (the complexity of eliminating one row is then
 * <tt>O(m*log(m))</tt>, since we only have to search @p m rows for the
 * respective element of the column). Given the fact that @p m is
 * roughly constant, irrespective of the discretization, and that the
 * number of boundary nodes is <tt>sqrt(N)</tt> in 2d, the algorithm for
 * symmetric sparsity patterns is <tt>O(sqrt(N)*m*log(m))</tt>, while it
 * would be <tt>O(N*sqrt(N)*log(m))</tt> for the general case; the latter
 * is too expensive to be performed.
 *
 * It seems as if we had to make clear not to overwrite the lines of
 * other boundary nodes when doing the Gauss elimination
 * step. However, since we reset the right hand side when passing such
 * a node, it is not a problem to change the right hand side values of
 * other boundary nodes not yet processed. It would be a problem to
 * change those entries of nodes already processed, but since the
 * matrix entry of the present column on the row of an already
 * processed node is zero, the Gauss step does not change the right
 * hand side. We need therefore not take special care of other
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
 * In some cases, it might be interesting to solve several times with
 * the same matrix, but for different right hand sides or boundary
 * values. However, since the modification for boundary values of the
 * right hand side vector depends on the original matrix, this is not
 * possible without storing the original matrix somewhere and applying
 * the @p apply_boundary_conditions function to a copy of it each
 * time we want to solve. In that case, you can use the
 * FilteredMatrix class in the @p LAC sublibrary. There you can
 * also find a formal (mathematical) description of the process of
 * modifying the matrix and right hand side vectors for boundary
 * values.
 *
 *
 * @sect3{Local elimination}
 *
 * The second way of handling boundary values is to modify the local matrix
 * and vector contributions appropriately before transferring them into the
 * global sparse matrix and vector. This is what local_apply_boundary_values()
 * does. The advantage is that we save the call to the apply_boundary_values
 * function (which is expensive because it has to work on sparse data
 * structures). On the other hand, the local_apply_boundary_values() function
 * is called many times, even if we only have a very small number of fixed
 * boundary nodes.
 *
 * However, since we do not have access to the data structures of some sparse
 * matrix formats (e.g. the PETSc matrix classes), this may be the only way to
 * get rid of boundary nodes for these matrix formats. In general, this
 * function should not be a loss in efficiency over the global one.
 * 
 * @author Wolfgang Bangerth, 1998, 2000, 2004
 */
class MatrixTools : public MatrixCreator
{
  public:
				     /**
				      * Apply dirichlet boundary conditions
				      * to the system matrix and vectors
				      * as described in the general
				      * documentation.
				      *
				      * For a replacement function,
				      * see the documentation of the
				      * FilteredMatrix class in
				      * the @p LAC sublibrary, or use the
				      * local_apply_boundary_values()
				      * function..
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
				      *
				      * For a replacement function, see the
				      * documentation of the
				      * FilteredMatrix class in the
				      * @p LAC sublibrary, or use the
				      * local_apply_boundary_values()
				      * function.
				      */
    static void
    apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
			   BlockSparseMatrix<double>           &matrix,
			   BlockVector<double>                 &solution,
			   BlockVector<double>                 &right_hand_side,
			   const bool           eliminate_columns = true);

				     /**
				      * Apply dirichlet boundary conditions to
				      * the system matrix and vectors as
				      * described in the general
				      * documentation. This function works on
				      * the classes that are used to wrap
				      * PETSc objects.
				      *
				      * Note that this function is not very
				      * efficient: it needs to alternatingly
				      * read and write into the matrix, a
				      * situation that PETSc does not handle
				      * too well. In addition, we only get rid
				      * of rows corresponding to boundary
				      * nodes, but the corresponding case of
				      * deleting the respective columns
				      * (i.e. if @p eliminate_columns is @p
				      * true) is not presently implemented,
				      * and probably will never because it is
				      * too expensive without direct access to
				      * the PETSc data structures. A third
				      * reason against this function is that
				      * it doesn't handle the case where the
				      * matrix is distributed across an MPI
				      * system.
				      *
				      * In order to still be able to eliminate
				      * boundary values, it is better to get
				      * rid of them before the local matrices
				      * and vectors are distributed to the
				      * global ones, because then we don't
				      * have to mess with the sparse data
				      * structures. The
				      * local_apply_boundary_values() function
				      * does that, and is recommended for use
				      * instead of the global one for PETSc
				      * matrices and vectors.
				      */
#ifdef DEAL_II_USE_PETSC
    static void
    apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
			   PETScWrappers::SparseMatrix  &matrix,
			   PETScWrappers::Vector  &solution,
			   PETScWrappers::Vector  &right_hand_side,
			   const bool             eliminate_columns = true);
#endif

                                     /**
                                      * Rather than applying boundary values
                                      * to the global matrix and vector after
                                      * assembly, this function does so @em
                                      * before assembling, by modifying the
                                      * local matrix and vector
                                      * contributions. If you call this
                                      * function on all local contributions,
                                      * the resulting matrix will have the
                                      * same entries, and the final call to
                                      * apply_boundary_values() on the global
                                      * system will not be necessary.
                                      *
                                      * Since this function does not have to
                                      * work on the complicated data
                                      * structures of sparse matrices, it is
                                      * relatively cheap. It may therefore be
                                      * a win if you have many fixed degrees
                                      * of freedom (e.g. boundary nodes), or
                                      * if access to the sparse matrix is
                                      * expensive (e.g. for block sparse
                                      * matrices, or for PETSc matrices).
                                      */
    static void
    local_apply_boundary_values (const std::map<unsigned int,double> &boundary_values,
                                 const std::vector<unsigned int> &local_dof_indices,
                                 FullMatrix<double> &local_matrix,
                                 Vector<double>     &local_rhs,
                                 const bool          eliminate_columns);
    
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
