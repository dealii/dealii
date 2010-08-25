//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__constraint_matrix_h
#define __deal2__constraint_matrix_h

#include <base/config.h>
#include <base/exceptions.h>
#include <base/index_set.h>
#include <base/subscriptor.h>
#include <base/template_constraints.h>

#include <lac/vector.h>

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <complex>

#include <boost/scoped_ptr.hpp>


DEAL_II_NAMESPACE_OPEN

template<int dim, class T> class Table;
template <typename> class FullMatrix;
class SparsityPattern;
class CompressedSparsityPattern;
class CompressedSetSparsityPattern;
class CompressedSimpleSparsityPattern;
class BlockSparsityPattern;
class BlockCompressedSparsityPattern;
class BlockCompressedSetSparsityPattern;
class BlockCompressedSimpleSparsityPattern;
template <typename number> class SparseMatrix;
template <typename number> class BlockSparseMatrix;
class BlockIndices;
template <typename MatrixType> struct IsBlockMatrix;

namespace internals
{
  struct GlobalRowsFromLocal;
}


/**
 * This class implements dealing with linear (possibly inhomogeneous)
 * constraints on degrees of freedom. In particular, it handles constraints of
 * the form $x_{i_1} = \sum_{j=2}^M a_{i_j} x_{i_j} + b_i$. In the context of
 * adaptive finite elements, such constraints appear most frequently as
 * "hanging nodes" and for implementing Dirichlet boundary conditions in
 * strong form. The class is meant to deal with a limited number of
 * constraints relative to the total number of degrees of freedom, for example
 * a few per cent up to maybe 30 per cent; and with a linear combination of
 * $M$ other degrees of freedom where $M$ is also relatively small (no larger
 * than at most around the average number of entries per row of a linear
 * system). It is <em>not</em> meant to describe full rank linear systems.
 *
 * The algorithms used in the implementation of this class are described in
 * some detail in the @ref hp_paper "hp paper".
 *
 *
 * <h3>Using the %ConstraintMatrix for hanging nodes</h3>
 *
 * For example, when using Q1 and Q2 elements (i.e. using
 * FE_Q&lt;dim,spacedim&gt;(1) and FE_Q&lt;dim,spacedim&gt;(2)) on the two
 * marked cells of the mesh
 *
 * @image html hp-refinement-simple.png
 *
 * there are three constraints: first $x_2=\frac 12 x_0 + \frac 12 x_1$,
 * then $x_4=\frac 14 x_0 + \frac 34 x_1$, and finally the identity
 * $x_3=x_1$. All three constraints fit the form given above. Similar
 * constraints occur as hanging nodes even if all used finite elements are
 * identical. While they are most frequent for hanging nodes, constraints of
 * the given form appear also in other contexts, see for example the
 * application the step-11 tutorial program.
 *
 * Homogenous constraints of this form also arise in the context of vector-valued
 * fields, for example if one wants to enforce boundary conditions of the form
 * $\vec{v}\cdot\vec{n}=0$. For example, the
 * VectorTools::compute_no_normal_flux_constraints function computes
 * such constraints.
 *
 *
 * <h3>Using the %ConstraintMatrix for Dirichlet boundary conditions</h3>
 *
 * The ConstraintMatrix provides an alternative for implementinging
 * Dirichlet boundary conditions (the standard way that is extensively
 * discussed in the tutorial programs is to use
 * MatrixTools::apply_boundary_values). The general principle of Dirichlet
 * conditions are algebraic constraints of the form $x_{i} = b_i$, which
 * fits into the form described above.
 *
 *
 * <h3>Description of constraints</h3>
 *
 * Each "line" in objects of this class corresponds to one constrained degree
 * of freedom, with the number of the line being $i_1$, and the entries in
 * this line being pairs $(i_j,a_{i_j})$. Note that the constraints are linear
 * in the $x_i$, and that there might be a constant (non-homogeneous) term in
 * the constraint. This is exactly the form we need for hanging node
 * constraints, where we need to constrain one degree of freedom in terms of
 * others. There are other conditions of this form possible, for example for
 * implementing mean value conditions as is done in the step-11
 * tutorial program. The name of the class stems from the fact that these
 * constraints can be represented in matrix form as $X x = b$, and this object
 * then describes the matrix $X$ (as well as, incidentally, the vector $b$ --
 * originally, the ConstraintMatrix class was only meant to handle homogenous
 * constraints where $b=0$, thus the name). The most frequent way to
 * create/fill objects of this type is using the
 * DoFTools::make_hanging_node_constraints() function. The use of these
 * objects is first explained in step-6.
 *
 * Matrices of the present type are organized in lines (rows), but only those
 * lines are stored where constraints are present. New constraints are added
 * by adding new lines using the add_line() function, and then populating it
 * using the add_entry() function to a given line, or add_entries() to add
 * more than one entry at a time. The right hand side element, if nonzero, can
 * be set using the set_inhomogeneity() function. After all constraints have
 * been added, you need to call close(), which compresses the storage format
 * and sorts the entries.
 *
 * <h3>Eliminating constraints</h3>
 *
 * Constraint matrices are used to handle hanging nodes and other constrained
 * degrees of freedom. When building the global system matrix and the right
 * hand sides, one can build them without taking care of the constraints,
 * purely on a topological base, i.e. by a loop over cells. In order to do
 * actual calculations, you have to 'condense' the linear system: eliminate
 * constrained degrees of freedom and distribute the appropriate values to the
 * unconstrained dofs. This changes the sparsity pattern of the sparse
 * matrices used in finite element calculations and is thus a quite expensive
 * operation. The general scheme of things is then that you build your system,
 * you eliminate (condense) away constrained nodes using the condense()
 * functions of this class, then you solve the remaining system, and finally
 * you compute the values of constrained nodes from the values of the
 * unconstrained ones using the distribute() function. Note that the
 * condense() function is applied to matrix and right hand side of the linear
 * system, while the distribute() function is applied to the solution
 * vector.
 *
 * This scheme of first building a linear system and then eliminating
 * constrained degrees of freedom is inefficient, and a bottleneck if there
 * are many constraints and matrices are full, i.e. especially for 3d and/or
 * higher order or hp finite elements. We therefore offer a second way of
 * building linear systems, using the add_entries_local_to_global() and
 * distribute_local_to_global() functions discussed below. The resulting
 * linear systems are equivalent to what one gets after calling the condense()
 * functions.
 *
 *
 * <h4>Condensing matrices and sparsity patterns</h4>
 *
 * As mentioned above, the first way of using constraints is to build linear
 * systems without regards to constraints and then "condensing" them away.
 * Condensation of a matrix is done in four steps: first one builds the
 * sparsity pattern (e.g. using DoFTools::create_sparsity_pattern()); then the
 * sparsity pattern of the condensed matrix is made out of the original
 * sparsity pattern and the constraints; third, the global matrix is
 * assembled; and fourth, the matrix is finally condensed. To do these steps,
 * you have (at least) two possibilities:
 *
 * <ul>
 * <li> Use two different sparsity patterns and two different matrices: you
 * may eliminate the lines and rows connected with a constraint and create a
 * totally new sparsity pattern and a new system matrix. This has the
 * advantage that the resulting system of equations is smaller and free from
 * artifacts of the condensation process and is therefore faster in the
 * solution process since no unnecessary multiplications occur (see
 * below). However, there are two major drawbacks: keeping two matrices at the
 * same time can be quite unacceptable if you're short of memory. Secondly,
 * the condensation process is expensive, since <em>all</em> entries of the
 * matrix have to be copied, not only those which are subject to constraints.
 *
 * This procedure is therefore not advocated and not discussed in the @ref
 * Tutorial.
 *
 * <li> Use only one sparsity pattern and one matrix: doing it this way, the
 * condense functions add nonzero entries to the sparsity pattern of the large
 * matrix (with constrained nodes in it) where the condensation process of the
 * matrix will create additional nonzero elements. In the condensation process
 * itself, lines and rows subject to constraints are distributed to the lines
 * and rows of unconstrained nodes. The constrained lines remain in place,
 * however, unlike in the first possibility described above. In order not to
 * disturb the solution process, these lines and rows are filled with zeros
 * and an appropriate positive value on the main diagonal (we choose an
 * average of the magnitudes of the other diagonal elements, so as to make
 * sure that the new diagonal entry has the same order of magnitude as the
 * other entries; this preserves the scaling properties of the matrix). The
 * corresponding value in the right hand sides is set to zero. This way, the
 * constrained node will always get the value zero upon solution of the
 * equation system and will not couple to other nodes any more.
 *
 * This method has the advantage that only one matrix and sparsity pattern is
 * needed thus using less memory. Additionally, the condensation process is
 * less expensive, since not all but only constrained values in the matrix
 * have to be copied. On the other hand, the solution process will take a bit
 * longer, since matrix vector multiplications will incur multiplications with
 * zeroes in the lines subject to constraints. Additionally, the vector size
 * is larger than in the first possibility, resulting in more memory
 * consumption for those iterative solution methods using a larger number of
 * auxiliary vectors (e.g. methods using explicit orthogonalization
 * procedures).
 *
 * Nevertheless, this process is overall more efficient due to its lower
 * memory consumption and the one among the two discussed here that is
 * exclusively discussed in the @ref Tutorial.
 * </ul>
 *
 * This class provides two sets of @p condense functions: those taking two
 * arguments refer to the first possibility above, those taking only one do
 * their job in-place and refer to the second possibility.
 *
 * The condensation functions exist for different argument types. The
 * in-place functions (i.e. those following the second way) exist for
 * arguments of type SparsityPattern, SparseMatrix and
 * BlockSparseMatrix. Note that there are no versions for arguments of type
 * PETScWrappers::SparseMatrix() or any of the other PETSc or Trilinos
 * matrix wrapper classes. This is due to the fact that it is relatively
 * hard to get a representation of the sparsity structure of PETSc matrices,
 * and to modify them; this holds in particular, if the matrix is actually
 * distributed across a cluster of computers. If you want to use
 * PETSc/Trilinos matrices, you can either copy an already condensed deal.II
 * matrix, or build the PETSc/Trilinos matrix in the already condensed form,
 * see the discussion below.
 *
 *
 * <h5>Condensing vectors</h5>
 *
 * Condensing vectors works exactly as described above for matrices. Note that
 * condensation is an idempotent operation, i.e. doing it more than once on a
 * vector or matrix yields the same result as doing it only once: once an
 * object has been condensed, further condensation operations don't change it
 * any more.
 *
 * In contrast to the matrix condensation functions, the vector condensation
 * functions exist in variants for PETSc and Trilinos vectors. However,
 * using them is typically expensive, and should be avoided. You should use
 * the same techniques as mentioned above to avoid their use.
 *
 *
 * <h5>Treatment of inhomogeneous constraints</h5>
 *
 * In case some constraint lines have inhomogeneities (which is typically
 * the case if the constraint comes from implementation of inhomogeneous
 * boundary conditions), the situation is a bit more complicated. This is
 * because the elimination of the non-diagonal values in the matrix generate
 * contributions in the eliminated rows in the vector. This means that
 * inhomogeneities can only be handled with functions that act
 * simultaneously on a matrix and a vector. This means that all
 * inhomogeneities are ignored in case the respective condense function is
 * called without any matrix (or if the matrix has already been condensed
 * before).
 *
 * The use of ConstraintMatrix for implementing Dirichlet boundary conditions
 * is discussed in the step-22 tutorial program.
 *
 *
 * <h3>Avoiding explicit condensation</h3>
 *
 * Sometimes, one wants to avoid explicit condensation of a linear system
 * after it has been built at all. There are two main reasons for wanting to
 * do so:
 *
 * <ul>
 * <li>
 * Condensation is an expensive operation, in particular if there
 * are many constraints and/or if the matrix has many nonzero entries. Both
 * is typically the case for 3d, or high polynomial degree computations, as
 * well as for hp finite element methods, see for example the @ref hp_paper
 * "hp paper". This is the case discussed in the hp tutorial program, @ref
 * step_27 "step-27", as well as in step-22 and @ref step_31
 * "step-31".
 *
 * <li>
 * There may not be a condense() function for the matrix you use (this
 * is, for example, the case for the PETSc and Trilinos wrapper classes,
 * where we have no access to the underlying representation of the matrix,
 * and therefore cannot efficiently implement the condense()
 * operation). This is the case discussed in step-17, @ref
 * step_18 "step-18", step-31, and step-32.
 * </ul>
 *
 * In this case, one possibility is to distribute local entries to the final
 * destinations right at the moment of transferring them into the global
 * matrices and vectors, and similarly build a sparsity pattern in the
 * condensed form at the time it is set up originally.
 *
 * This class offers support for these operations as well. For example, the
 * add_entries_local_to_global() function adds nonzero entries to a sparsity
 * pattern object. It not only adds a given entry, but also all entries that
 * we will have to write to if the current entry corresponds to a constrained
 * degree of freedom that will later be eliminated. Similarly, one can use the
 * distribute_local_to_global() functions to directly distribute entries in
 * vectors and matrices when copying local contributions into a global matrix
 * or vector. These calls make a subsequent call to condense() unnecessary.
 *
 * Note that, despite their name which describes what the function really
 * does, the distribute_local_to_global() functions has to be applied to
 * matrices and right hand side vectors, whereas the distribute() function
 * discussed below is applied to the solution vector after solving the linear
 * system.
 *
 *
 * <h3>Distributing constraints</h3>
 *
 * After solving the condensed system of equations, the solution vector has
 * to be redistributed. This is done by the two distribute() functions, one
 * working with two vectors, one working in-place. The operation of
 * distribution undoes the condensation process in some sense, but it should
 * be noted that it is not the inverse operation. Basically, distribution
 * sets the values of the constrained nodes to the value that is computed
 * from the constraint given the values of the unconstrained nodes plus
 * possible inhomogeneities. This is usually necessary since the condensed
 * linear systems only describe the equations for unconstrained nodes, and
 * constrained nodes need to get their values in a second step.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, Martin Kronbichler, 1998, 2004, 2008, 2009
 */
class ConstraintMatrix : public Subscriptor
{
  public:
				     /**
				      * Constructor
				      */
    ConstraintMatrix (const IndexSet & local_constraints = IndexSet());

				     /**
				      * Copy constructor
				      */
    ConstraintMatrix (const ConstraintMatrix &constraint_matrix);

				     /**
				      * Reinit the ConstraintMatrix object and
				      * supply an IndexSet with lines that may
				      * be constrained. This function is only
				      * relevant in the distributed case, to
				      * supply a different IndexSet. Otherwise
				      * this routine is equivalent to calling
				      * clear(). Normally an IndexSet with all
				      * locally_active_dofs should be supplied
				      * here.
				      */
    void reinit (const IndexSet & local_constraints = IndexSet());

				     /**
				      * Determines if we can store a
				      * constraint for the given @p
				      * line_index. This routine only matters
				      * in the distributed case and checks if
				      * the IndexSet allows storage of this
				      * line. Always returns true if not in
				      * the distributed case.
				      */
    bool can_store_line(unsigned int line_index) const;

				     /**
				      * This function copies the content of @p
				      * constraints_in with DoFs that are
				      * element of the IndexSet @p
				      * filter. Elements that are not present
				      * in the IndexSet are ignored. All DoFs
				      * will be transformed to local index
				      * space of the filter, both the
				      * constrained DoFs and the other DoFs
				      * these entries are constrained to. The
				      * local index space of the filter is a
				      * contiguous numbering of all (global)
				      * DoFs that are elements in the
				      * filter.
				      *
				      * If, for example, the filter represents
				      * the range <tt>[10,20)</tt>, and the
				      * constraint matrix @p constraints_in
				      * includes the global indices
				      * <tt>{7,13,14}</tt>, the indices
				      * <tt>{3,4}</tt> are added to the
				      * calling constraint matrix (since 13
				      * and 14 are elements in the filter and
				      * element 13 is the fourth element in
				      * the index, and 14 is the fifth).
				      *
				      * This function provides an easy way to
				      * create a ConstraintMatrix for certain
				      * vector components in a vector-valued
				      * problem from a full ConstraintMatrix,
				      * i.e. extracting a diagonal subblock
				      * from a larger ConstraintMatrix. The
				      * block is specified by the IndexSet
				      * argument.
				      */
    void add_selected_constraints (const ConstraintMatrix &constraints_in,
				   const IndexSet         &filter);

				     /**
				      * @name Adding constraints
				      * @{
				      */

				     /**
				      * Add a new line to the
				      * matrix. If the line already
				      * exists, then the function
				      * simply returns without doing
				      * anything.
				      */
    void add_line (const unsigned int line);

				     /**
				      * Call the first add_line() function for
				      * every index <code>i</code> for which
				      * <code>lines[i]</code> is true.
				      *
				      * This function essentially exists to
				      * allow adding several constraints of
				      * the form $x_i=0$ all at once, where
				      * the set of indices $i$ for which these
				      * constraints should be added are given
				      * by the argument of this function. On
				      * the other hand, just as if the
				      * single-argument add_line() function
				      * were called repeatedly, the
				      * constraints can later be modified to
				      * include linear dependencies using the
				      * add_entry() function as well as
				      * inhomogeneities using
				      * set_inhomogeneity().
				      */
    void add_lines (const std::vector<bool> &lines);

				     /**
				      * Call the first add_line() function for
				      * every index <code>i</code> that
				      * appears in the argument.
				      *
				      * This function essentially exists to
				      * allow adding several constraints of
				      * the form $x_i=0$ all at once, where
				      * the set of indices $i$ for which these
				      * constraints should be added are given
				      * by the argument of this function. On
				      * the other hand, just as if the
				      * single-argument add_line() function
				      * were called repeatedly, the
				      * constraints can later be modified to
				      * include linear dependencies using the
				      * add_entry() function as well as
				      * inhomogeneities using
				      * set_inhomogeneity().
				      */
    void add_lines (const std::set<unsigned int> &lines);

				     /**
				      * Call the first add_line() function for
				      * every index <code>i</code> that
				      * appears in the argument.
				      *
				      * This function essentially exists to
				      * allow adding several constraints of
				      * the form $x_i=0$ all at once, where
				      * the set of indices $i$ for which these
				      * constraints should be added are given
				      * by the argument of this function. On
				      * the other hand, just as if the
				      * single-argument add_line() function
				      * were called repeatedly, the
				      * constraints can later be modified to
				      * include linear dependencies using the
				      * add_entry() function as well as
				      * inhomogeneities using
				      * set_inhomogeneity().
				      */
    void add_lines (const IndexSet &lines);

				     /**
				      * Add an entry to a given
				      * line. The list of lines is
				      * searched from the back to the
				      * front, so clever programming
				      * would add a new line (which is
				      * pushed to the back) and
				      * immediately afterwards fill
				      * the entries of that line. This
				      * way, no expensive searching is
				      * needed.
				      *
				      * If an entry with the same
				      * indices as the one this
				      * function call denotes already
				      * exists, then this function
				      * simply returns provided that
				      * the value of the entry is the
				      * same. Thus, it does no harm to
				      * enter a constraint twice.
				      */
    void add_entry (const unsigned int line,
		    const unsigned int column,
		    const double value);

				     /**
				      * Add a whole series of entries,
				      * denoted by pairs of column indices
				      * and values, to a line of
				      * constraints. This function is
				      * equivalent to calling the preceeding
				      * function several times, but is
				      * faster.
				      */
    void add_entries (const unsigned int                                  line,
		      const std::vector<std::pair<unsigned int,double> > &col_val_pairs);

				     /**
				      * Set an imhomogeneity to the
				      * constraint line <i>i</i>, according
				      * to the discussion in the general
				      * class description.
				      */
    void set_inhomogeneity (const unsigned int line,
			    const double       value);

				     /**
				      * Close the filling of entries. Since
				      * the lines of a matrix of this type
				      * are usually filled in an arbitrary
				      * order and since we do not want to
				      * use associative constainers to store
				      * the lines, we need to sort the lines
				      * and within the lines the columns
				      * before usage of the matrix.  This is
				      * done through this function.
				      *
				      * Also, zero entries are discarded,
				      * since they are not needed.
				      *
				      * After closing, no more entries are
				      * accepted. If the object was already
				      * closed, then this function returns
				      * immediately.
				      *
				      * This function also resolves chains
				      * of constraints. For example, degree
				      * of freedom 13 may be constrained to
				      * $u_{13}=u_3/2+u_7/2$ while degree of
				      * freedom 7 is itself constrained as
				      * $u_7=u_2/2+u_4/2$. Then, the
				      * resolution will be that
				      * $u_{13}=u_3/2+u_2/4+u_4/4$. Note,
				      * however, that cycles in this graph
				      * of constraints are not allowed,
				      * i.e. for example $u_4$ may not be
				      * constrained, directly or indirectly,
				      * to $u_{13}$ again.
				      */
    void close ();

				     /**
				      * Merge the constraints represented by
				      * the object given as argument into
				      * the constraints represented by this
				      * object. Both objects may or may not
				      * be closed (by having their function
				      * @p close called before). If this
				      * object was closed before, then it
				      * will be closed afterwards as
				      * well. Note, however, that if the
				      * other argument is closed, then
				      * merging may be significantly faster.
				      *
				      * Note that the constraints in each of
				      * the two objects (the old one
				      * represented by this object and the
				      * argument) may not refer to the same
				      * degree of freedom, i.e. a degree of
				      * freedom that is constrained in one
				      * object may not be constrained in the
				      * second. If this is nevertheless the
				      * case, an exception is thrown.
				      *
				      * However, the following is possible:
				      * if DoF @p x is constrained to dofs
				      * @p x_i for some set of indices @p i,
				      * then the DoFs @p x_i may be further
				      * constrained by the constraints
				      * object given as argument, although
				      * not to other DoFs that are
				      * constrained in either of the two
				      * objects. Note that it is not
				      * possible that the DoFs @p x_i are
				      * constrained within the present
				      * object.
				      *
				      * Because of simplicity of
				      * implementation, and also to avoid
				      * cycles, this operation is not
				      * symmetric: degrees of freedom that
				      * are constrained in the given
				      * argument object may not be
				      * constrained to DoFs that are
				      * themselves constrained within the
				      * present object.
				      *
				      * The aim of these merging operations
				      * is that if, for example, you have
				      * hanging nodes that are constrained
				      * to the degrees of freedom adjacent
				      * to them, you cannot originally,
				      * i.e. within one object, constrain
				      * these adjacent nodes
				      * further. However, that may be
				      * desirable in some cases, for example
				      * if they belong to a symmetry
				      * boundary for which the nodes on one
				      * side of the domain should have the
				      * same values as those on the other
				      * side. In that case, you would first
				      * construct a costraints object
				      * holding the hanging nodes
				      * constraints, and a second one that
				      * contains the constraints due to the
				      * symmetry boundary. You would then
				      * finally merge this second one into
				      * the first, possibly eliminating
				      * constraints of hanging nodes to
				      * adjacent boundary nodes by
				      * constraints to nodes at the opposite
				      * boundary.
				      */
    void merge (const ConstraintMatrix &other_constraints);

				     /**
				      * Shift all entries of this matrix
				      * down @p offset rows and over @p
				      * offset columns.
				      *
				      * This function is useful if you are
				      * building block matrices, where all
				      * blocks are built by the same @p
				      * DoFHandler object, i.e. the matrix
				      * size is larger than the number of
				      * degrees of freedom. Since several
				      * matrix rows and columns correspond
				      * to the same degrees of freedom,
				      * you'd generate several constraint
				      * objects, then shift them, and
				      * finally @p merge them together
				      * again.
				      */
    void shift (const unsigned int offset);

				     /**
				      * Clear all entries of this
				      * matrix. Reset the flag determining
				      * whether new entries are accepted or
				      * not.
				      *
				      * This function may be called also on
				      * objects which are empty or already
				      * cleared.
				      */
    void clear ();

				     /**
				      * @}
				      */


				     /**
				      * @name Querying constraints
				      * @{
				      */

				     /**
				      * Return number of constraints stored in
				      * this matrix.
				      */
    unsigned int n_constraints () const;

				     /**
				      * Return whether the degree of freedom
				      * with number @p index is a
				      * constrained one.
				      *
				      * Note that if @p close was called
				      * before, then this function is
				      * significantly faster, since then the
				      * constrained degrees of freedom are
				      * sorted and we can do a binary
				      * search, while before @p close was
				      * called, we have to perform a linear
				      * search through all entries.
				      */
    bool is_constrained (const unsigned int index) const;

				     /**
				      * Return whether the dof is
				      * constrained, and whether it is
				      * constrained to only one other degree
				      * of freedom with weight one. The
				      * function therefore returns whether
				      * the degree of freedom would simply
				      * be eliminated in favor of exactly
				      * one other degree of freedom.
				      *
				      * The function returns @p false if
				      * either the degree of freedom is not
				      * constrained at all, or if it is
				      * constrained to more than one other
				      * degree of freedom, or if it is
				      * constrained to only one degree of
				      * freedom but with a weight different
				      * from one.
				      */
    bool is_identity_constrained (const unsigned int index) const;

				     /**
				      * Return the maximum number of other
				      * dofs that one dof is constrained
				      * to. For example, in 2d a hanging
				      * node is constrained only to its two
				      * neighbors, so the returned value
				      * would be @p 2. However, for higher
				      * order elements and/or higher
				      * dimensions, or other types of
				      * constraints, this number is no more
				      * obvious.
				      *
				      * The name indicates that within the
				      * system matrix, references to a
				      * constrained node are indirected to
				      * the nodes it is constrained to.
				      */
    unsigned int max_constraint_indirections () const;

				     /**
				      * Returns <tt>true</tt> in case the
				      * dof is constrained and there is a
				      * non-trivial inhomogeneous valeus set
				      * to the dof.
				      */
    bool is_inhomogeneously_constrained (const unsigned int index) const;

				     /**
				      * Returns <tt>false</tt> if all
				      * constraints in the ConstraintMatrix
				      * are homogeneous ones, and
				      * <tt>true</tt> if there is at least
				      * one inhomogeneity.
				      */
    bool has_inhomogeneities () const;

                                     /**
                                      * Returns a pointer to the the vector of
                                      * entries if a line is constrained, and a
                                      * zero pointer in case the dof is not
                                      * constrained.
                                      */
    const std::vector<std::pair<unsigned int,double> >*
    get_constraint_entries (unsigned int line) const;

				     /**
				      * Print the constraint lines. Mainly
				      * for debugging purposes.
				      *
				      * This function writes out all entries
				      * in the constraint matrix lines with
				      * their value in the form <tt>row col
				      * : value</tt>. Unconstrained lines
				      * containing only one identity entry
				      * are not stored in this object and
				      * are not printed.
				      */
    void print (std::ostream &) const;

				     /**
				      * Write the graph of constraints in
				      * 'dot' format. 'dot' is a program
				      * that can take a list of nodes and
				      * produce a graphical representation
				      * of the graph of constrained degrees
				      * of freedom and the degrees of
				      * freedom they are constrained to.
				      *
				      * The output of this function can be
				      * used as input to the 'dot' program
				      * that can convert the graph into a
				      * graphical representation in
				      * postscript, png, xfig, and a number
				      * of other formats.
				      *
				      * This function exists mostly for
				      * debugging purposes.
				      */
    void write_dot (std::ostream &) const;

				     /**
				      * Determine an estimate for the memory
				      * consumption (in bytes) of this
				      * object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * @}
				      */

				     /**
				      * @name Eliminating constraints from linear systems after their creation
				      * @{
				      */

				     /**
				      * Condense a given sparsity
				      * pattern. This function assumes the
				      * uncondensed matrix struct to be
				      * compressed and the one to be filled
				      * to be empty. The condensed structure
				      * is compressed afterwards.
				      *
				      * The constraint matrix object must be
				      * closed to call this function.
				      *
				      * @note The hanging nodes are
				      * completely eliminated from the
				      * linear system refering to
				      * <tt>condensed</tt>. Therefore, the
				      * dimension of <tt>condensed</tt> is
				      * the dimension of
				      * <tt>uncondensed</tt> minus the
				      * number of constrained degrees of
				      * freedom.
				      */
    void condense (const SparsityPattern &uncondensed,
		   SparsityPattern       &condensed) const;


				     /**
				      * This function does much the same as
				      * the above one, except that it
				      * condenses the matrix struct
				      * 'in-place'. It does not remove
				      * nonzero entries from the matrix but
				      * adds those needed for the process of
				      * distribution of the constrained
				      * degrees of freedom.
				      *
				      * Since this function adds new nonzero
				      * entries to the sparsity pattern, the
				      * argument must not be
				      * compressed. However the constraint
				      * matrix must be closed. The matrix
				      * struct is compressed at the end of
				      * the function.
				      */
    void condense (SparsityPattern &sparsity) const;

				     /**
				      * Same function as above, but
				      * condenses square block sparsity
				      * patterns.
				      */
    void condense (BlockSparsityPattern &sparsity) const;

				     /**
				      * Same function as above, but
				      * condenses square compressed sparsity
				      * patterns.
				      *
				      * Given the data structure used by
				      * CompressedSparsityPattern, this
				      * function becomes quadratic in the
				      * number of degrees of freedom for
				      * large problems and can dominate
				      * setting up linear systems when
				      * several hundred thousand or millions
				      * of unknowns are involved and for
				      * problems with many nonzero elements
				      * per row (for example for
				      * vector-valued problems or hp finite
				      * elements). In this case, it is
				      * advisable to use the
				      * CompressedSetSparsityPattern class
				      * instead, see for example @ref
				      * step_27 "step-27", or to use the
				      * CompressedSimpleSparsityPattern
				      * class, see for example @ref step_31
				      * "step-31".
				      */
    void condense (CompressedSparsityPattern &sparsity) const;

				     /**
				      * Same function as above, but
				      * condenses compressed sparsity
				      * patterns, which are based on the
				      * std::set container.
				      */
    void condense (CompressedSetSparsityPattern &sparsity) const;

				     /**
				      * Same function as above, but
				      * condenses compressed sparsity
				      * patterns, which are based on the
				      * ''simple'' aproach.
				      */
    void condense (CompressedSimpleSparsityPattern &sparsity) const;

				     /**
				      * Same function as above, but
				      * condenses square compressed sparsity
				      * patterns.
				      *
				      * Given the data structure used by
				      * BlockCompressedSparsityPattern, this
				      * function becomes quadratic in the
				      * number of degrees of freedom for
				      * large problems and can dominate
				      * setting up linear systems when
				      * several hundred thousand or millions
				      * of unknowns are involved and for
				      * problems with many nonzero elements
				      * per row (for example for
				      * vector-valued problems or hp finite
				      * elements). In this case, it is
				      * advisable to use the
				      * BlockCompressedSetSparsityPattern
				      * class instead, see for example @ref
				      * step_27 "step-27" and @ref step_31
				      * "step-31".
				      */
    void condense (BlockCompressedSparsityPattern &sparsity) const;

				     /**
				      * Same function as above, but
				      * condenses square compressed sparsity
				      * patterns.
				      */
    void condense (BlockCompressedSetSparsityPattern &sparsity) const;

				     /**
				      * Same function as above, but
				      * condenses square compressed sparsity
				      * patterns.
				      */
    void condense (BlockCompressedSimpleSparsityPattern &sparsity) const;


				     /**
				      * Condense a given matrix. The
				      * associated matrix struct should be
				      * condensed and compressed. It is the
				      * user's responsibility to guarantee
				      * that all entries in the @p condensed
				      * matrix be zero!
				      *
				      * The constraint matrix object must be
				      * closed to call this function.
				      */
    template<typename number>
    void condense (const SparseMatrix<number> &uncondensed,
		   SparseMatrix<number>       &condensed) const;

				     /**
				      * This function does much the same as
				      * the above one, except that it
				      * condenses the matrix 'in-place'. See
				      * the general documentation of this
				      * class for more detailed information.
				      */
    template<typename number>
    void condense (SparseMatrix<number> &matrix) const;

				     /**
				      * Same function as above, but
				      * condenses square block sparse
				      * matrices.
				      */
    template <typename number>
    void condense (BlockSparseMatrix<number> &matrix) const;

				     /**
				      * Condense the given vector @p
				      * uncondensed into @p condensed. It is
				      * the user's responsibility to
				      * guarantee that all entries of @p
				      * condensed be zero. Note that this
				      * function does not take any
				      * inhomogeneity into account and
				      * throws an exception in case there
				      * are any inhomogeneities. Use
				      * the function using both a matrix and
				      * vector for that case.
				      *
				      * The @p VectorType may be a
				      * Vector<float>, Vector<double>,
				      * BlockVector<tt><...></tt>, a PETSc
				      * or Trilinos vector wrapper class, or
				      * any other type having the same
				      * interface.
				      */
    template <class VectorType>
    void condense (const VectorType &uncondensed,
		   VectorType       &condensed) const;

				     /**
				      * Condense the given vector
				      * in-place. The @p VectorType may be a
				      * Vector<float>, Vector<double>,
				      * BlockVector<tt><...></tt>, a PETSc
				      * or Trilinos vector wrapper class, or
				      * any other type having the same
				      * interface. Note that this function
				      * does not take any inhomogeneity into
				      * account and throws an exception in
				      * case there are any
				      * inhomogeneities. Use the function
				      * using both a matrix and vector for
				      * that case.
				      */
    template <class VectorType>
    void condense (VectorType &vec) const;

				     /**
				      * Condense a given matrix and a given
				      * vector. The associated matrix struct
				      * should be condensed and
				      * compressed. It is the user's
				      * responsibility to guarantee that all
				      * entries in the @p condensed matrix
				      * and vector be zero! This function is
				      * the appropriate choice for applying
				      * inhomogeneous constraints.
				      *
				      * The constraint matrix object must be
				      * closed to call this function.
				      */
    template<typename number, class VectorType>
    void condense (const SparseMatrix<number> &uncondensed_matrix,
		   const VectorType           &uncondensed_vector,
		   SparseMatrix<number>       &condensed_matrix,
		   VectorType                 &condensed_vector) const;

				     /**
				      * This function does much the same as
				      * the above one, except that it
				      * condenses matrix and vector
				      * 'in-place'. See the general
				      * documentation of this class for more
				      * detailed information.
				      */
    template<typename number, class VectorType>
    void condense (SparseMatrix<number> &matrix,
		   VectorType           &vector) const;

				     /**
				      * Same function as above, but
				      * condenses square block sparse
				      * matrices and vectors.
				      */
    template <typename number, class BlockVectorType>
    void condense (BlockSparseMatrix<number> &matrix,
		   BlockVectorType           &vector) const;

				     /**
				      * Delete hanging nodes in a
				      * vector. Sets all hanging node
				      * values to zero. The @p
				      * VectorType may be a
				      * Vector<float>, Vector<double>,
				      * BlockVector<tt><...></tt>, a
				      * PETSc or Trilinos vector
				      * wrapper class, or any other
				      * type having the same
				      * interface.
				      */
    template <class VectorType>
    void set_zero (VectorType &vec) const;

				     /**
				      * @}
				      */

				     /**
				      * @name Eliminating constraints from linear systems during their creation
				      * @{
				      */

                                     /**
                                      * This function takes a vector of
                                      * local contributions (@p
                                      * local_vector) corresponding to the
                                      * degrees of freedom indices given in
                                      * @p local_dof_indices and distributes
                                      * them to the global vector. In most
                                      * cases, these local contributions
                                      * will be the result of an integration
                                      * over a cell or face of a
                                      * cell. However, as long as @p
                                      * local_vector and @p
                                      * local_dof_indices have the same
                                      * number of elements, this function is
                                      * happy with whatever it is
                                      * given.
                                      *
                                      * In contrast to the similar function
                                      * in the DoFAccessor class, this
                                      * function also takes care of
                                      * constraints, i.e. if one of the
                                      * elements of @p local_dof_indices
                                      * belongs to a constrained node, then
                                      * rather than writing the
                                      * corresponding element of @p
                                      * local_vector into @p global_vector,
                                      * the element is distributed to the
                                      * entries in the global vector to
                                      * which this particular degree of
                                      * freedom is constrained.
                                      *
                                      * Thus, by using this function to
                                      * distribute local contributions to the
                                      * global object, one saves the call to
                                      * the condense function after the
                                      * vectors and matrices are fully
                                      * assembled. On the other hand, by
                                      * consequence, the function does not
                                      * only write into the entries enumerated
                                      * by the @p local_dof_indices array, but
                                      * also (possibly) others as necessary.
				      *
				      * Note that this function will apply all
                                      * constraints as if they were
                                      * homogeneous. For correctly setting
                                      * inhomogeneous constraints, use the
                                      * similar function with a matrix
                                      * argument or the function with both
                                      * matrix and vector arguments.
				      *
				      * Note: This function is not
				      * thread-safe, so you will need to make
				      * sure that only on process at a time
				      * calls this function.
                                      */
    template <class InVector, class OutVector>
    void
    distribute_local_to_global (const InVector                  &local_vector,
                                const std::vector<unsigned int> &local_dof_indices,
                                OutVector                       &global_vector) const;

                                     /**
                                      * This function takes a vector of
                                      * local contributions (@p
                                      * local_vector) corresponding to the
                                      * degrees of freedom indices given in
                                      * @p local_dof_indices and distributes
                                      * them to the global vector. In most
                                      * cases, these local contributions
                                      * will be the result of an integration
                                      * over a cell or face of a
                                      * cell. However, as long as @p
                                      * local_vector and @p
                                      * local_dof_indices have the same
                                      * number of elements, this function is
                                      * happy with whatever it is
                                      * given.
                                      *
                                      * In contrast to the similar function in
                                      * the DoFAccessor class, this function
                                      * also takes care of constraints,
                                      * i.e. if one of the elements of @p
                                      * local_dof_indices belongs to a
                                      * constrained node, then rather than
                                      * writing the corresponding element of
                                      * @p local_vector into @p global_vector,
                                      * the element is distributed to the
                                      * entries in the global vector to which
                                      * this particular degree of freedom is
                                      * constrained.
                                      *
                                      * Thus, by using this function to
                                      * distribute local contributions to the
                                      * global object, one saves the call to
                                      * the condense function after the
                                      * vectors and matrices are fully
                                      * assembled. On the other hand, by
                                      * consequence, the function does not
                                      * only write into the entries enumerated
                                      * by the @p local_dof_indices array, but
                                      * also (possibly) others as
                                      * necessary. This includes writing into
                                      * diagonal elements of the matrix if the
                                      * corresponding degree of freedom is
                                      * constrained.
				      *
				      * The fourth argument
				      * <tt>local_matrix</tt> is intended to
				      * be used in case one wants to apply
				      * inhomogeneous constraints on the
				      * vector only. Such a situation could be
				      * where one wants to assemble of a right
				      * hand side vector on a problem with
				      * inhomogeneous constraints, but the
				      * global matrix has been assembled
				      * previously. A typical example of this
				      * is a time stepping algorithm where the
				      * stiffness matrix is assembled once,
				      * and the right hand side updated every
				      * time step. Note that, however, the
				      * entries in the columns of the local
				      * matrix have to be exactly the same as
				      * those that have been written into the
				      * global matrix. Otherwise, this
				      * function will not be able to correctly
				      * handle inhomogeneities.
				      *
				      * Note: This function is not
				      * thread-safe, so you will need to make
				      * sure that only on process at a time
				      * calls this function.
                                      */
    template <typename VectorType>
    void
    distribute_local_to_global (const Vector<double>            &local_vector,
                                const std::vector<unsigned int> &local_dof_indices,
                                VectorType                      &global_vector,
				const FullMatrix<double>        &local_matrix) const;


                                     /**
                                      * This function takes a pointer to a
                                      * vector of local contributions (@p
                                      * local_vector) corresponding to the
                                      * degrees of freedom indices given in
                                      * @p local_dof_indices and distributes
                                      * them to the global vector. In most
                                      * cases, these local contributions
                                      * will be the result of an integration
                                      * over a cell or face of a
                                      * cell. However, as long as the
                                      * entries in @p local_dof_indices
                                      * indicate reasonable global vector
                                      * entries, this function is happy with
                                      * whatever it is given.
                                      *
                                      * If one of the elements of @p
                                      * local_dof_indices belongs to a
                                      * constrained node, then rather than
                                      * writing the corresponding element of
                                      * @p local_vector into @p
                                      * global_vector, the element is
                                      * distributed to the entries in the
                                      * global vector to which this
                                      * particular degree of freedom is
                                      * constrained.
                                      *
                                      * Thus, by using this function to
                                      * distribute local contributions to
                                      * the global object, one saves the
                                      * call to the condense function after
                                      * the vectors and matrices are fully
                                      * assembled. Note that this function
                                      * completely ignores inhomogeneous
                                      * constraints.
				      *
				      * Note: This function is not
				      * thread-safe, so you will need to
				      * make sure that only on process at a
				      * time calls this function.
                                      */
    template <typename ForwardIteratorVec, typename ForwardIteratorInd,
              class VectorType>
    void
    distribute_local_to_global (ForwardIteratorVec local_vector_begin,
                                ForwardIteratorVec local_vector_end,
				ForwardIteratorInd local_indices_begin,
                                VectorType        &global_vector) const;

                                     /**
                                      * This function takes a matrix of
                                      * local contributions (@p
                                      * local_matrix) corresponding to the
                                      * degrees of freedom indices given in
                                      * @p local_dof_indices and distributes
                                      * them to the global matrix. In most
                                      * cases, these local contributions
                                      * will be the result of an integration
                                      * over a cell or face of a
                                      * cell. However, as long as @p
                                      * local_matrix and @p
                                      * local_dof_indices have the same
                                      * number of elements, this function is
                                      * happy with whatever it is given.
                                      *
                                      * In contrast to the similar function
                                      * in the DoFAccessor class, this
                                      * function also takes care of
                                      * constraints, i.e. if one of the
                                      * elements of @p local_dof_indices
                                      * belongs to a constrained node, then
                                      * rather than writing the
                                      * corresponding element of @p
                                      * local_matrix into @p global_matrix,
                                      * the element is distributed to the
                                      * entries in the global matrix to
                                      * which this particular degree of
                                      * freedom is constrained.
                                      *
                                      * With this scheme, we never write
                                      * into rows or columns of constrained
                                      * degrees of freedom. In order to make
                                      * sure that the resulting matrix can
                                      * still be inverted, we need to do
                                      * something with the diagonal elements
                                      * corresponding to constrained
                                      * nodes. Thus, if a degree of freedom
                                      * in @p local_dof_indices is
                                      * constrained, we distribute the
                                      * corresponding entries in the matrix,
                                      * but also add the absolute value of
                                      * the diagonal entry of the local
                                      * matrix to the corresponding entry in
                                      * the global matrix. Since the exact
                                      * value of the diagonal element is not
                                      * important (the value of the
                                      * respective degree of freedom will be
                                      * overwritten by the distribute() call
                                      * later on anyway), this guarantees
                                      * that the diagonal entry is always
                                      * non-zero, positive, and of the same
                                      * order of magnitude as the other
                                      * entries of the matrix.
                                      *
                                      * Thus, by using this function to
                                      * distribute local contributions to
                                      * the global object, one saves the
                                      * call to the condense function after
                                      * the vectors and matrices are fully
                                      * assembled.
				      *
				      * Note: This function is not
				      * thread-safe, so you will need to
				      * make sure that only on process at a
				      * time calls this function.
                                      */
    template <typename MatrixType>
    void
    distribute_local_to_global (const FullMatrix<double>        &local_matrix,
                                const std::vector<unsigned int> &local_dof_indices,
                                MatrixType                      &global_matrix) const;

                                     /**
                                      * Does the same as the function
                                      * above but can treat non
                                      * quadratic matrices.
                                      */
    template <typename MatrixType>
    void
    distribute_local_to_global (const FullMatrix<double>        &local_matrix,
                                const std::vector<unsigned int> &row_indices,
                                const std::vector<unsigned int> &col_indices,
                                MatrixType                      &global_matrix) const;

				     /**
				      * This function simultaneously writes
				      * elements into matrix and vector,
				      * according to the constraints
				      * specified by the calling
				      * ConstraintMatrix. This function can
				      * correctly handle inhomogeneous
				      * constraints as well.
				      *
				      * Note: This function is not
				      * thread-safe, so you will need to
				      * make sure that only on process at a
				      * time calls this function.
				      */
    template <typename MatrixType, typename VectorType>
    void
    distribute_local_to_global (const FullMatrix<double>        &local_matrix,
				const Vector<double>            &local_vector,
                                const std::vector<unsigned int> &local_dof_indices,
                                MatrixType                      &global_matrix,
				VectorType                      &global_vector) const;

				     /**
				      * Do a similar operation as the
				      * distribute_local_to_global() function
				      * that distributed writing entries into
				      * a matrix for constrained degrees of
				      * freedom, except that here we don't
				      * write into a matrix but only allocate
				      * sparsity pattern entries.
				      *
				      * As explained in the
				      * @ref hp_paper "hp paper"
				      * and in step-27,
				      * first allocating a sparsity pattern
				      * and later coming back and allocating
				      * additional entries for those matrix
				      * entries that will be written to due to
				      * the elimination of constrained degrees
				      * of freedom (using
				      * ConstraintMatrix::condense() ), can be
				      * a very expensive procedure. It is
				      * cheaper to allocate these entries
				      * right away without having to do a
				      * second pass over the sparsity pattern
				      * object. This function does exactly
				      * that.
				      *
				      * Because the function only allocates
				      * entries in a sparsity pattern, all it
				      * needs to know are the degrees of
				      * freedom that couple to each
				      * other. Unlike the previous function,
				      * no actual values are written, so the
				      * second input argument is not necessary
				      * here.
				      *
				      * The third argument to this function,
				      * keep_constrained_entries determines
				      * whether the function shall allocate
				      * entries in the sparsity pattern at
				      * all for entries that will later be
				      * set to zero upon condensation of the
				      * matrix. These entries are necessary
				      * if the matrix is built
				      * unconstrained, and only later
				      * condensed. They are not necessary if
				      * the matrix is built using the
				      * distribute_local_to_global()
				      * function of this class which
				      * distributes entries right away when
				      * copying a local matrix into a global
				      * object. The default of this argument
				      * is true, meaning to allocate the few
				      * entries that may later be set to
				      * zero.
				      *
				      * By default, the function adds
				      * entries for all pairs of indices
				      * given in the first argument to the
				      * sparsity pattern (unless
				      * keep_constrained_entries is
				      * false). However, sometimes one would
				      * like to only add a subset of all of
				      * these pairs. In that case, the last
				      * argument can be used which specifies
				      * a boolean mask which of the pairs of
				      * indices should be considered. If the
				      * mask is false for a pair of indices,
				      * then no entry will be added to the
				      * sparsity pattern for this pair,
				      * irrespective of whether one or both
				      * of the indices correspond to
				      * constrained degrees of freedom.
				      *
				      * This function is not typically called
				      * from user code, but is used in the
				      * DoFTools::make_sparsity_pattern()
				      * function when passed a constraint
				      * matrix object.
				      */
    template <typename SparsityType>
    void
    add_entries_local_to_global (const std::vector<unsigned int> &local_dof_indices,
				 SparsityType                    &sparsity_pattern,
				 const bool                       keep_constrained_entries = true,
				 const Table<2,bool>             &dof_mask = default_empty_table) const;

				     /**
				      * Similar to the other function,
				      * but for non-quadratic sparsity
				      * patterns.
				      */

    template <typename SparsityType>
    void
    add_entries_local_to_global (const std::vector<unsigned int> &row_indices,
				 const std::vector<unsigned int> &col_indices,
				 SparsityType                    &sparsity_pattern,
				 const bool                       keep_constrained_entries = true,
				 const Table<2,bool>             &dof_mask = default_empty_table) const;

                                     /**
                                      * This function imports values from a
                                      * global vector (@p global_vector) by
                                      * applying the constraints to a vector
                                      * of local values, expressed in
                                      * iterator format.  In most cases, the
                                      * local values will be identified by
                                      * the local dof values on a
                                      * cell. However, as long as the
                                      * entries in @p local_dof_indices
                                      * indicate reasonable global vector
                                      * entries, this function is happy with
                                      * whatever it is given.
                                      *
                                      * If one of the elements of @p
                                      * local_dof_indices belongs to a
                                      * constrained node, then rather than
                                      * writing the corresponding element of
                                      * @p global_vector into @p
                                      * local_vector, the constraints are
                                      * resolved as the respective
                                      * distribute function does, i.e., the
                                      * local entry is constructed from the
                                      * global entries to which this
                                      * particular degree of freedom is
                                      * constrained.
				      *
                                      * In contrast to the similar function
                                      * get_dof_values in the DoFAccessor
                                      * class, this function does not need
                                      * the constrained values to be
                                      * correctly set (i.e., distribute to
                                      * be called).
                                      */
    template <typename ForwardIteratorVec, typename ForwardIteratorInd,
              class VectorType>
    void
    get_dof_values (const VectorType  &global_vector,
		    ForwardIteratorInd local_indices_begin,
		    ForwardIteratorVec local_vector_begin,
		    ForwardIteratorVec local_vector_end) const;

				     /**
				      * @}
				      */

				     /**
				      * @name Dealing with constraints after solving a linear system
				      * @{
				      */

				     /**
				      * Re-distribute the elements of the
				      * vector @p condensed to @p
				      * uncondensed. It is the user's
				      * responsibility to guarantee that all
				      * entries of @p uncondensed be zero!
				      *
				      * This function undoes the action of
				      * @p condense somehow, but it should
				      * be noted that it is not the inverse
				      * of @p condense.
				      *
				      * The @p VectorType may be a
				      * Vector<float>, Vector<double>,
				      * BlockVector<tt><...></tt>, a PETSc
				      * or Trilinos vector wrapper class, or
				      * any other type having the same
				      * interface.
				      */
    template <class VectorType>
    void distribute (const VectorType &condensed,
		     VectorType       &uncondensed) const;

				     /**
				      * Re-distribute the elements of the
				      * vector in-place. The @p VectorType
				      * may be a Vector<float>,
				      * Vector<double>,
				      * BlockVector<tt><...></tt>, a PETSc
				      * or Trilinos vector wrapper class, or
				      * any other type having the same
				      * interface.
				      */
    template <class VectorType>
    void distribute (VectorType &vec) const;

				     /**
				      * @}
				      */

				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcMatrixIsClosed);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException0 (ExcMatrixNotClosed);
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException1 (ExcLineInexistant,
		    unsigned int,
		    << "The specified line " << arg1
		    << " does not exist.");
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException4 (ExcEntryAlreadyExists,
		    int, int, double, double,
		    << "The entry for the indices " << arg1 << " and "
		    << arg2 << " already exists, but the values "
		    << arg3 << " (old) and " << arg4 << " (new) differ "
                    << "by " << (arg4-arg3) << ".");
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException2 (ExcDoFConstrainedToConstrainedDoF,
		    int, int,
		    << "You tried to constrain DoF " << arg1
		    << " to DoF " << arg2
		    << ", but that one is also constrained. This is not allowed!");
				     /**
				      * Exception.
				      *
				      * @ingroup Exceptions
				      */
    DeclException1 (ExcDoFIsConstrainedFromBothObjects,
		    int,
		    << "Degree of freedom " << arg1
		    << " is constrained from both object in a merge operation.");
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException1 (ExcDoFIsConstrainedToConstrainedDoF,
		    int,
		    << "In the given argument a degree of freedom is constrained "
		    << "to another DoF with number " << arg1
		    << ", which however is constrained by this object. This is not"
		    << " allowed.");
				     /**
				      * Exception
				      *
				      * @ingroup Exceptions
				      */
    DeclException1 (ExcRowNotStoredHere,
		    int,
		    << "The index set given to this constraint matrix indicates "
		    << "constraints for degree of freedom " << arg1
		    << " should not be stored by this object, but a constraint "
		    << "is being added.");

  private:

				     /**
				      * This class represents one line of a
				      * constraint matrix.
				      */
    struct ConstraintLine
    {
				         /**
					  * A data type in which we store the list
					  * of entries that make up the homogenous
					  * part of a constraint.
					  */
	typedef std::vector<std::pair<unsigned int,double> > Entries;

					 /**
					  * Number of this line. Since only
					  * very few lines are stored, we
					  * can not assume a specific order
					  * and have to store the line
					  * number explicitly.
					  */
	unsigned int line;

					 /**
					  * Row numbers and values of the
					  * entries in this line.
					  *
					  * For the reason why we use a
					  * vector instead of a map and the
					  * consequences thereof, the same
					  * applies as what is said for
					  * ConstraintMatrix::lines.
					  */
	Entries entries;

				         /**
					  * Value of the inhomogeneity.
					  */
	double inhomogeneity;

                                         /**
					  * This operator is a bit weird and
					  * unintuitive: it compares the
					  * line numbers of two lines. We
					  * need this to sort the lines; in
					  * fact we could do this using a
					  * comparison predicate.  However,
					  * this way, it is easier, albeit
					  * unintuitive since two lines
					  * really have no god-given order
					  * relation.
					  */
	bool operator < (const ConstraintLine &) const;

                                         /**
					  * This operator is likewise weird:
					  * it checks whether the line
					  * indices of the two operands are
					  * equal, irrespective of the fact
					  * that the contents of the line
					  * may be different.
					  */
	bool operator == (const ConstraintLine &) const;

					 /**
					  * Determine an estimate for the
					  * memory consumption (in bytes) of
					  * this object.
					  */
	unsigned int memory_consumption () const;
    };

				     /**
				      * Store the lines of the matrix.
				      * Entries are usually appended in an
				      * arbitrary order and insertion into a
				      * vector is done best at the end, so
				      * the order is unspecified after all
				      * entries are inserted. Sorting of the
				      * entries takes place when calling the
				      * <tt>close()</tt> function.
				      *
				      * We could, instead of using a vector,
				      * use an associative array, like a map
				      * to store the lines. This, however,
				      * would mean a much more fractioned
				      * heap since it allocates many small
				      * objects, and would additionally make
				      * usage of this matrix much slower.
				      */
    std::vector<ConstraintLine> lines;

				     /**
				      * A list of unsigned integers that
				      * contains the position of the
				      * ConstraintLine of a constrained degree
				      * of freedom, or @p
				      * numbers::invalid_unsigned_int if the
				      * degree of freedom is not
				      * constrained. The @p invalid_unsigned
				      * int return value returns thus whether
				      * there is a constraint line for a given
				      * degree of freedom index. Note that
				      * this class has no notion of how many
				      * degrees of freedom there really are,
				      * so if we check whether there is a
				      * constraint line for a given degree of
				      * freedom, then this vector may actually
				      * be shorter than the index of the DoF
				      * we check for.
				      *
				      * This field exists since when adding a
				      * new constraint line we have to figure
				      * out whether it already
				      * exists. Previously, we would simply
				      * walk the unsorted list of constraint
				      * lines until we either hit the end or
				      * found it. This algorithm is O(N) if N
				      * is the number of constraints, which
				      * makes it O(N^2) when inserting all
				      * constraints. For large problems with
				      * many constraints, this could easily
				      * take 5-10 per cent of the total run
				      * time. With this field, we can save
				      * this time since we find any constraint
				      * in O(1) time or get to know that it a
				      * certain degree of freedom is not
				      * constrained.
				      *
				      * To make things worse, traversing the
				      * list of existing constraints requires
				      * reads from many different places in
				      * memory. Thus, in large 3d
				      * applications, the add_line() function
				      * showed up very prominently in the
				      * overall compute time, mainly because
				      * it generated a lot of cache
				      * misses. This should also be fixed by
				      * using the O(1) algorithm to access the
				      * fields of this array.
				      *
				      * The field is useful in a number of
				      * other contexts as well, e.g. when one
				      * needs random access to the constraints
				      * as in all the functions that apply
				      * constraints on the fly while add cell
				      * contributions into vectors and
				      * matrices.
				      */
    std::vector<unsigned int> lines_cache;

				     /**
				      * This IndexSet is used to limit the
				      * lines to save in the ContraintMatrix
				      * to a subset. This is necessary,
				      * because the lines_cache vector would
				      * become too big in a distributed
				      * calculation.
				      */
    IndexSet local_lines;

				     /**
				      * Store whether the arrays are sorted.
				      * If so, no new entries can be added.
				      */
    bool sorted;

				     /**
				      * Internal function to calculate the
				      * index of line @p line in the vector
				      * lines_cache using local_lines.
				      */
    unsigned int calculate_line_index (const unsigned int line) const;

				     /**
				      * Return @p true if the weight of an
				      * entry (the second element of the
				      * pair) equals zero. This function is
				      * used to delete entries with zero
				      * weight.
				      */
    static bool check_zero_weight (const std::pair<unsigned int, double> &p);

 				     /**
				      * Dummy table that serves as default
				      * argument for function
				      * <tt>add_entries_local_to_global()</tt>.
				      */
    static const Table<2,bool> default_empty_table;

				     /**
				      * This function actually implements
				      * the local_to_global function for
				      * standard (non-block) matrices.
				      */
    template <typename MatrixType, typename VectorType>
    void
    distribute_local_to_global (const FullMatrix<double>        &local_matrix,
				const Vector<double>            &local_vector,
                                const std::vector<unsigned int> &local_dof_indices,
                                MatrixType                      &global_matrix,
				VectorType                      &global_vector,
				internal::bool2type<false>) const;

				     /**
				      * This function actually implements
				      * the local_to_global function for
				      * block matrices.
				      */
    template <typename MatrixType, typename VectorType>
    void
    distribute_local_to_global (const FullMatrix<double>        &local_matrix,
				const Vector<double>            &local_vector,
                                const std::vector<unsigned int> &local_dof_indices,
                                MatrixType                      &global_matrix,
				VectorType                      &global_vector,
				internal::bool2type<true>) const;

				     /**
				      * This function actually implements
				      * the local_to_global function for
				      * standard (non-block) sparsity types.
				      */
    template <typename SparsityType>
    void
    add_entries_local_to_global (const std::vector<unsigned int> &local_dof_indices,
				 SparsityType                    &sparsity_pattern,
				 const bool                       keep_constrained_entries,
				 const Table<2,bool>             &dof_mask,
				 internal::bool2type<false>) const;

				     /**
				      * This function actually implements
				      * the local_to_global function for
				      * block sparsity types.
				      */
    template <typename SparsityType>
    void
    add_entries_local_to_global (const std::vector<unsigned int> &local_dof_indices,
				 SparsityType                    &sparsity_pattern,
				 const bool                       keep_constrained_entries,
				 const Table<2,bool>             &dof_mask,
				 internal::bool2type<true>) const;

				     /**
				      * Internal helper function for
				      * distribute_local_to_global
				      * function.
				      *
				      * Creates a list of affected
				      * global rows for distribution,
				      * including the local rows where
				      * the entries come from.
				      */
    void
    make_sorted_row_list (const std::vector<unsigned int> &local_dof_indices,
			  internals::GlobalRowsFromLocal  &global_rows) const;

				     /**
				      * Internal helper function for
				      * add_entries_local_to_global
				      * function.
				      *
				      * Creates a list of affected
				      * rows for distribution without
				      * any additional information.
				      */
    void
    make_sorted_row_list (const std::vector<unsigned int> &local_dof_indices,
			  std::vector<unsigned int>       &active_dofs) const;

				     /**
				      * Internal helper function for
				      * distribute_local_to_global function.
				      */
    template <typename MatrixType>
    void
    make_sorted_row_list (const FullMatrix<double>        &local_matrix,
			  const std::vector<unsigned int> &local_dof_indices,
			  MatrixType                      &global_matrix,
			  internals::GlobalRowsFromLocal  &global_rows) const;

				     /**
				      * Internal helper function for
				      * add_entries_local_to_global function.
				      */
    template <typename SparsityType>
    void
    make_sorted_row_list (const std::vector<unsigned int> &local_dof_indices,
			  const bool                       keep_constrained_entries,
			  SparsityType                    &sparsity_pattern,
			  std::vector<unsigned int>       &active_dofs) const;

				     /**
				      * Internal helper function for
				      * add_entries_local_to_global function.
				      */
    template <typename SparsityType>
    void
    make_sorted_row_list (const Table<2,bool>             &dof_mask,
			  const std::vector<unsigned int> &local_dof_indices,
			  const bool                       keep_constrained_entries,
			  SparsityType                    &sparsity_pattern,
			  internals::GlobalRowsFromLocal  &global_rows) const;

				     /**
				      * Internal helper function for
				      * distribute_local_to_global function.
				      */
    double
    resolve_vector_entry (const unsigned int                    i,
			  const internals::GlobalRowsFromLocal &global_rows,
			  const Vector<double>                 &local_vector,
			  const std::vector<unsigned int>      &local_dof_indices,
			  const FullMatrix<double>             &local_matrix) const;
};



/* ---------------- template and inline functions ----------------- */

inline
ConstraintMatrix::ConstraintMatrix (const IndexSet &local_constraints)
		 :
		 lines (),
		 local_lines (local_constraints),
		 sorted (false)
{}



inline
ConstraintMatrix::ConstraintMatrix (const ConstraintMatrix &constraint_matrix)
		 :
		 Subscriptor (),
		 lines (constraint_matrix.lines),
		 lines_cache (constraint_matrix.lines_cache),
		 local_lines (constraint_matrix.local_lines),
		 sorted (constraint_matrix.sorted)
{}


inline
void
ConstraintMatrix::add_line (const unsigned int line)
{
  Assert (sorted==false, ExcMatrixIsClosed());

				   // the following can happen when we
				   // compute with distributed meshes
				   // and dof handlers and we
				   // constrain a degree of freedom
				   // whose number we don't have
				   // locally. if we don't abort here
				   // the program will try to allocate
				   // several terabytes of memory to
				   // resize the various arrays below
				   // :-)
  Assert (line != numbers::invalid_unsigned_int,
	  ExcInternalError());
  const unsigned int line_index = calculate_line_index (line);

				   // check whether line already exists; it
				   // may, in which case we can just quit
  if (is_constrained(line))
    return;

				   // if necessary enlarge vector of
				   // existing entries for cache
  if (line_index >= lines_cache.size())
    lines_cache.resize (std::max(2*static_cast<unsigned int>(lines_cache.size()),
				 line_index+1),
			numbers::invalid_unsigned_int);

				   // push a new line to the end of the
				   // list
  lines.push_back (ConstraintLine());
  lines.back().line = line;
  lines.back().inhomogeneity = 0.;
  lines_cache[line_index] = lines.size()-1;
}



inline
void
ConstraintMatrix::add_entry (const unsigned int line,
                             const unsigned int column,
                             const double       value)
{
  Assert (sorted==false, ExcMatrixIsClosed());
  Assert (line != column,
	  ExcMessage ("Can't constrain a degree of freedom to itself"));

				   // if in debug mode, check whether an
				   // entry for this column already
				   // exists and if its the same as
				   // the one entered at present
				   //
				   // in any case: exit the function if an
				   // entry for this column already exists,
				   // since we don't want to enter it twice
  Assert (lines_cache[calculate_line_index(line)] != numbers::invalid_unsigned_int,
	  ExcInternalError());
  ConstraintLine* line_ptr = &lines[lines_cache[calculate_line_index(line)]];
  Assert (line_ptr->line == line, ExcInternalError());
  for (ConstraintLine::Entries::const_iterator
         p=line_ptr->entries.begin();
       p != line_ptr->entries.end(); ++p)
    if (p->first == column)
      {
	Assert (p->second == value,
		ExcEntryAlreadyExists(line, column, p->second, value));
	return;
      }

  line_ptr->entries.push_back (std::make_pair(column,value));
}



inline
void
ConstraintMatrix::set_inhomogeneity (const unsigned int line,
				     const double       value)
{
  ConstraintLine* line_ptr = &lines[lines_cache[calculate_line_index(line)]];
  line_ptr->inhomogeneity = value;
}



inline
bool
ConstraintMatrix::is_constrained (const unsigned int index) const
{
  const unsigned int line_index = calculate_line_index(index);
  return ((line_index < lines_cache.size())
          &&
          (lines_cache[line_index] != numbers::invalid_unsigned_int));
}



inline
bool
ConstraintMatrix::is_inhomogeneously_constrained (const unsigned int index) const
{
  return (is_constrained(index) &&
	  lines[lines_cache[calculate_line_index(index)]].inhomogeneity != 0);
}



inline
const std::vector<std::pair<unsigned int,double> >*
ConstraintMatrix::get_constraint_entries (unsigned int line) const
{
  if (is_constrained(line))
    return &lines[lines_cache[calculate_line_index(line)]].entries;
  else
    return 0;
}



inline unsigned int
ConstraintMatrix::calculate_line_index (const unsigned int line) const
{
				   //IndexSet is unused (serial case)
  if (!local_lines.size())
    return line;

  Assert(local_lines.is_element(line),
	 ExcRowNotStoredHere(line));

  return local_lines.index_within_set(line);
}



inline bool
ConstraintMatrix::can_store_line(unsigned int line_index) const
{
  return !local_lines.size() || local_lines.is_element(line_index);
}



template <class InVector, class OutVector>
inline
void
ConstraintMatrix::
distribute_local_to_global (const InVector                  &local_vector,
                            const std::vector<unsigned int> &local_dof_indices,
                            OutVector                       &global_vector) const
{
  Assert (local_vector.size() == local_dof_indices.size(),
          ExcDimensionMismatch(local_vector.size(), local_dof_indices.size()));
  distribute_local_to_global (local_vector.begin(), local_vector.end(),
			      local_dof_indices.begin(), global_vector);
}



template <typename ForwardIteratorVec, typename ForwardIteratorInd,
          class VectorType>
inline
void ConstraintMatrix::
  distribute_local_to_global (ForwardIteratorVec local_vector_begin,
			      ForwardIteratorVec local_vector_end,
			      ForwardIteratorInd local_indices_begin,
			      VectorType        &global_vector) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  for ( ; local_vector_begin != local_vector_end;
	++local_vector_begin, ++local_indices_begin)
    {
      if (is_constrained(*local_indices_begin) == false)
	global_vector(*local_indices_begin) += *local_vector_begin;
      else
	{
	  const ConstraintLine& position =
	    lines[lines_cache[calculate_line_index(*local_indices_begin)]];
	  for (unsigned int j=0; j<position.entries.size(); ++j)
	    global_vector(position.entries[j].first)
	      += *local_vector_begin * position.entries[j].second;
	}
    }
}



template <typename ForwardIteratorVec, typename ForwardIteratorInd,
          class VectorType>
inline
void ConstraintMatrix::get_dof_values (const VectorType  &global_vector,
				       ForwardIteratorInd local_indices_begin,
				       ForwardIteratorVec local_vector_begin,
				       ForwardIteratorVec local_vector_end) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  for ( ; local_vector_begin != local_vector_end;
	++local_vector_begin, ++local_indices_begin)
    {
      if (is_constrained(*local_indices_begin) == false)
	*local_vector_begin = global_vector(*local_indices_begin);
      else
	{
	  const ConstraintLine & position =
	    lines[lines_cache[calculate_line_index(*local_indices_begin)]];
	  typename VectorType::value_type value = position.inhomogeneity;
	  for (unsigned int j=0; j<position.entries.size(); ++j)
	    value += (global_vector(position.entries[j].first) *
		      position.entries[j].second);
	  *local_vector_begin = value;
	}
    }
}



template <typename MatrixType>
inline
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<double>        &local_matrix,
                            const std::vector<unsigned int> &local_dof_indices,
                            MatrixType                      &global_matrix) const
{
                                   // create a dummy and hand on to the
                                   // function actually implementing this
                                   // feature in the cm.templates.h file.
  Vector<double> dummy(0);
  distribute_local_to_global (local_matrix, dummy, local_dof_indices,
                              global_matrix, dummy,
                              internal::bool2type<IsBlockMatrix<MatrixType>::value>());
}



template <typename MatrixType, typename VectorType>
inline
void
ConstraintMatrix::
distribute_local_to_global (const FullMatrix<double>        &local_matrix,
                            const Vector<double>            &local_vector,
                            const std::vector<unsigned int> &local_dof_indices,
                            MatrixType                      &global_matrix,
                            VectorType                      &global_vector) const
{
                                   // enter the internal function with the
                                   // respective block information set, the
                                   // actual implementation follows in the
                                   // cm.templates.h file.
  distribute_local_to_global (local_matrix, local_vector, local_dof_indices,
                              global_matrix, global_vector,
                              internal::bool2type<IsBlockMatrix<MatrixType>::value>());
}



template <typename SparsityType>
inline
void
ConstraintMatrix::
add_entries_local_to_global (const std::vector<unsigned int> &local_dof_indices,
                             SparsityType                    &sparsity_pattern,
                             const bool                       keep_constrained_entries,
                             const Table<2,bool>             &dof_mask) const
{
                                   // enter the internal function with the
                                   // respective block information set, the
                                   // actual implementation follows in the
                                   // cm.templates.h file.
  add_entries_local_to_global (local_dof_indices, sparsity_pattern,
                               keep_constrained_entries, dof_mask,
                               internal::bool2type<IsBlockMatrix<SparsityType>::value>());
}


DEAL_II_NAMESPACE_CLOSE

#endif

