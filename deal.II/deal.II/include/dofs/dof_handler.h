/*----------------------------   dof.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __dof_H
#define __dof_H
/*----------------------------   dof.h     ---------------------------*/

#include <vector>
#include <map>
#include <base/exceptions.h>
#include <base/forward-declarations.h>
#include <base/smartpointer.h>
#include <basic/forward-declarations.h>
#include <lac/forward-declarations.h>







/**
 * Define some types which differ between the dimensions. This class
 * is analogous to the \Ref{TriaDimensionInfo} class hierarchy.
 * 
 * @see DoFDimensionInfo<1>
 * @see DoFDimensionInfo<2>
 */
template <int dim>
class DoFDimensionInfo;





/**
 * Define some types for the DoF handling in one dimension.
 *
 * The types have the same meaning as those declared in \Ref{TriaDimensionInfo<2>}.
 */
class DoFDimensionInfo<1> {
  public:
    typedef TriaRawIterator<1,DoFCellAccessor<1> >    raw_line_iterator;
    typedef TriaIterator<1,DoFCellAccessor<1> >       line_iterator;
    typedef TriaActiveIterator<1,DoFCellAccessor<1> > active_line_iterator;

    typedef void * raw_quad_iterator;
    typedef void * quad_iterator;
    typedef void * active_quad_iterator;

    typedef void * raw_hex_iterator;
    typedef void * hex_iterator;
    typedef void * active_hex_iterator;

    typedef raw_line_iterator    raw_cell_iterator;
    typedef line_iterator        cell_iterator;
    typedef active_line_iterator active_cell_iterator;

    typedef void * raw_face_iterator;
    typedef void * face_iterator;
    typedef void * active_face_iterator;
};






/**
 * Define some types for the DoF handling in two dimensions.
 *
 * The types have the same meaning as those declared in \Ref{TriaDimensionInfo<2>}.
 */
class DoFDimensionInfo<2> {
  public:
    typedef TriaRawIterator<2,DoFLineAccessor<2,LineAccessor<2> > >    raw_line_iterator;
    typedef TriaIterator<2,DoFLineAccessor<2,LineAccessor<2> > >       line_iterator;
    typedef TriaActiveIterator<2,DoFLineAccessor<2,LineAccessor<2> > > active_line_iterator;
    
    typedef TriaRawIterator<2,DoFCellAccessor<2> >               raw_quad_iterator;
    typedef TriaIterator<2,DoFCellAccessor<2> >                  quad_iterator;
    typedef TriaActiveIterator<2,DoFCellAccessor<2> >            active_quad_iterator;

    typedef void * raw_hex_iterator;
    typedef void * hex_iterator;
    typedef void * active_hex_iterator;

    typedef raw_quad_iterator    raw_cell_iterator;
    typedef quad_iterator        cell_iterator;
    typedef active_quad_iterator active_cell_iterator;

    typedef raw_line_iterator    raw_face_iterator;
    typedef line_iterator        face_iterator;
    typedef active_line_iterator active_face_iterator;    
};



/**
 * Define some types for the DoF handling in two dimensions.
 *
 * The types have the same meaning as those declared in \Ref{TriaDimensionInfo<3>}.
 */
class DoFDimensionInfo<3> {
  public:
    typedef TriaRawIterator<3,DoFLineAccessor<3,LineAccessor<3> > >    raw_line_iterator;
    typedef TriaIterator<3,DoFLineAccessor<3,LineAccessor<3> > >       line_iterator;
    typedef TriaActiveIterator<3,DoFLineAccessor<3,LineAccessor<3> > > active_line_iterator;

    typedef TriaRawIterator<3,DoFQuadAccessor<3,QuadAccessor<3> > >    raw_quad_iterator;
    typedef TriaIterator<3,DoFQuadAccessor<3,QuadAccessor<3> > >       quad_iterator;
    typedef TriaActiveIterator<3,DoFQuadAccessor<3,QuadAccessor<3> > > active_quad_iterator;

    typedef TriaRawIterator<3,DoFCellAccessor<3> >               raw_hex_iterator;
    typedef TriaIterator<3,DoFCellAccessor<3> >                  hex_iterator;
    typedef TriaActiveIterator<3,DoFCellAccessor<3> >            active_hex_iterator;

    typedef raw_hex_iterator    raw_cell_iterator;
    typedef hex_iterator        cell_iterator;
    typedef active_hex_iterator active_cell_iterator;

    typedef raw_quad_iterator    raw_face_iterator;
    typedef quad_iterator        face_iterator;
    typedef active_quad_iterator active_face_iterator;    
};







/**
 * Manage the distribution and numbering of the degrees of freedom for
 * non-multigrid algorithms.
 *
 * We store a list of numbers for each vertex, line, quad, etc
 * denoting the mapping between the degrees of freedom on this object
 * and the global number of this degree of freedom. The numbers refer
 * to the unconstrained degrees of freedom, i.e. constrained degrees
 * of freedom are numbered in the same way as unconstrained ones.
 * This leads to the fact that indices in global vectors and matrices
 * also refer to all degrees of freedom and some kind of condensation
 * is needed to restrict the systems of equations to the unconstrained
 * degrees of freedom only. The actual layout of storage of the indices
 * is described in the \Ref{DoFLevel} class documentation.
 *
 * Additionally, the DoFHandler is able to generate the condensation
 * matrix which connects constrained and unconstrained matrices and
 * vectors.
 *
 * Finally it offers a starting point for the assemblage of the matrices
 * by offering #begin()# and #end()# functions which return iterators
 * to walk on the DoF structures as well as the triangulation data.
 * These iterators work much like those described in the documentation
 * of the #Triangulation# class and of the iterator classes themselved,
 * but offer more functionality than pure triangulation iterators. The
 * order in which dof iterators are presented by the #++# and #--# operators
 * is the same as that for the alike triangulation iterators.
 *
 * This class also provides functions to create the sparsity patterns of
 * global matrices as well as matrices living on (parts of) the boundary
 * and handles some simple forms of transfer of data from one grid to another.
 * 
 * 
 * \subsection{Distribution of indices for degrees of freedom}
 *
 * The degrees of freedom (`dofs') are distributed on the given triangulation
 * by the function #distribute_dofs()#. It gets passed a finite element object
 * describing how many degrees of freedom are located on vertices, lines, etc.
 * It traverses the triangulation cell by cell and numbers the dofs of that
 * cell if not yet numbered. For non-multigrid algorithms, only active cells
 * are considered. Active cells are defined to be those cells which have no
 * children, i.e. they are the most refined ones.
 *
 * Since the triangulation is traversed starting with the cells of the coarsest
 * active level and going to more refined levels, the lowest numbers for dofs
 * are given to the largest cells as well as their bounding lines and vertices,
 * with the dofs of more refined cells getting higher numbers.
 *
 * This numbering implies very large bandwiths of the resulting matrices and
 * is thus vastly suboptimal for some solution algorithms. For this reason,
 * the #DoFRenumbering# class offers the function #renumber_dofs# which reorders
 * the dof numbering according to some scheme. See there for a discussion of
 * the implemented algorithms.
 *
 *
 * \subsection{User defined renumbering schemes}
 *
 * The #DoFRenumbering# class offers a fixed number of renumbering
 * schemes like the Cuthill-McKey scheme. Basically, the function sets
 * up an array in which for each degree of freedom the index is stored
 * which is to be assigned by the renumbering. Using this array, the
 * #renumber_dofs(vector<int>)# function is called, which actually
 * does the change from old DoF indices to the ones given in the
 * array. In some cases, however, a user may want to compute her own
 * renumbering order; in this case, allocate an array with one element
 * per degree of freedom and fill it with the number that the
 * respective degree of freedom shall be assigned. This number may,
 * for example, be obtained by sorting the support points of the
 * degrees of freedom in downwind direction.  Then call the
 * #renumber_dofs(vector<int>)# with the array, which converts old
 * into new degree of freedom indices.
 *
 *
 * \subsection{Setting up sparsity patterns}
 *
 * When assembling system matrices, the entries are usually of the form
 * $a_{ij} = a(\phi_i, \phi_j)$, where $a$ is a bilinear functional, often an
 * integral. When using sparse matrices, we therefore only need to reserve space
 * for those $a_{ij}$ only, which are nonzero, which is the same as to say that
 * that the basis functions $\phi_i$ and $\phi_j$ have a nonempty intersection of
 * their support. Since the support of basis functions is bound only on cells
 * on which they are located or to which they are adjacent, to determine the
 * sparsity pattern it is sufficient to loop over all cells and connect all
 * basis functions on each cell with all other basis functions on that cell.
 * There may be finite elements for which not all basis functions on a cell
 * connect with each other, but no use of this case is made since no examples
 * where this occurs are known to the author.
 *
 * When setting up sparsity patterns for matrices on the boundary, the same
 * procedure is done, except for the fact that the loop only goes over faces
 * on the boundary and the basis functions thereon. It is assumed that all
 * other basis functions on a cell adjacent to the boundary vanish at the
 * boundary itself, except for those which are located on the boundary.
 *
 *
 * \subsection{Boundaries}
 *
 * When projecting the traces of functions to the boundary or parts thereof, one
 * needs to built matrices and vectors with the degrees of freedom on the
 * boundary. What is needed in this case is a numbering of the boundary degrees
 * of freedom, starting from zero on and not considering the degrees of freedom
 * in the interior. The #map_dof_to_boundary_indices# function does exactly
 * this, by providing a vector with as many entries as there are degrees of
 * freedom on the whole domain, with each entry being the number in the
 * numbering of the boundary or #-1# if the dof is not on the boundary. You
 * should always use this function to get the mapping between local (boundary)
 * and the global numbers, for example to build the mass matrix on the
 * boundary, or to get the global index of a degree of freedom if we want to
 * use the solution of the projection onto the boundary to eliminate the
 * boundary degrees of freedom from the global matrix.
 *
 * The algorithm to provide this numbering mapping is simple, but you should
 * not rely on it since it may be changed sometimes: we loop over all faces,
 * check whether it is on the boundary, if so get the global numbers of the
 * degrees of freedom on that face, and for each of these we give a
 * subsequent boundary number if none has already been given to this dof.
 * But it should be emphasized again that you should not try to use this
 * internal knowledge about the used algorithm, you are better off if you
 * just accept the mapping `as is'.
 *
 * Actually, there are two #map_dof_to_boundary_indices# functions, one
 * producing a numbering for all boundary degrees of freedom and one producing
 * a numbering for only parts of the boundary, namely those parts for which
 * the boundary indicator is listed in a set of indicators given to the
 * function. The latter case is needed if, for example, we would only want to
 * project the boundary values for the Dirichlet part of the boundary, not for
 * the other boundary conditions. You then give the function a list of boundary
 * indicators referring to Dirichlet parts on which the projection is to be
 * performed. The parts of the boundary on which you want to project need not
 * be contiguous; however, it is not guaranteed that the indices of each of the
 * boundary parts are continuous, i.e. the indices of degrees of freedom on
 * different parts may be intermixed.
 *
 * Degrees of freedom on the boundary but not on one of the specified boundary
 * parts are given the index #-1#, as if they were in the interior. If no
 * boundary indicator was given or if no face of a cell has a boundary indicator
 * contained in the given list, the vector of new indices consists solely of
 * #-1#s.
 *
 * The question what a degree of freedom on the boundary is, is not so easy.
 * It should really be a degree of freedom of which the respective basis
 * function has nonzero values on the boundary. At least for Lagrange elements
 * this definition is equal to the statement that the off-point of the trial
 * function, i.e. the point where the function assumes its nominal value (for
 * Lagrange elements this is the point where it has the function value #1#), is
 * located on the boundary. We do not check this directly, the criterion is
 * rather defined through the information the finite element class gives: the
 * #FiniteElementBase# class defines the numbers of basis functions per vertex,
 * per line, and so on and the basis functions are numbered after this
 * information; a basis function is to be considered to be on the face of a
 * cell (and thus on the boundary if the cell is at the boundary) according
 * to its belonging to a vertex, line, etc but not to the cell. The finite
 * element uses the same cell-wise numbering so that we can say that if a
 * degree of freedom was numbered as one of the dofs on lines, we assume that
 * it is located on the line. Where the off-point actually is, is a secret of
 * the finite element (well, you can ask it, but we don't do it here) and not
 * relevant in this context.
 *
 *
 * \subsection{Data transfer between grids}
 *
 * {\bf Note: The functionality described in this section has never been tested
 * and is expected to be broken or at least incomplete. Contact the author if
 * you need these functions for more information.}
 *
 * The #DoFHandler# class offers two functions #make_transfer_matrix# which create
 * a matrix to transform the data of one grid to another. The functions assumes the
 * coarsest mesh of the two grids to be the same. However there are few ways to
 * check this (only the number of cells on the coarsest grid is compared). Also,
 * the selected finite element type of the two degree of freedom handler objects
 * must be the same.
 *
 * The algorithm goes recursively from the coarse mesh cells to their children
 * until the grids differ at this level. It then tries to prolong or restrict the
 * old cell(s) to the new cell(s) and makes up a matrix of these prolongations and
 * restrictions. This matrix multiplied with a vector on the old grid yields an
 * approximation of the projection of the function on the old grid to the new one.
 *
 * Building and using the transfer matrix is usually quite an expensive operation,
 * since we have to perform two runs over all cells (one for building the sparsity
 * structure, one to build the entries) and because of the memory consumption.
 * It may, however, pay if you have many
 * equations, since then the entries in the matrix can be considered as block
 * entries which are then applied to all function values at a given degree of
 * freedom.
 *
 * To build the matrix, you first have to call
 * #make_transfer_matrix (old_dof_object, sparsity_pattern);#, then create a
 * sparse matrix out of this pattern, e.g. by #SparseMatrix<double> m(sparsity_pattern);#
 * and finally give this to the second run:
 * #make_transfer_matrix (old_dof_object, m);#. The spasity pattern created
 * by the first run is automatically compressed.
 *
 * When creating the #SparseMatrixStruct# sparsity pattern, you have to give the
 * dimension and the maximum number of entries per row. Obviously the image
 * dimension is the number of dofs on the new grid (you can get this using the
 * #n_dofs()# function), while the range dimension is the number of dofs on the
 * old grid. The maximum number of entries per row is determined by the maximum
 * number of levels $d$ which we have to cross upon transferring from one cell to
 * another (presently, transfer of one cell is only possible for #d=0,1#, i.e.
 * the two cells match or one is refined once more than the other, the
 * number of degrees of freedom per vertex $d_v$, those on lines $d_l$, those
 * on quads $d_q$ and the number of subcells a cell is
 * refined to, which is $2^{dim}$. The maximum number of entries per row in one
 * dimension is then given by $(2*d_l+d_v)*2+1$ if $d=1$. For example, a one
 * dimensional linear element would need two entries per row.
 * In two dimensions, the maxmimum number is $(4*d_q+12*d_l+5*d_v)*4+1$ if $d=1$.
 * You can get these numbers by drawing little pictures and counting, there is
 * no mystique behind this. You can also get the right number by calling the
 * #max_transfer_entries (max_level_difference)# function. The actual number
 * depends on the finite element selected and may be much less, especially in
 * higher dimensions.
 *
 *
 *
 * @author Wolfgang Bangerth, 1998 */
template <int dim>
class DoFHandler
  :
  public Subscriptor,
  public DoFDimensionInfo<dim>
{
  public:
    typedef typename DoFDimensionInfo<dim>::raw_line_iterator raw_line_iterator;
    typedef typename DoFDimensionInfo<dim>::line_iterator line_iterator;
    typedef typename DoFDimensionInfo<dim>::active_line_iterator active_line_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_quad_iterator raw_quad_iterator;
    typedef typename DoFDimensionInfo<dim>::quad_iterator quad_iterator;
    typedef typename DoFDimensionInfo<dim>::active_quad_iterator active_quad_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_hex_iterator raw_hex_iterator;
    typedef typename DoFDimensionInfo<dim>::hex_iterator hex_iterator;
    typedef typename DoFDimensionInfo<dim>::active_hex_iterator active_hex_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_cell_iterator raw_cell_iterator;
    typedef typename DoFDimensionInfo<dim>::cell_iterator cell_iterator;
    typedef typename DoFDimensionInfo<dim>::active_cell_iterator active_cell_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_face_iterator raw_face_iterator;
    typedef typename DoFDimensionInfo<dim>::face_iterator face_iterator;
    typedef typename DoFDimensionInfo<dim>::active_face_iterator active_face_iterator;

				     /**
				      *	Declare a data type which denotes a
				      *	mapping between a boundary indicator
				      *	and the function denoting the boundary
				      *	values on this part of the boundary.
				      *	Only one boundary function may be given
				      *	for each boundary indicator, which is
				      *	guaranteed by the #map# data type.
				      */
    typedef map<unsigned char,const Function<dim>*> FunctionMap;

    
				     /**
				      * Constructor. Take #tria# as the
				      * triangulation to work on.
				      */
    DoFHandler (Triangulation<dim> *tria);

				     /**
				      * Destructor.
				      */
    virtual ~DoFHandler ();
    
				     /**
				      * Go through the triangulation and
				      * distribute the degrees of freedoms
				      * needed for the given finite element
				      * according to the given distribution
				      * method.
				      *
				      * A pointer of the transferred finite
				      * element is stored. Therefore, the
				      * lifetime of the finite element object
				      * shall be longer than that of this
				      * object. If you don't want this
				      * behaviour, you may want to call
				      * the #clear# member function which also
				      * releases the lock of this object to the
				      * finite element.
				      *
				      * This function uses the user flags of the
				      * triangulation object, so make sure you
				      * don't need them after calling this
				      * function, or if so store them.
				      */
    virtual void distribute_dofs (const FiniteElement<dim> &);

				     /**
				      * Clear all data of this object and
				      * especially delete the lock this object
				      * has to the finite element used the last
				      * time when #distribute_dofs# was called.
				      */
    virtual void clear ();
    
				     /**
				      * Actually do the renumbering based on
				      * a list of new dof numbers for all the
				      * dofs.
				      *
				      * #new_numbers# is an array of integers
				      * with size equal to the number of dofs
				      * on the present grid. It stores the new
				      * indices after renumbering in the
				      * order of the old indices.
				      *
				      * This function is called by the
				      * #renumber_dofs# function after computing
				      * the ordering of the degrees of freedom.
				      * However, you can call this function
				      * yourself, which is necessary if a user
				      * wants to implement an ordering scheme
				      * herself, for example downwind numbering.
				      */
    void renumber_dofs (const vector<int> &new_numbers);

				     /**
				      * Make up the constraints which
				      * is result from the use of hanging
				      * nodes. The object into which these
				      * are inserted is later
				      * used to condensate the global
				      * system matrices and to prolong
				      * the solution vectors from the true
				      * degrees of freedom also to the
				      * constraint nodes.
				      *
				      * Since this method does not make sense in
				      * one dimension, the function returns
				      * immediately. The object is not cleared
				      * before use, so you should make sure
				      * it containts only constraints you still
				      * want; otherwise call the #clear#
				      * function.
				      *
				      * To condense a given sparsity pattern,
				      * use #ConstraintMatrix::condense#.
				      * Before doing so, you need to close
				      * the constraint object, which must be
				      * done after all constraints are entered.
				      * This function does not close the object
				      * since you may want to enter other
				      * constraints later on yourself.
				      *
				      * This function uses the user flags for
				      * the faces. If you need the user flags,
				      * store them beforehand.
				      */
    void make_hanging_node_constraints (ConstraintMatrix &) const;


				     /**
				      * Write the sparsity structure of the
				      * full matrix including constrained
				      * degrees of freedom into the
				      * matrix structure. The sparsity pattern
				      * does not include entries introduced by
				      * the elimination of constrained nodes.
				      * The sparsity
				      * pattern is not compressed, since if
				      * you want to call
				      * #ConstraintMatrix::condense(1)#
				      * afterwards, new entries have to be
				      * added. However, if you don't want to call
				      * #ConstraintMatrix::condense(1)#, you
				      * have to compress the matrix yourself,
				      * using #SparseMatrixStruct::compress()#.
				      */
    void make_sparsity_pattern (SparseMatrixStruct &) const; 

				     /**
				      * This function does mistly the same as
				      * the above one, but it is specialized for
				      * mixed finite elements and allows to
				      * specify which variables couple in which
				      * equation. For example, if wanted to solve
				      * the Stokes equations,
				      * \begin{verbatim}
				      *    -\Delta \vec u + \nabla p = 0,
				      *    \div u                    = 0
				      * \end{verbatim}
				      * in, say, two space dimensions, using
				      * nonstabilized Q2/Q2/Q1 mixed elements
				      * (i.e. using the #FESystem# class), then
				      * you don't want all degrees of freedom
				      * to couple in each equation. You rather
				      * may want to give the following pattern
				      * of couplings:
				      * \begin{verbatim}
				      *    1  0  1
				      *    0  1  1
				      *    1  1  0
				      * \end{verbatim}
				      * where "1" indicates that two variables
				      * (i.e. components of the #FESystem#)
				      * couple in the respective equation, and
				      * a "0" means no coupling, in which case
				      * it is not necessary to allocate space
				      * in the matrix structure. Obviously, the
				      * mask refers to components of the
				      * composed #FESystem#, rather than to the
				      * degrees of freedom contained in there.
				      *
				      * This function is designed to accept
				      * a mask, like the one shown above,
				      * through the #mask# parameter, which
				      * contains boolean values. It builds
				      * the matrix structure just like the
				      * previous function, but does not create
				      * elements if not specified by the mask.
				      */
    void make_sparsity_pattern (const vector<vector<bool> > &mask,
				SparseMatrixStruct          &) const;
    
    				     /**
				      * Write the sparsity structure of the
				      * matrix composed of the basis functions
				      * on the boundary into the
				      * matrix structure. The sparsity pattern
				      * does not include entries introduced by
				      * the elimination of constrained nodes.
				      * The sparsity
				      * pattern is not compressed, since if
				      * you want to call
				      * #ConstraintMatrix::condense(1)#
				      * afterwards, new entries have to be
				      * added. However, if you don't want to call
				      * #ConstraintMatrix::condense(1)#, you
				      * have to compress the matrix yourself,
				      * using #SparseMatrixStruct::compress()#.
				      *
				      * Since this function is obviously useless
				      * in one spatial dimension, it is not
				      * implemented.
				      */
    void make_boundary_sparsity_pattern (const vector<int> &dof_to_boundary_mapping,
					 SparseMatrixStruct &) const; 

				     /**
				      * Write the sparsity structure of the
				      * matrix composed of the basis functions
				      * on the boundary into the
				      * matrix structure. In contrast to the
				      * previous function, only those parts
				      * of the boundary are considered of which
				      * the boundary indicator is listed in the
				      * set of numbers passed to this function.
				      *
				      * In fact, rather than a #set# of boundary
				      * indicators, a #map# needs to be passed,
				      * since most of the functions handling with
				      * boundary indicators take a mapping of
				      * boundary indicators and the respective
				      * boundary functions. The boundary function,
				      * however, is ignored in this function.
				      * If you have no functions at hand, but only
				      * the boundary indicators, set the function
				      * pointers to null pointers.
				      *
				      * Since this function is obviously useless
				      * in one spatial dimension, it is not
				      * implemented.
				      */
    void make_boundary_sparsity_pattern (const FunctionMap  &boundary_indicators,
					 const vector<int>  &dof_to_boundary_mapping,
					 SparseMatrixStruct &sparsity) const; 

    
				     /**
				      * Make up the transfer matrix which
				      * transforms the data vectors from one
				      * triangulation to the present one.
				      * You have to pass references to the old
				      * dof handler object and to a matrix
				      * sparsity object. This function therefore
				      * only makes up the sparsity pattern.
				      *
				      * The given sparsity pattern is
				      * compressed at the end.
				      *
				      * In the matrix, row indices belong to
				      * new dof numbers, column indices to the
				      * ones on the old grid. Therefore,
				      * multiplying this matrix by a vector
				      * of the old grid yields the vector on
				      * the new one.
				      *
				      * For more details see the general
				      * documentation for this class.
				      */
    void make_transfer_matrix (const DoFHandler<dim> &transfer_from,
			       SparseMatrixStruct    &transfer_pattern) const;

				     /**
				      * Make up the transfer matrix which
				      * transforms the data vectors from one
				      * triangulation to the present one.
				      * You have to pass references to the old
				      * dof handler object and to a matrix
				      * object. This function therefore
				      * builds the matrix itself
				      *
				      * The matrix object should be
				      * associated with the sparsity pattern
				      * constructed by the other
				      * #make_transfer_matrix# object.
				      *
				      * For more details see the general
				      * documentation for this class.
				      */
    void make_transfer_matrix (const DoFHandler<dim> &transfer_from,
			       SparseMatrix<double>  &transfer_matrix) const;

				     /**
				      * Return the maximum number of
				      * degrees of freedom a degree of freedom
				      * in the given triangulation with the
				      * given finite element may couple with.
				      * This is the maximum number of entries
				      * per line in the system matrix; this
				      * information can therefore be used upon
				      * construction of the #SparseMatrixStruct#
				      * object.
				      *
				      * The returned number is not really the
				      * maximum number but an estimate based
				      * on the finite element and the maximum
				      * number of cells meeting at a vertex.
				      * The number holds for the constrained
				      * matrix also.
				      *
				      * The determination of the number of
				      * couplings can be done by simple
				      * picture drawing. An example can be
				      * found in the implementation of this
				      * function.
				      */
    unsigned int max_couplings_between_dofs () const;

				     /**
				      * Return the number of degrees of freedom
				      * located on the boundary another dof on
				      * the boundary can couple with.
				      *
				      * The number is the same as for
				      * #max_coupling_between_dofs# in one
				      * dimension less.
				      */
    unsigned int max_couplings_between_boundary_dofs () const;
    
				     /**
				      * Return the maximum number of entries
				      * a row in a transfer matrix may contain
				      * if any two cells of which the dofs are
				      * to be transferred differ in refinement
				      * level at most by #max_level_diff#.
				      * It is assumed that the finite element
				      * selected by the last call to
				      * #distribute_dofs# is used also for
				      * the transfer process.
				      */
    unsigned int max_transfer_entries (const unsigned int max_level_diff) const;

				     /**
				      * Take a vector of values which live on
				      * cells (e.g. an error per cell) and
				      * distribute it to the dofs in such a
				      * way that a finite element field results,
				      * which can then be further processed,
				      * e.g. for output. You should note that
				      * the resulting field will not be
				      * continuous at hanging nodes. This can,
				      * however, easily be arranged by calling
				      * the appropraite #distribute# function
				      * of a #ConstraintMatrix# object created
				      * for this #DoFHandler# object.
				      *
				      * It is assumed that the number of
				      * elements in #cell_data# equals the
				      * number of active cells. The size of
				      * #dof_data# is adjusted to the right
				      * size.
				      *
				      * Note that the input vector may be
				      * a vector of any data type as long
				      * as it is convertible to #double#.
				      * The output vector, being a data
				      * vector on the grid, always consists
				      * of elements of type #double#.
				      *
				      * In case the finite element used by
				      * this DoFHandler consists of more than
				      * one component, you should give which
				      * component in the output vector should
				      * be used to store the finite element
				      * field in; the default is zero (no other
				      * value is allowed if the finite element
				      * consists only of one component). All
				      * other components of the vector remain
				      * untouched, i.e. their contents are
				      * not changed.
				      *
				      * It is assumed that the output vector
				      * #dof_data# already has the right size,
				      * i.e. #n_dofs()# elements.
				      */
    template <typename Number>
    void distribute_cell_to_dof_vector (const Vector<Number> &cell_data,
					Vector<double>       &dof_data,
					const unsigned int    component = 0) const;

				     /**
				      * Create a mapping from degree of freedom
				      * indices to the index of that degree
				      * of freedom on the boundary. After this
				      * operation, #mapping[dof]# gives the
				      * index of the the degree of freedom with
				      * global number #dof# in the list of
				      * degrees of freedom on the boundary.
				      * If the degree of freedom requested is
				      * not on the boundary, the value of
				      * #mapping[dof]# is #-1#. This function is
				      * mainly used when setting up matrices and
				      * vectors on the boundary from the trial
				      * functions, which have global numbers,
				      * while the matrices and vectors use
				      * numbers of the trial functions local
				      * to the boundary.
				      *
				      * Prior content of #mapping# is deleted.
				      *
				      * This function is not implemented for
				      * one dimension. See the general doc
				      * for more information on boundary
				      * treatment.
				      */
    void map_dof_to_boundary_indices (vector<int> &mapping) const;

				     /**
				      * Same as the previous function, except
				      * that only selected parts of the
				      * boundary are considered.
				      *
				      * See the general doc of this class for
				      * more information.
				      */
    void map_dof_to_boundary_indices (const FunctionMap &boundary_indicators,
				      vector<int> &mapping) const;

				     /**
				      *  @name Cell iterator functions
				      */
				     /*@{*/
				     /**
				      *  Return iterator to the first cell, used
				      *  or not, on level #level#. If a level
				      *  has no cells, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls #begin_raw_line#
				      *  in 1D and #begin_raw_quad# in 2D.
				      */
    raw_cell_iterator    begin_raw   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used cell
				      *  on level #level#.
				      *
				      *  This function calls #begin_line#
				      *  in 1D and #begin_quad# in 2D.
				      */
    cell_iterator        begin       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  cell on level #level#.
				      *
				      *  This function calls #begin_active_line#
				      *  in 1D and #begin_active_quad# in 2D.
				      */
    active_cell_iterator begin_active(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls #end_line#
				      *  in 1D and #end_quad# in 2D.
				      */
    raw_cell_iterator    end () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    cell_iterator        end (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_cell_iterator    end_raw (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_cell_iterator end_active (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last cell, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 1D and #last_raw_quad# in 2D.
				      */
    raw_cell_iterator    last_raw () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  cell of the level #level#, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 1D and #last_raw_quad# in 2D.
				      */
    raw_cell_iterator    last_raw (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell.
				      *
				      *  This function calls #last_line#
				      *  in 1D and #last_quad# in 2D.
				      */
    cell_iterator        last () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used cell on level #level#.
				      *
				      *  This function calls #last_line#
				      *  in 1D and #last_quad# in 2D.
				      */
    cell_iterator        last (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active cell.
				      *
				      *  This function calls #last_active_line#
				      *  in 1D and #last_active_quad# in 2D.
				      */
    active_cell_iterator last_active () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active cell on level #level#.
				      *
				      *  This function calls #last_active_line#
				      *  in 1D and #last_active_quad# in 2D.
				      */
    active_cell_iterator last_active (const unsigned int level) const;
				     //@}

    				     /*---------------------------------------*/

    				     /**
				      *  @name Face iterator functions
				      */
				     /*@{*/
				     /**
				      *  Return iterator to the first face, used
				      *  or not, on level #level#. If a level
				      *  has no faces, a past-the-end iterator
				      *  is returned.
				      *
				      *  This function calls #begin_raw_line#
				      *  in 2D and #begin_raw_quad# in 3D.
				      */
    raw_face_iterator    begin_raw_face   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used face
				      *  on level #level#.
				      *
				      *  This function calls #begin_line#
				      *  in 2D and #begin_quad# in 3D.
				      */
    face_iterator        begin_face       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  face on level #level#.
				      *
				      *  This function calls #begin_active_line#
				      *  in 2D and #begin_active_quad# in 3D.
				      */
    active_face_iterator begin_active_face(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      *
				      *  This function calls #end_line#
				      *  in 2D and #end_quad# in 3D.
				      */
    raw_face_iterator    end_face () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    face_iterator        end_face (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_face_iterator    end_raw_face (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_face_iterator end_active_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last face, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 2D and #last_raw_quad# in 3D.
				      */
    raw_face_iterator    last_raw_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  face of the level #level#, used or not.
				      *
				      *  This function calls #last_raw_line#
				      *  in 2D and #last_raw_quad# in 3D.
				      */
    raw_face_iterator    last_raw_face (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face.
				      *
				      *  This function calls #last_line#
				      *  in 2D and #last_quad# in 3D.
				      */
    face_iterator        last_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used face on level #level#.
				      *
				      *  This function calls #last_line#
				      *  in 2D and #last_quad# in 3D.
				      */
    face_iterator        last_face (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active face.
				      *
				      *  This function calls #last_active_line#
				      *  in 2D and #last_active_quad# in 3D.
				      */
    active_face_iterator last_active_face () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active face on level #level#.
				      *
				      *  This function calls #last_active_line#
				      *  in 2D and #last_active_quad# in 3D.
				      */
    active_face_iterator last_active_face (const unsigned int level) const;
				     //@}

    
				     /*---------------------------------------*/

				     /**
				      *  @name Line iterator functions
				      */
				     /*@{*/
				     /**
				      *  Return iterator to the first line, used
				      *  or not, on level #level#. If a level
				      *  has no lines, a past-the-end iterator
				      *  is returned.
				      */
    raw_line_iterator    begin_raw_line   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used line
				      *  on level #level#.
				      */
    line_iterator        begin_line       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  line on level #level#.
				      */
    active_line_iterator begin_active_line(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_line_iterator    end_line () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    line_iterator        end_line (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_line_iterator    end_raw_line (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_line_iterator end_active_line (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last line, used or not.
				      */
    raw_line_iterator    last_raw_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  line of the level #level#, used or not.

				      */
    raw_line_iterator    last_raw_line (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line.
				      */
    line_iterator        last_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line on level #level#.
				      */
    line_iterator        last_line (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active line.
				      */
    active_line_iterator last_active_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active line on level #level#.
				      */
    active_line_iterator last_active_line (const unsigned int level) const;
				     /*@}*/	  

				     /*---------------------------------------*/

				     /**
				      *  @name Quad iterator functions*/
    				     /*@{
				      */
    				     /**
				      *  Return iterator to the first quad, used
				      *  or not, on level #level#. If a level
				      *  has no quads, a past-the-end iterator
				      *  is returned.
				      */
    raw_quad_iterator    begin_raw_quad   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used quad
				      *  on level #level#.
				      */
    quad_iterator        begin_quad       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  quad on level #level#.
				      */
    active_quad_iterator begin_active_quad(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_quad_iterator    end_quad () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    quad_iterator        end_quad (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_quad_iterator    end_raw_quad (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_quad_iterator end_active_quad (const unsigned int level) const;


				     /**
				      *  Return an iterator pointing to the
				      *  last quad, used or not.
				      */
    raw_quad_iterator    last_raw_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  quad of the level #level#, used or not.

				      */
    raw_quad_iterator    last_raw_quad (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad.
				      */
    quad_iterator        last_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad on level #level#.
				      */
    quad_iterator        last_quad (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active quad.
				      */
    active_quad_iterator last_active_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active quad on level #level#.
				      */
    active_quad_iterator last_active_quad (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/

				     /**
				      *  @name Hex iterator functions*/
    				     /*@{
				      */
    				     /**
				      *  Return iterator to the first hex, used
				      *  or not, on level #level#. If a level
				      *  has no hexs, a past-the-end iterator
				      *  is returned.
				      */
    raw_hex_iterator
    begin_raw_hex   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used hex
				      *  on level #level#.
				      */
    hex_iterator
    begin_hex       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  hex on level #level#.
				      */
    active_hex_iterator
    begin_active_hex(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_hex_iterator
    end_hex () const;

				     /**
				      * Return an iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    hex_iterator        end_hex (const unsigned int level) const;
    
				     /**
				      * Return a raw iterator which is the first
				      * iterator not on level. If #level# is
				      * the last level, then this returns
				      * #end()#.
				      */
    raw_hex_iterator    end_raw_hex (const unsigned int level) const;

    				     /**
				      * Return an active iterator which is the
				      * first iterator not on level. If #level#
				      * is the last level, then this returns
				      * #end()#.
				      */
    active_hex_iterator end_active_hex (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the
				      *  last hex, used or not.
				      */
    raw_hex_iterator
    last_raw_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  hex of the level #level#, used or not.

				      */
    raw_hex_iterator
    last_raw_hex (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used hex.
				      */
    hex_iterator
    last_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used hex on level #level#.
				      */
    hex_iterator
    last_hex (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active hex.
				      */
    active_hex_iterator
    last_active_hex () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active hex on level #level#.
				      */
    active_hex_iterator
    last_active_hex (const unsigned int level) const;
				     /*@}*/

				     /*---------------------------------------*/


				     /**
				      * Return number of degrees of freedom.
				      * Included in this number are those
				      * DoFs which are constrained by
				      * hanging nodes.
				      */
    unsigned int n_dofs () const;

				     /**
				      * Return the number of degrees of freedom
				      * located on the boundary.
				      */
    unsigned int n_boundary_dofs () const;

    				     /**
				      * Return the number of degrees of freedom
				      * located on those parts of the boundary
				      * which have a boundary indicator listed
				      * in the given set. The reason that a
				      * #map# rather than a #set# is used is the
				      * same as descibed in the section on the
				      * #make_boundary_sparsity_pattern# function.
				      */
    unsigned int n_boundary_dofs (const FunctionMap &boundary_indicators) const;

				     /**
				      * Return a constant reference to the
				      * selected finite element object. This
				      * function is inline, so it should
				      * be reasonably fast.
				      */
    const FiniteElement<dim> & get_fe () const;

				     /**
				      * Return a constant reference to the
				      * triangulation underlying this object.
				      */
    const Triangulation<dim> & get_tria () const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidTriangulation);
				     /**
				      * Exception
				      */
    DeclException2 (ExcDifferentDimensions,
		    int, int,
		    << "One dimension of the matrices is differing: "
		    << arg1 << " vs. " << arg2);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoFESelected);
    				     /**
				      * Exception
				      */
    DeclException0 (ExcRenumberingIncomplete);
				     /**
				      * Exception
				      */
    DeclException0 (ExcGridsDoNotMatch);
				     /**
				      * Exception
				      */
    DeclException0 (ExcOnlyOnelevelTransferImplemented);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidBoundaryIndicator);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException1 (ExcMatrixHasWrongSize,
		    int,
		    << "The matrix has the wrong dimension " << arg1);
				     /**
				      *  Exception
				      */
    DeclException0 (ExcFunctionNotUseful);
				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongSize,
		    int, int,
		    << "The dimension " << arg1 << " of the vector is wrong. "
		    << "It should be " << arg2);
				     /**
				      * Exception
				      */
    DeclException1 (ExcNewNumbersNotConsecutive,
		    int,
		    << "The given list of new dof indices is not consecutive: "
		    << "the index " << arg1 << " does not exist.");
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidMaskDimension,
		    int, int,
		    << "The dimension of the mask " << arg1 << " does not match"
		    << " the number of components in the finite element object.");    
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidComponent,
		    int, int,
		    << "The component you gave (" << arg1 << ") "
		    << "is invalid with respect to the number "
		    << "of components in the finite element "
		    << "(" << arg2 << ")");
    
  protected:
    
				     /**
				      * Address of the triangulation to
				      * work on.
				      */
    Triangulation<dim> * const tria;

				     /**
				      * Store a pointer to the finite element
				      * given latest for the distribution of
				      * dofs. In order to avoid destruction of
				      * the object before the lifetime of
				      * the DoF handler, we subscribe to
				      * the finite element object. To unlock
				      * the FE before the end of the lifetime
				      * of this DoF handler, use the #clear()#
				      * function (this clears all data of
				      * this object as well, though).
				      */
    SmartPointer<const FiniteElement<dim> > selected_fe;

  private:
				     /**
				      * Reserve enough space in the 
				      * #levels[]# objects to store the
				      * numbers of the degrees of freedom
				      * needed for the given element. The
				      * given element is that one which
				      * was selected when calling
				      * #distribute_dofs# the last time.
				      */
    void reserve_space ();

				     /**
				      * Free all used memory.
				      */
    void clear_space ();
    
				     /**
				      * Distribute dofs on the given cell,
				      * with new dofs starting with index
				      * #next_free_dof#. Return the next
				      * unused index number. The finite
				      * element used is the one given to
				      * #distribute_dofs#, which is copied
				      * to #selected_fe#.
				      *
				      * This function is excluded from the
				      * #distribute_dofs# function since
				      * it can not be implemented dimension
				      * independent.
				      */
    unsigned int distribute_dofs_on_cell (active_cell_iterator &cell,
					  unsigned int next_free_dof);

				     /**
				      * Make up part of the sparsity pattern of
				      * the transfer matrix by looking at the
				      * two cells given.
				      */
    void transfer_cell (const cell_iterator &old_cell,
			const cell_iterator &new_cell,
			SparseMatrixStruct  &transfer_pattern) const;

    				     /**
				      * Make up part of the transfer matrix by
				      * looking at the two cells given.
				      */
    void transfer_cell (const cell_iterator  &old_cell,
			const cell_iterator  &new_cell,
			SparseMatrix<double> &transfer_matrix) const;

				     /**
				      * Space to store the DoF numbers for the
				      * different levels. Analogous to the
				      * #levels[]# tree of the \Ref{Triangulation}
				      * objects.
				      */
    vector<DoFLevel<dim>*>    levels;

				     /**
				      * Store the number of dofs created last
				      * time.
				      */
    unsigned int              used_dofs;

				     /**
				      * Array to store the indices for degrees
				      * of freedom located at vertices.
				      */
    vector<int>               vertex_dofs;
    
    
    friend class DoFLineAccessor<dim, LineAccessor<dim> >;
    friend class DoFLineAccessor<dim, CellAccessor<dim> >;
    friend class DoFQuadAccessor<dim, QuadAccessor<dim> >;
    friend class DoFQuadAccessor<dim, CellAccessor<dim> >;
    friend class DoFHexAccessor<dim, HexAccessor<dim> >;
    friend class DoFHexAccessor<dim, CellAccessor<dim> >;
};






/* ----------------------- Inline functions ---------------------------------- */


template <int dim>
inline
const FiniteElement<dim> & DoFHandler<dim>::get_fe () const {
  return *selected_fe;
};






/*----------------------------   dof.h     ---------------------------*/
/* end of #ifndef __dof_H */
#endif
/*----------------------------   dof.h     ---------------------------*/
