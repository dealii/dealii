/*----------------------------   dof.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dof_H
#define __dof_H
/*----------------------------   dof.h     ---------------------------*/

#include <vector>
#include <base/exceptions.h>




template <int dim> class TriaAccessor;
template <int dim> class LineAccessor;
template <int dim> class QuadAccessor;
template <int dim> class CellAccessor;

template <int dim, class BaseClass> class DoFLineAccessor;
template <int dim, class BaseClass> class DoFQuadAccessor;
template <int dim> class DoFCellAccessor;


template <int dim, class Accessor> class TriaRawIterator;
template <int dim, class Accessor> class TriaIterator;
template <int dim, class Accessor> class TriaActiveIterator;

template <int dim> class Triangulation;

template <int dim> class FiniteElementBase;

class dVector;
class dSMatrix;
class dSMatrixStruct;
class ConstraintMatrix;




/**
  Store the indices of the degrees of freedom which are located on the lines.
  Declare it to have a template parameter, but do not actually declare
  other types than those explicitely instantiated.
  */
template <int N>
class DoFLevel;





/**
  Store the indices of the degrees of freedom which are located on the lines.

  \subsection{Information for all #DoFLevel# classes}

  The #DoFLevel<N># classes 
  store the global indices of the degrees of freedom for each cell on a
  certain level. The index or number of a degree of freedom is the zero-based
  index of the according value in the solution vector and the row and column
  index in the global matrix or the multigrid matrix for this level. These
  indices refer to the unconstrained vectors and matrices, where we have not
  taken account of the constraints introduced by hanging nodes. If more than
  one value corresponds to one basis function, for example for vector equations
  where the solution is vector valued and thus has several degrees of freedom
  for each basis function, we nonetheless store only one index. This can then
  be viewed as the index into a block vector, where each block contains the
  different values according to a degree of freedom. It is left to the derived
  classes, whether the values in a block are stored consecutively or distributed
  (e.g. if the solution function is $u=(u_1, u_2)$, we could store the values
  in the solution vector like
  $\ldots, u_1^m, u_2^m, u_1^{m+1}, u_2^{m+1},\ldots$ with $m$ denoting the
  $m$th basis function, or $\ldots, u_1^m, u_1^{m+1}, u_1^{m+2}, \ldots,
  u_2^m, u_2^{m+1}, u_2^{m+2}, \ldots$, respectively). Likewise, the
  constraint matrix returned by #DoFHandler::make_constraint_matrix ()# is then
  to be understood as a block matrix.

  The storage format of the degrees of freedom indices (short: DoF indices) is
  somewhat like a mirror of the data structures of the triangulation classes.
  There is a hierarchy of #DoFLevel<dim># classes for the different dimensions
  which have objects named #line_dofs#, #quad_dofs# and so on, in which the
  indices of DoFs located on lines and quads, respectively, are stored. The
  indices are stored levelwise. The layout in
  these arrays is as follows: if for a selected finite element (use
  #DoFHandler::distribute_dofs()# to select a finite element) the number of
  DoFs on each line (without those in the vertices) is #N#, then the length
  of the #line_dofs# array is #N# times the number of lines on this level. The
  DoF indices for the #i#th line are at the positions #N*i...(N+1)*i-1#.

  The DoF indices for vertices are not stored this way, since they need
  different treatment in multigrid environments. If no multigrid is used, the
  indices are stored in the #vertex_dofs# array of the #DoFHandler# class.
  */
class DoFLevel<1> {
  public:
				     /**
				      * Store the global indices of the degrees
				      * of freedom. See \Ref{DoFLevel} for
				      * detailed information.
				      */
    vector<int> line_dofs;
};




/**
  Store the indices of the degrees of freedom which are located on quads.
  See \Ref{DoFLevel<1>} for more information.
  */
class DoFLevel<2> : public DoFLevel<1> {
  public:
				     /**
				      * Store the global indices of the degrees
				      * of freedom. See \Ref{DoFLevel} for
				      * detailed information.
				      */
    vector<int> quad_dofs;
};






/**
  Define some types which differ between the dimensions. This class
  is analogous to the \Ref{TriaDimensionInfo} class hierarchy.
  
  @see DoFDimensionInfo<1>
  @see DoFDimensionInfo<2>
  */
template <int dim>
class DoFDimensionInfo;





/**
  Define some types for the DoF handling in one dimension.

  The types have the same meaning as those declared in \Ref{TriaDimensionInfo<2>}.
  */
class DoFDimensionInfo<1> {
  public:
    typedef TriaRawIterator<1,DoFCellAccessor<1> >    raw_line_iterator;
    typedef TriaIterator<1,DoFCellAccessor<1> >       line_iterator;
    typedef TriaActiveIterator<1,DoFCellAccessor<1> > active_line_iterator;

    typedef void * raw_quad_iterator;
    typedef void * quad_iterator;
    typedef void * active_quad_iterator;

    typedef raw_line_iterator    raw_cell_iterator;
    typedef line_iterator        cell_iterator;
    typedef active_line_iterator active_cell_iterator;

    typedef void * raw_face_iterator;
    typedef void * face_iterator;
    typedef void * active_face_iterator;
};






/**
  Define some types for the DoF handling in two dimensions.

  The types have the same meaning as those declared in \Ref{TriaDimensionInfo<2>}.
  */
class DoFDimensionInfo<2> {
  public:
    typedef TriaRawIterator<2,DoFLineAccessor<2,LineAccessor<2> > >    raw_line_iterator;
    typedef TriaIterator<2,DoFLineAccessor<2,LineAccessor<2> > >       line_iterator;
    typedef TriaActiveIterator<2,DoFLineAccessor<2,LineAccessor<2> > > active_line_iterator;
    
    typedef TriaRawIterator<2,DoFCellAccessor<2> >               raw_quad_iterator;
    typedef TriaIterator<2,DoFCellAccessor<2> >                  quad_iterator;
    typedef TriaActiveIterator<2,DoFCellAccessor<2> >            active_quad_iterator;

    typedef raw_quad_iterator    raw_cell_iterator;
    typedef quad_iterator        cell_iterator;
    typedef active_quad_iterator active_cell_iterator;

    typedef raw_line_iterator    raw_face_iterator;
    typedef line_iterator        face_iterator;
    typedef active_line_iterator active_face_iterator;    
};






/**
  Give names to the different possibilities of renumbering the degrees
  of freedom.

  \begin{itemize}
  \item #Cuthill_McKee# and #reverse_Cuthill_McKee# traverse the triangulation
    in a diagonal, advancing front like method and produce matrices with an
    almost minimal bandwidth.
  \item #reverse_Cuthill_McKey# does the same thing, but numbers the dofs in
    the reverse order.
  \end{itemize}

  For a description of the algorithms see the book of Schwarz (H.R.Scharz:
  Methode der finiten Elemente).
  */
enum RenumberingMethod {
      Cuthill_McKee,
      reverse_Cuthill_McKee
};





/**
  Manage the distribution and numbering of the degrees of freedom for
  non-multigrid algorithms.

  We store a list of numbers for each cells
  denoting the mapping between the degrees of freedom on this cell
  and the global number of this degree of freedom; the number of a
  degree of freedom lying on the interface of two cells is thus stored
  twice, but is the same. The numbers refer to the unconstrained
  matrices and vectors. The layout of storage of these indices is
  described in the \Ref{DoFLevel} class documentation.

  Additionally, the DoFHandler is able to generate the condensation
  matrix which connects constrained and unconstrained matrices and
  vectors.

  Finally it offers a starting point for the assemblage of the matrices
  by offering #begin()# and #end()# functions which return iterators
  to walk on the DoF structures as well as the triangulation data.
  These iterators work much like those described in the documentation
  of the #Triangulation# class and of the iterator classes themselved,
  but offer more functionality than pure triangulation iterators. The
  order in which dof iterators are presented by the #++# and #--# operators
  is the same as that for the alike triangulation iterators.
  
  
  \subsection{Distribution of degrees of freedom}

  The degrees of freedom (`dofs') are distributed on the given triangulation
  by the function #distribute_dofs()#. It gets passed a finite element object
  describing how many degrees of freedom are located on vertices, lines, etc.
  It traverses the triangulation cell by cell and numbers the dofs of that
  cell if not yet numbered. For non-multigrid algorithms, only active cells
  are considered.

  Since the triangulation is traversed starting with the cells of the coarsest
  active level and going to more refined levels, the lowest numbers for dofs
  are given to the largest cells as well as their bounding lines and vertices,
  with the dofs of more refined cells getting higher numbers.

  This numbering implies very large bandwiths of the resulting matrices and
  is thus vastly suboptimal for some solution algorithms. For this reason,
  the #DoFHandler# class offers the function #renumber_dofs# which reorders
  the dof numbering according to some scheme. Presently available are the
  Cuthill-McKey (CM) and the Reverse Cuthill-McKey algorithm. These algorithms
  have one major drawback: they require a good starting point, i.e. the degree
  of freedom index afterwards to be numbered zero. This can thus be given by
  the user, e.g. by exploiting knowledge of the actual topology of the
  domain. It is also possible to given several starting indices, which may
  be used to simulate a simple upstream numbering (by giving the inflow
  dofs as starting values) or to make preconditioning faster (by letting
  the dirichlet boundary indices be starting points).

  If no starting index is given, one is chosen by the program, namely one
  with the smallest coordination number (the coordination number is the
  number of other dofs this dof couples with). This dof is usually located
  on the boundary of the domain. There is, however, large ambiguity in this
  when using the hierarchical meshes used in this library, since in most
  cases the computational domain is not approximated by tilting and deforming
  elements and by plugging together variable numbers of elements at vertices,
  but rather by hierarchical refinement. There is therefore a large number
  of dofs with equal coordination numbers. The renumbering algorithms will
  therefore not give optimal results.

  In the book of Schwarz (H.R.Schwarz: Methode der finiten Elemente), it is
  advised to test many starting points, if possible all with the smallest
  coordination number and also those with slightly higher numbers. However,
  this seems only possible for meshes with at most several dozen or a few
  hundred elements found in small engineering problems of the early 1980s
  (the second edition was published in 1984), but certainly not with those
  used in this library, featuring several 10,000 to a few 100,000 elements.

  On the other hand, the need to reduce the bandwidth has decreased since
  with the mentioned number of cells, only iterative solution methods are
  able to solve the resulting matrix systems. These, however, are not so
  demanding with respect to the bandwidth as direct solvers used for
  smaller problems. Things like upstream numbering become much more important
  in recent times, so the suboptimality of the renumbering algorithms is
  not that important any more.

  
  \subsection{Implementation of renumbering schemes}

  The renumbering algorithms need quite a lot of memory, since they have
  to store for each dof with which other dofs it couples. This is done
  using a #dSMatrixStruct# object used to store the sparsity pattern. It
  is not useful for the user to do anything between distributing the dofs
  and renumbering, i.e. the calls to #DoFHandler::distribute_dofs# and
  #DoFHandler::renumber_dofs# should follow each other immediately. If
  you try to create a sparsity pattern or anything else in between, these
  will be invalid afterwards.

  The renumbering may take care of dof-to-dof couplings only induced by
  eliminating constraints. In addition to the memory consumption mentioned
  above, this also takes quite some computational time, but it may be
  switched of upon calling the #renumber_dofs# function. This will then
  give inferior results, since knots in the graph (representing dofs)
  are not found to be neighbors even if they would be after condensation.
  
  The renumbering algorithms work on a purely algebraic basis, due to the
  isomorphism between the graph theoretical groundwork underlying the
  algorithms and binary matrices (matrices of which the entries are binary
  values) represented by the sparsity patterns. In special, the algorithms
  do not try to exploit topological knowledge (e.g. corner detection) to
  find appropriate starting points. This way, however, they work in
  arbitrary space dimension.

  If you want to give starting points, you may give a list of dof indices
  which will form the first step of the renumbering. The dofs of the list
  will be consecutively numbered starting with zero, i.e. this list is not
  renumbered according to the coordination number of the nodes. Indices not
  in the allowed range are deleted. If no index is allowed, the algorithm
  will search for its own starting point.

  
  \subsection{Results of renumbering}

  The renumbering schemes mentioned above do not lead to optimal results.
  However, after all there is no algorithm that accomplishes this within
  reasonable time. There are situations where the lack of optimality even
  leads to worse results than with the original, crude, levelwise numering
  scheme; one of these examples is a mesh of four cells of which always
  those cells are refined which are neighbors to the center (you may call
  this mesh a `zoom in' mesh). In one such example the bandwidth was
  increased by about 50 per cent.

  In most other cases, the bandwith is reduced significantly. The reduction
  is the better the less structured the grid is. With one grid where the
  cells were refined according to a random driven algorithm, the bandwidth
  was reduced by a factor of six.

  Using the constraint information usually leads to reductions in bandwidth
  of 10 or 20 per cent, but may for some very unstructured grids also lead
  to an increase. You have to weigh the decrease in your case with the time
  spent to use the constraint information, which usually is several times
  longer than the `pure' renumbering algorithm.

  In almost all cases, the renumbering scheme finds a corner to start with.
  Since there is more than one corner in most grids and since even an
  interior degree of freedom may be a better starting point, giving the
  starting point by the user may be a viable way if you have a simple
  scheme to derive a suitable point (e.g. by successively taking the
  third child of the cell top left of the coarsest level, taking its
  third vertex and the dof index thereof, if you want the top left corner
  vertex). If you do not know beforehand what your grid will look like
  (e.g. when using adaptive algorithms), searching a best starting point
  may be difficult, however, and in many cases will not justify the effort.

  \subsection{Data transfer between grids}

  The #DoFHandler# class offers two functions #make_transfer_matrix# which create
  a matrix to transform the data of one grid to another. The functions assumes the
  coarsest mesh of the two grids to be the same. However there are few ways to
  check this (only the number of cells on the coarsest grid is compared). Also,
  the selected finite element type of the two degree of freedom handler objects
  must be the same.

  The algorithm goes recursively from the coarse mesh cells to their children
  until the grids differ at this level. It then tries to prolong or restrict the
  old cell(s) to the new cell(s) and makes up a matrix of these prolongations and
  restrictions. This matrix multiplied with a vector on the old grid yields an
  approximation of the projection of the function on the old grid to the new one.

  Building and using the transfer matrix is usually a quite expensive operation,
  since we have to perform two runs over all cells (one for building the sparsity
  structure, one to build the entries) and because of the memory consumption.
  It may, however, pay if you have many
  equations, since then the entries in the matrix can be considered as block
  entries which are then applied to all function values at a given degree of
  freedom.

  To build the matrix, you have to call first
  #make_transfer_matrix (old_dof_object, sparsity_pattern);#, then create a
  sparse matrix out of this pattern, e.g. by #dSMatrix m(sparsity_pattern);#
  and finally give this to the second run:
  #make_transfer_matrix (old_dof_object, m);#. The spasity pattern created
  by the first run is automatically compressed.

  When creating the #dSMatrixStruct# sparsity pattern, you have to give the
  dimension and the maximum number of entries per row. Obviously the image
  dimension is the number of dofs on the new grid (you can get this using the
  #n_dofs()# function), while the range dimension is the number of dofs on the
  old grid. The maximum number of entries per row is determined by the maximum
  number of levels $d$ which we have to cross upon transferring from one cell to
  another (presently, transfer of one cell is only possible for #d=0,1#, i.e.
  the two cells match or one is refined once more than the other, the
  number of degrees of freedom per per vertex $d_v$, those on lines $d_l$, those
  on quads $d_q$ and the number of subcells a cell is
  refined to, which is $2**dim$. The maximum number of entries per row in one
  dimension is then given by $(2*d_l+d_v)*2+1$ if $d=1$. For example, a one
  dimensional linear element would need two entries per row.
  In two dimensions, the maxmimum number is $(4*d_q+12*d_l+5*d_v)*4+1$ if $d=1$.
  You can get these numbers by drawing little pictures and counting, there is
  no mystique behind this. You can also get the right number by calling the
  #max_transfer_entries (max_level_difference)# function. The actual number
  depends on the finite element selected and may be much less, especially in
  higher dimensions.
  
  If you do not have multiple equations and do not really use the matrix but still
  have to transfer an arbitrary number of vectors to transfer, you can use the
  #transfer()# function, which is able to transfer any number of vectors in only
  one loop over all cells and without the memory consumption of the matrix. The
  matrix seems only useful when trying to transfer whole matrices instead of
  rebuilding them on the new grid.
  
  @author Wolfgang Bangerth, February 1998
  */
template <int dim>
class DoFHandler : public DoFDimensionInfo<dim> {
  public:
				     // insert these definitions for gcc2.8,
				     // since it can't inherit typedefs (I
				     // believe it should, but it can't)
    typedef typename DoFDimensionInfo<dim>::raw_line_iterator raw_line_iterator;
    typedef typename DoFDimensionInfo<dim>::line_iterator line_iterator;
    typedef typename DoFDimensionInfo<dim>::active_line_iterator active_line_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_quad_iterator raw_quad_iterator;
    typedef typename DoFDimensionInfo<dim>::quad_iterator quad_iterator;
    typedef typename DoFDimensionInfo<dim>::active_quad_iterator active_quad_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_cell_iterator raw_cell_iterator;
    typedef typename DoFDimensionInfo<dim>::cell_iterator cell_iterator;
    typedef typename DoFDimensionInfo<dim>::active_cell_iterator active_cell_iterator;

    typedef typename DoFDimensionInfo<dim>::raw_face_iterator raw_face_iterator;
    typedef typename DoFDimensionInfo<dim>::face_iterator face_iterator;
    typedef typename DoFDimensionInfo<dim>::active_face_iterator active_face_iterator;

    
				     /**
				      * Constructor. Take #tria# as the
				      * triangulation to work on.
				      */
    DoFHandler (Triangulation<dim> *tria);

				     /**
				      * Destructor
				      */
    ~DoFHandler ();
    
				     /**
				      * Go through the triangulation and
				      * distribute the degrees of freedoms
				      * needed for the given finite element
				      * according to the given distribution
				      * method.
				      *
				      * A copy of the transferred finite
				      * element is stored.
				      *
				      * This function uses the user flags of the
				      * triangulation object, so make sure you
				      * don't need them after calling this
				      * function, or if so store them.
				      */
    void distribute_dofs (const FiniteElementBase<dim> &);

				     /**
				      * Renumber the degrees of freedom according
				      * to the given scheme, eventually using
				      * constraint information and the given
				      * starting points. The starting points
				      * default to an empty list, the use of
				      * constraint information defaults to
				      * false.
				      *
				      * See the general documentation of this
				      * class for more details.
				      */
    void renumber_dofs (const RenumberingMethod method,
			bool use_constraints               = false,
			const vector<int> &starting_points = vector<int>());
    
				     /**
				      * Make up the constraint matrix which
				      * is used to condensate the global
				      * system matrices and to prolong
				      * the solution vectors from the true
				      * degrees of freedom also to the
				      * constraint nodes.
				      *
				      * Since this method does not make sense in
				      * one dimension, the functions returns
				      * immediately after clearing the
				      * constraint matrix.
				      * For more than one dimension, the matrix
				      * is cleared before usage. The constraint
				      * matrix is closed anyway, no matter of the
				      * dimension.
				      *
				      * To condense a given sparsity pattern,
				      * use #ConstraintMatrix::condense#.
				      */
    void make_constraint_matrix (ConstraintMatrix &) const;


				     /**
				      * Write the sparsity structure of the
				      * full matrix (including constrained
				      * degrees of freedom) into the
				      * matrix structure. The sparsity
				      * pattern is not compressed, since if
				      * you want to call
				      * #ConstraintMatrix::condense(1)#
				      * afterwards, new entries have to be
				      * added. However, if you want to call
				      * #ConstraintMatrix::condense(1)#, you
				      * have to compress the matrix yourself,
				      * using #dSMatrixStruct::compress()#.
				      */
    void make_sparsity_pattern (dSMatrixStruct &) const; 

				     /**
				      * Make up the transfer matrix which
				      * transforms the data vectors from one
				      * triangulation to the present one.
				      * You have to pass references to the old
				      * dof handler object and to a matrix
				      * sparsity object. This function therefore
				      * only makes up the sparsity pattern.
				      *
				      * The matrix given sparsity pattern is
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
			       dSMatrixStruct        &transfer_pattern) const;

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
			       dSMatrix              &transfer_matrix) const;

				     /**
				      * Return the maximum number of
				      * degrees of freedom a degree of freedom
				      * in the given triangulation with the
				      * given finite element may couple with.
				      * This is the maximum number of entries
				      * per line in the system matrix; this
				      * information can therefore be used upon
				      * construction of the #dSMatrixStruct#
				      * object.
				      *
				      * The returned number is not really the
				      * maximum number but an estimate based
				      * on the finite element and the maximum
				      * number of cells meeting at a vertex.
				      * The number holds for the constrained
				      * matrix also.
				      */
    unsigned int max_couplings_between_dofs () const;

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
				      * e.g. for output.
				      *
				      * It is assumed that the number of
				      * elements in #cell_data# equals the
				      * number of active cells. The size of
				      * #dof_data# is adjusted to the right
				      * size.
				      */
    void distribute_cell_to_dof_vector (const dVector &cell_data,
					dVector       &dof_data) const;

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
    raw_line_iterator
    begin_raw_line   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used line
				      *  on level #level#.
				      */
    line_iterator
    begin_line       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  line on level #level#.
				      */
    active_line_iterator
    begin_active_line(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_line_iterator
    end_line () const;

				     /**
				      *  Return an iterator pointing to the
				      *  last line, used or not.
				      */
    raw_line_iterator
    last_raw_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  line of the level #level#, used or not.

				      */
    raw_line_iterator
    last_raw_line (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line.
				      */
    line_iterator
    last_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used line on level #level#.
				      */
    line_iterator
    last_line (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active line.
				      */
    active_line_iterator
    last_active_line () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active line on level #level#.
				      */
    active_line_iterator
    last_active_line (const unsigned int level) const;
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
    raw_quad_iterator
    begin_raw_quad   (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first used quad
				      *  on level #level#.
				      */
    quad_iterator
    begin_quad       (const unsigned int level = 0) const;

				     /**
				      *  Return iterator to the first active
				      *  quad on level #level#.
				      */
    active_quad_iterator
    begin_active_quad(const unsigned int level = 0) const;

				     /**
				      *  Return iterator past the end; this
				      *  iterator serves for comparisons of
				      *  iterators with past-the-end or
				      *  before-the-beginning states.
				      */
    raw_quad_iterator
    end_quad () const;

				     /**
				      *  Return an iterator pointing to the
				      *  last quad, used or not.
				      */
    raw_quad_iterator
    last_raw_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  quad of the level #level#, used or not.

				      */
    raw_quad_iterator
    last_raw_quad (const unsigned int level) const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad.
				      */
    quad_iterator
    last_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  used quad on level #level#.
				      */
    quad_iterator
    last_quad (const unsigned int level) const;

    				     /**
				      *  Return an iterator pointing to the last
				      *  active quad.
				      */
    active_quad_iterator
    last_active_quad () const;

				     /**
				      *  Return an iterator pointing to the last
				      *  active quad on level #level#.
				      */
    active_quad_iterator
    last_active_quad (const unsigned int level) const;
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
				      * Return a constant reference to the
				      * selected finite element object.
				      */
    const FiniteElementBase<dim> & get_selected_fe () const;

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
    
  protected:
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
				      * Distribute dofs on the given cell,
				      * with new dofs starting with index
				      * #next_free_dof#. Return the next
				      * unused index number. The finite
				      * element used is the one given to
				      * #distribute_dofs#, which is copied
				      * to #selected_fe#.
				      *
				      */
    int distribute_dofs_on_cell (active_cell_iterator &cell,
				 unsigned int next_free_dof);

				     /**
				      * Actually do the renumbering prepared
				      * by the #renumber_dofs# function. Since
				      * this is dimension specific, we
				      * need to have another function.
				      */
    void do_renumbering (const vector<int> &);

				     /**
				      * Make up part of the sparsity pattern of
				      * the transfer matrix by looking at the
				      * two cells given.
				      */
    void transfer_cell (const cell_iterator &old_cell,
			const cell_iterator &new_cell,
			dSMatrixStruct      &transfer_pattern) const;

    				     /**
				      * Make up part of the transfer matrix by
				      * looking at the two cells given.
				      */
    void transfer_cell (const cell_iterator &old_cell,
			const cell_iterator &new_cell,
			dSMatrix            &transfer_matrix) const;

				     /**
				      * Address of the triangulation to
				      * work on.
				      */
    Triangulation<dim>       *tria;

				     /**
				      * Space to store the DoF numbers for the
				      * different levels. Analogous to the
				      * #levels[]# tree of the \Ref{Triangulation}
				      * objects.
				      */
    vector<DoFLevel<dim>*>    levels;

				     /**
				      * Store a copy of the finite element given
				      * latest for the distribution of dofs. In
				      * fact, since the FE base class itself has
				      * not much functionality, this object only
				      * stores numbers such as the number of
				      * dofs per line, but also restriction
				      * and prolongation matrices, etc.
				      */
    const FiniteElementBase<dim>   *selected_fe;

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
};






/*----------------------------   dof.h     ---------------------------*/
/* end of #ifndef __dof_H */
#endif
/*----------------------------   dof.h     ---------------------------*/
