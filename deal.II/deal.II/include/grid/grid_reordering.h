//----------------------------  grid_reordering.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_reordering.h  ---------------------------
#ifndef __deal2__grid_reordering_h
#define __deal2__grid_reordering_h


#include <map>
#include <vector>
#include <grid/tria.h>


/**
 * Class declaring some dimension dependent numbers which are needed
 * for the grid reordering class.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
class GridReorderingInfo
{
};


/**
 * Class declaring some dimension dependent numbers which are needed
 * for the grid reordering class. This is the specialization for the
 * 2d case.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <>
class GridReorderingInfo<2>
{
  public:
				     /**
				      * Number of possible valid
				      * orientations of a cell. They
				      * are the state in which it was
				      * delivered and three possible
				      * rotations in counter-clockwise
				      * sense, thus a total of four.
				      */
    static const unsigned int rotational_states_of_cells = 4;

				     /**
				      * Number of possible
				      * orientations of a face in
				      * 2d. It is the face and the
				      * face with vertices exchanged,
				      * thus two.
				      */
    static const unsigned int rotational_states_of_faces = 2;
};



/**
 * Class declaring some dimension dependent numbers which are needed
 * for the grid reordering class. This is the specialization for the
 * 2d case.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <>
class GridReorderingInfo<3>
{
  public:
				     /**
				      * ???
				      */
    static const unsigned int rotational_states_of_cells = static_cast<unsigned int>(-1);
    static const unsigned int rotational_states_of_faces = 4;
};


  


/**
 * This class reorders the vertices of cells such that they meet the
 * requirements of the @ref{Triangulation} class when creating
 * grids. This class is mainly used when reading in grids from files
 * and converting them to deal.II triangulations.
 *
 *
 * @sect3{Statement of problems}
 *
 * Triangulations in deal.II have a special structure, in that there
 * are not only cells, but also faces, and in 3d also edges, that are
 * objects of their own right. Faces and edges have unique
 * orientations, and they have a specified orientation also with
 * respect to the cells that are adjacent. Thus, a line that separates
 * two cells in two space dimensions does not only have a direction,
 * but it must also have a well-defined orientation with respect to
 * the other lines bounding the two quadrilaterals adjacent to the
 * first line. Likewise definitions hold for three dimensional cells
 * and the objects (lines, quads) that separate them.
 *
 * For example, in two dimensions, a quad consists of four lines which
 * have a direction, which is by definition as follows:
 * @begin{verbatim}
 *   3-->--2
 *   |     |
 *   ^     ^
 *   |     |
 *   0-->--1
 * @end{verbatim}
 * Now, two adjacent cells must have a vertex numbering such that the direction
 * of the common side is the same. For example, the following two quads
 * @begin{verbatim}
 *   3---4---5
 *   |   |   |
 *   0---1---2
 * @end{verbatim}
 * may be characterised by the vertex numbers @p{(0 1 4 3)} and
 * @p{(1 2 5 4)}, since the middle line would get the direction @p{1->4}
 * when viewed from both cells.  The numbering @p{(0 1 4 3)} and
 * @p{(5 4 1 2)} would not be allowed, since the left quad would give the
 * common line the direction @p{1->4}, while the right one would want
 * to use @p{4->1}, leading to an ambiguity.
 *
 * As a sidenote, we remark that if one adopts the idea that having
 * directions of faces is useful, then the orientation of the four
 * faces of a cell as shown above is almost necessary. In particular,
 * it is not possible to orient them such that they represent a
 * (counter-)clockwise sense, since then we couldn't already find a
 * valid orientation of the following patch of three cells:
 * @begin{verbatim}
 *       o
 *     /   \
 *   o       o 
 *   | \   / |
 *   |   o   |    
 *   |   |   |
 *   o---o---o
 * @end{verbatim}
 * (The reader is aked to try to find a conforming choice of line
 * directions; it will soon be obvious that there can't exists such a
 * thing, even if we allow that there might be cells with clockwise
 * and counterclockwise orientation of the lines at the same time.)
 * 
 * One might argue that the definition of unique directions for faces
 * and edges, and the definition of directions relative to the cells
 * they bound, is a misfeature of deal.II. In fact, it makes reading
 * in grids created by mesh generators rather difficult, as they
 * usually don't follow these conventions when generating their
 * output. On the other hand, there are good reasons to introduce such
 * conventions, as they can make programming much simpler in many
 * cases, leading to an increase in speed of some computations as one
 * can avoid expensive checks in many places because the orientation
 * of faces is known by assumption that it is guaranteed by the
 * triangulation.
 *
 * The purpose of this class is now to find an ordering for a given
 * set of cells such that the generated triangulation satisfies all
 * the requirements stated above. To this end, we will first show some
 * examples why this is a difficult problem, and then develop an
 * algorithm that finds such a reordering. Note that the algorithm
 * operates on a set of @ref{CellData} objects that are used to
 * describe a mesh to the triangulation class. These objects are, for
 * example, generated by the @ref{GridIn} class, when reading in grids
 * from input files.
 *
 *
 * @sect3{Examples of problems}
 *
 * As noted, reordering the vertex lists of cells such that the
 * resulting grid is not a trivial problem. In particular, it is often
 * not sufficient to only look at the neighborhood of a cell that
 * cannot be added to a set of other cells without violating the
 * requirements stated above. We will show two examples where this is
 * obvious.
 *
 * The first such example is the following, which we will call the
 * ``four cells at the end'' because of the four cells that close of
 * the right end of a row of three vertical cells each (in the
 * following picture we only show one such column of three cells at
 * the left, but we will indicate what happens if we prolong this
 * list):
 * @begin{verbatim}
 *   9---10-----11
 *   |   |    / |
 *   6---7---8  |
 *   |   |   |  |
 *   3---4---5  |
 *   |   |    \ |
 *   0---1------2
 * @end{verbatim}
 * Assume that you had numbered the vertices in the cells at the left boundary
 * in a way, that the following line directions are induced:
 * @begin{verbatim}
 *   9->-10-----11
 *   ^   ^    / |
 *   6->-7---8  |
 *   ^   ^   |  |
 *   3->-4---5  |
 *   ^   ^    \ |
 *   0->-1------2
 * @end{verbatim}
 * (This could for example be done by using the indices @p{(0 1 4 3)}, @p{(3 4 7 6)},
 * @p{(6 7 10 9)} for the three cells). Now, you will not find a way of giving
 * indices for the right cells, without introducing either ambiguity for
 * one line or other, or without violating that within each cells, there must be
 * one vertex from which both lines are directed away and the opposite one to
 * which both adjacent lines point to.
 *
 * The solution in this case is to renumber one of the three left cells, e.g.
 * by reverting the sense of the line between vertices 7 and 10 by numbering
 * the top left cell by @p{(9 6 7 10)}:
 * @begin{verbatim}
 *   9->-10-----11
 *   v   v    / |
 *   6->-7---8  |
 *   ^   ^   |  |
 *   3->-4---5  |
 *   ^   ^    \ |
 *   0->-1------2
 * @end{verbatim}
 *
 * The point here is the following: assume we wanted to prolong the grid to 
 * the left like this:
 * @begin{verbatim}
 *   o---o---o---o---o------o
 *   |   |   |   |   |    / |
 *   o---o---o---o---o---o  |
 *   |   |   |   |   |   |  |
 *   o---o---o---o---o---o  |
 *   |   |   |   |   |    \ |
 *   o---o---o---o---o------o
 * @end{verbatim}
 * Then we run into the same problem as above if we order the cells at
 * the left uniformly, thus forcing us to revert the ordering of one
 * cell (the one which we could order as @p{(9 6 7 10)}
 * above). However, since opposite lines have to have the same
 * direction, this in turn would force us to rotate the cell left of
 * it, and then the one left to that, and so on until we reach the
 * left end of the grid. This is therefore an example we we have to
 * track back right until the first column of three cells to find a
 * consistent ordering, if we had initially ordered them uniformly.
 *
 * As a second example, consider the following simple grid, where the
 * order in which the cells are numbered is important:
 * @begin{verbatim}
 *   3-----2-----o-----o ... o-----7-----6
 *   |     |     |     |     |     |     |
 *   |  0  |  N  | N-1 | ... |  2  |  1  |
 *   |     |     |     |     |     |     |
 *   0-----1-----o-----o ... o-----4-----5
 * @end{verbatim}
 * We have here only indicated the numbers of the vertices that are
 * relevant. Assume that the user had given the cells 0 and 1 by the
 * vertex indices @p{0 1 2 3} and @p{6 7 4 5}. Then, if we follow this
 * orientation, the grid after creating the lines for these two cells
 * would look like this:
 * @begin{verbatim}
 *   3-->--2-----o-----o ... o-----7--<--6
 *   |     |     |     |     |     |     |
 *   ^  0  ^  N  | N-1 | ... |  2  v  1  v
 *   |     |     |     |     |     |     |
 *   0-->--1-----o-----o ... o-----4--<--5
 * @end{verbatim}
 * Now, since opposite lines must point in the same direction, we can
 * only add the cells 2 through N-1 to cells 1 such that all vertical
 * lines point down. Then, however, we cannot add cell N in any
 * direction, as it would have two opposite lines that do not point in
 * the same direction. We would have to rotate either cell 0 or 1 in
 * order to be able to add all the other cells such that the
 * requirements of deal.II triangulations are met.
 *
 * These two examples demonstrate that if we have added a certain
 * number of cells in some oeirntation of faces and can't add the next
 * one without introducingfaces that had already been added in another
 * direction, then it might not be sufficient to only rotate cells in
 * the neighborhood of the the cell that we failed to add. It might be
 * necessary to go back a long way and rotate cells that have been
 * entered long ago.
 *
 *
 * @sect3{Solution}
 *
 * From the examples above, it is obvious that if we encounter a cell
 * that cannot be added to the cells which have already been entered,
 * we can not usually point to a cell that is the culprit and that
 * must be entered in a different oreintation. Furthermore, even if we
 * knew which cell, there might be large number of cells that would
 * then cease to fit into the grid and which we would have to find a
 * different orientation as well (in the second example above, if we
 * rotated cell 1, then we would have to rotate the cells 1 through
 * N-1 as well).
 *
 * The only solution to this problem seems to be the following: if
 * cell N can't be added, the try to rotate cell N-1. If we can't
 * rotate cell N-1 any more, then try to rotate cell N-2 and try to
 * add cell N with all orientations of cell N-1. And so
 * on. Algorithmically, we can visualize this by a tree structure,
 * where node N has as many children as there are possible
 * orientations of node N+1 (in two space dimensions, there are four
 * orientations in which each cell can be constructed from its four
 * vertices; for example, if the vertex indicaes are @p{(0 1 2 3)},
 * then the four possibilities would be @p{(0 1 2 3)}, @p{(1 2 3 0)},
 * @p{(2 3 0 1)}, and @p{(3 0 1 2)}). When adding one cell after the
 * other, we traverse this tree in a depth-first (pre-order)
 * fashion. When we encounter that one path from the root (cell 0) to
 * a leaf (the last cell) is not allowed (i.e. that the orientations
 * of the cells which are encoded in the path through the tree does
 * not lead to a valid triangulation), we have to track back and try
 * another path through the tree.
 *
 * In practice, of course, we do not follow each path to a final node
 * and then find out whether a path leads to a valid triangulation,
 * but rather use an inductive argument: if for all previously added
 * cells the triangulation is a valid one, then we can find out
 * whether a path through the tree can yield a valid triangulation by
 * checking whether entering the present cell would introduce any
 * faces that have a nonunique direction; if that is so, then we can
 * stop following all paths below this point and track back
 * immediately.
 *
 * Nevertheless, it is already obvious that the tree has @p{4**N}
 * leaves in two space dimensions, since each of the N cells can be
 * added in four orientations. Most of these nodes can be discarded
 * rapidly, since firstly the orientation of the first cell is
 * irrelevant, and secondly if we add one cell that has a neighbor
 * that has already been added, then there are already only two
 * possible orientations left, so the total number of checks we have
 * to make until we find a valid way is significantly smaller than
 * @p{4**N}. However, an algorithm is still exponential in time and
 * linear in memory (we only have to store the information for the
 * present path in form of a stack of orientations of cells that have
 * already been added).
 *
 * In fact, the two examples above show that the exponential estimate
 * is not a pessimized one: we indeed have to track back to one of the
 * very first cells there to find a way to add all cells in a
 * consistent fashion.
 *
 * This discouraging situation is geatly improved by the fact that we
 * can find an algorithm that in practice is usually only roughly
 * linear in time and memory. We will describe this algorithm in the
 * following.
 *
 * The first observation is that although there are counterexamples,
 * problems are usually local. For example, in the second example, if
 * we had numbered the cells in a way that neighboring cells have
 * similar cell numbers, then the amount pf backtracking needed is
 * greatly reduced. Therefore, in the implementation of the algorithm,
 * the first step is to renumber the cells in a Cuthill-McKee fashion:
 * start with the cell with the least number of neighbors and assign
 * to it the cell number zero. Then find all neighbors of this cell
 * and assign to them consecutive further numbers. Then find their
 * neighbors that have not yet been numbered and assign to them
 * numbers, and so on. Graphically, this represents finding zones of
 * cells consecutively further away from the initial cells and number
 * them in this front-marching way. This already greatly improves
 * locality of problems and consequently reduced the necessary amount
 * of backtracking.
 *
 * The second point is that we can use some methods to prune the tree,
 * which usually lead to a valid orientation of all cells very
 * quickly.
 *
 * The first such method is based on the observation that if we
 * fail to insert one cell with number N, then this may not be due to
 * cell N-1 unless N-1 is a direct neighbor of N. The reason is
 * abvious: the chosen orientation of cell M could only affect the
 * possibilities to add cell N if either it were a direct neighbor or
 * if there were a sequence of cells that were added after M and that
 * connected cells M and N. Clearly, for M=N-1, the latter cannot be
 * the case. Conversely, if we fail to add cell N, then it is not
 * necessary to track back to cell N-1, but we can track back to the
 * neighbor of N with the largest cell index and which has already
 * been added.
 *
 * Unfortunately, this method can fail to yield a valid path through
 * tree if not applied with care. Consider the following situation,
 * initially extracted from a mesh of 950 cells generated
 * automatically by the program BAMG (this program usually generates
 * meshes that are quite badly balanced, often have many -- somtimes
 * 10 or more -- neighbors of one vertex, and exposed several problems
 * in the initial algorithm):
 * @begin{verbatim}
 * 12----13----14----15
 * |     |     |     |
 * |  5  |  7  |  6  |
 * |     |     |     |
 * 8-----9-----10----11
 * |     |     |     |
 * |  3  |     |  4  |
 * |     |     |     |
 * 4-----5-----6-----7
 * |     |     |     |
 * |  0  |  1  |  2  |
 * |     |     |     |
 * 0-----1-----2-----3
 * @end{verbatim}
 */
/*
 * Note that there is a hole in the middle. Assume now that the user
 * described the first cell 0 by the vertex numbers @p{0 1 5 4}, cell
 * 5 by @p{12 8 9 13}, and cell 6 by @p{10 11 15 14}. All other cells
 * are numbered in the usual way, i.e. starting at the bottom left and
 * counting counterclockwise. Cell 5 therefore is the only one that
 * does not follow this order; however, note that the bottom line of
 * cell 5 given by this order of cell 5 does match with the top line
 * of cell 4 in that orientation. Given this description of cells, the
 * algorithm will start with cell zero and add one cell after the
 * other, up until the sixth one. Then the situation will be the
 * following:
 * @begin{verbatim}
 * 12->--13----14->--15
 * |     |     |     |
 * v  5  v  7  ^  6  ^
 * |     |     |     |
 * 8-->--9-----10->--11
 * |     |     |     |
 * ^  3  ^     ^  4  ^
 * |     |     |     |
 * 4-->--5-->--6-->--7
 * |     |     |     |
 * ^  0  ^  1  ^  2  ^
 * |     |     |     |
 * 0-->--1-->--2-->--3
 * @end{verbatim}
 *
 * Coming now to cell 7, we see that the two
 * opposite lines to its left and right have different directions; we
 * will therefore find no orientation of cell 7 in which it can be
 * added without violation of the consistency of the
 * triangulation. According to the rule stated above, we track back to
 * the neighbor with greatest index, which is cell 6......
 *
 * The second method to prune the tree is that usually we cannot add a
 * new cell since the orientation of one of its neighbors that have
 * already been added is wrong. Thus, if we may try to rotate one of
 * the neighbors (of course making sure that rotating that neighbor
 * does not violate the consistency of the triangulation) in order to
 * allow the present cell to be added.
 *
 * While the first method could be explained in terms of backtracking
 * in the tree of orientations more than one step at once, turning a
 * neighbor means jumping to a totally different place in the
 * tree. For both methods, one can find arguments that they will never
 * miss a path that is valid and only skip paths that are invalid
 * anyway.
 *
 * These two methods have proven extremely efficient. We have been
 * able to read very large grids (several ten thousands of cells)
 * without the need to backtrack much. In particular, the time to find
 * an ordering of the cells was found to be mostly linear in the
 * number of cells, and the time to reorder them is usually much
 * smaller (for example by one order of magnitude) than the time
 * needed to read the data from a file, and also to actually generate
 * the triangulation from this data using the
 * @ref{Triangulation}@p{<dim>::create_triangulation} function.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int dim>
class GridReordering : private GridReorderingInfo<dim>
{
  public:
				     /**
				      * This is the main function,
				      * doing what is announced in the
				      * general documentation of this
				      * class.
				      */
    static void reorder_cells (vector<CellData<dim> > &original_cells);
    
  private:

				     /**
				      * Forward declarations of local
				      * classes.
				      */
    class Cell;
    class Face;
    class FaceData;

				     /**
				      * Typedef for a stack type that
				      * describes the rotational
				      * states of all cells that have
				      * already been fitted into the
				      * grid. It is mostly used like a
				      * stack, but sometimes we need
				      * random access into values
				      * below the top, so we can't use
				      * the @p{stack} adaptor from
				      * STL.
				      */
    typedef vector<unsigned int> RotationStack;
    
				     /**
				      * Class that describes the
				      * properties of cells beyond
				      * what is provided by the data
				      * that is available from the
				      * calling functions of this
				      * class. In particular, several
				      * fields are available that
				      * describe connections of cells
				      * to faces and to
				      * neighbors. These fields are
				      * filled in a first pass before
				      * the actual reoordering starts,
				      * as they are needed for the
				      * latter purpose.
				      *
				      * Since this class is derived
				      * from the @ref{CellData} class,
				      * it also contains all the
				      * information available
				      * beforehand.
				      *
				      * @author Wolfgang Bangerth, 2000
				      */
    class Cell : public CellData<dim>
    {
      public:
					 /**
					  * Value to be used if a
					  * neighbor does not exist,
					  * i.e. if the cell is at the
					  * boundary of the domain
					  * with a certain face.
					  */
	static const unsigned int invalid_neighbor = static_cast<unsigned int>(-1);

					 /**
					  * Pointers to the faces of
					  * this cell and their
					  * rotations. If the first
					  * index is zero, then the
					  * faces denote the faces in
					  * their standard direction
					  * with respect to the
					  * ordering in this cell. If
					  * it is nonzero, then they
					  * denote the faces that
					  * would be needed if the
					  * cell were rotate so often.
					  *
					  * Note that the order in
					  * which the faces for a
					  * specific rotational state
					  * appear is not specified,
					  * as this is not
					  * important. It is only
					  * important that each face
					  * in one rotational state or
					  * other appears once for
					  * each orientation of the
					  * cell.
					  */
	typename map<Face,FaceData>::iterator
	faces[GridReorderingInfo<dim>::rotational_states_of_cells][GeometryInfo<dim>::faces_per_cell];

					 /**
					  * Cell indices of the
					  * neighbors of this cell in
					  * the global array of cells.
					  */
	unsigned int neighbors[GeometryInfo<dim>::faces_per_cell];

					 /**
					  * The index of this cell in
					  * the global array of cells.
					  */
	unsigned int cell_no;

					 /**
					  * If we fail to insert this
					  * cell, then we have to
					  * track back. We could track
					  * back right to the previous
					  * cell, but we can do better
					  * than that, by tracking
					  * back to the cell indicated
					  * by this field. Which value
					  * it has is described in the
					  * documentation of the
					  * @ref{GridReordering}
					  * class.
					  */
	unsigned int track_back_to_cell;

					 /**
					  * Default
					  * constructor. Invalidate
					  * all data.
					  */
	Cell ();
	
					 /**
					  * Constructor that copies
					  * the data of an object of
					  * the base class and
					  * requires to be given the
					  * index of this cell in the
					  * global array.
					  */
	Cell (const CellData<dim> &cd,
	      const unsigned int   cell_no);
	
					 /**
					  * Count the existing neighbors
					  * of this cell.
					  */
	unsigned int count_neighbors () const;

					 /**
					  * Insert the faces of the
					  * present cell into the map
					  * of all faces. This
					  * function inserts them in
					  * all orientations possible
					  * if the given cell is
					  * rotated. The function also
					  * takes care to fill in the
					  * @p{adjacent_cells} field
					  * of the inserted faces.
					  */
	void insert_faces (map<Face,FaceData > &global_faces);

					 /**
					  * Find out the neighbors of the
					  * given cell by looking at the
					  * @p{adjacent_cells} field of
					  * the faces of this cell. Store
					  * the neighbor indices in the
					  * present object.
					  */
	void fix_cell_neighbors ();

					 /**
					  * Compute back to which cell we
					  * have to backtrack in case we
					  * can't insert this cell in any
					  * orientation into the already
					  * existing part of the
					  * triangulation. The method of
					  * how to determine the point to
					  * which we have to backtrack is
					  * described in the documentation
					  * of the @ref{GridReordering}
					  * class.
					  */
	void find_backtracking_point ();

					 /**
					  * Find out whether the cell
					  * could be inserted into the
					  * already existing part of
					  * the triangulation with
					  * orientation given by the
					  * parameter, by checking
					  * that no face would be
					  * inserted twice in
					  * different orientations.
					  */
	bool check_consistency (const unsigned int rot) const;

					 /**
					  * Tag the faces of this cell in
					  * the given orientation as used
					  * by this cell.
					  */
	void mark_faces_used (const unsigned int rot);

					 /**
					  * Remove the use tags on the
					  * faces of this cell in the
					  * given orientation by this
					  * cell. Tags may remain if
					  * there is another cell that
					  * uses a given face.
					  */
	void mark_faces_unused (const unsigned int rot);
    };


				     /**
				      * Structure describing a face of
				      * a cell. This class is used as
				      * key in a map storing all faces
				      * possible in a triangulation,
				      * i.e. all faces between cells
				      * in all possible orientations.
				      *
				      * @author Wolfgang Bangerth, 2000
				      */
    struct Face
    {
					 /**
					  * Indices of the vertices of
					  * this face.
					  */
	unsigned int vertices[GeometryInfo<dim>::vertices_per_face];

					 /**
					  * Comparison operator. Use
					  * the vertex indices as
					  * primary, secondary,
					  * ... criteria for
					  * comparison.
					  */
	bool operator < (const Face &face) const;
    };


				     /**
				      * Class describing some data to
				      * be stored on each
				      * face. Objects of this type are
				      * used as values in a map
				      * containing all faces.
				      *
				      * @author Wolfgang Bangerth, 2000
				      */
    struct FaceData
    {
					 /**
					  * Value denoting
					  * non-existing adjacent
					  * cells of this face,
					  * i.e. when the face is at
					  * the boundary of the domain
					  * and has only one adjacent
					  * cell.
					  */
	static const unsigned int invalid_adjacent_cell = static_cast<unsigned int>(-1);

					 /**
					  * Pointers to the same face
					  * but in all other
					  * orientations. Storing
					  * these pointers makes it
					  * much easier to find out
					  * whether a given faces has
					  * already been used in
					  * another direction, thus
					  * forbidding the present
					  * face to be used.
					  *
					  * Note that the order in
					  * which the reverted faces
					  * appear here is not
					  * specified, as it is not
					  * important for the
					  * algorithm.
					  */
	typename map<Face,FaceData >::const_iterator
	reverse_faces[GridReorderingInfo<dim>::rotational_states_of_faces-1];

					 /**
					  * Indices of the one or two
					  * adjacent cells of this
					  * face in the global array
					  * of cells.
					  */
	unsigned int adjacent_cells[2];

					 /**
					  * Number of cells presently
					  * using this face in the
					  * orientation represented by
					  * this object. May be zero,
					  * one, or two.
					  */
	unsigned int use_count;

					 /**
					  * Default constructor.
					  */
	FaceData ();
    };


				     /**
				      * If we couldn't insert a cell
				      * into the already existing part
				      * of the mesh, then we need to
				      * track back a while. This
				      * function does so, given the
				      * array of cells, the stack of
				      * rotation states, and the
				      * position to which backtracking
				      * shall take place.
				      *
				      * In some cases, it is possible
				      * that from the place where we
				      * backtracked to, there is no
				      * more possibility to orient a
				      * cell. Then we will have to
				      * backtrack further until we
				      * come to a place where further
				      * work is possible; this
				      * recursive backtracking is also
				      * done by this function,
				      * although it is not implemented
				      * as recursive calls but rather
				      * as eliminated tail-recursion.
				      */
    static void track_back (vector<Cell>  &cells,
			    RotationStack &rotation_states,
			    unsigned int   track_back_to_cell);

    static bool try_rotate_single_neighbors (vector<Cell>  &cells,
					     RotationStack &rotation_states);
    
				     /**
				      * This is the main function that
				      * does the main work. It is
				      * called by the
				      * @p{reorder_cells} function
				      * after all the preparations
				      * have been completed and
				      * operates on the @p{cells}
				      * array. After a way to reorder
				      * the cells has been found, the
				      * @p{original_cells} are reorder
				      * accordingly, where the
				      * @p{new_cell_numbers} array is
				      * needed to find the connection
				      * between original cells and
				      * presorted cells.
				      */
    static void find_reordering (vector<Cell>               &cells,
				 vector<CellData<dim> >     &original_cells,
				 const vector<unsigned int> &new_cell_numbers);

				     /**
				      * Preorder the incoming cells by
				      * some kind of Cuthill-McKee
				      * algorithm. The reason for the
				      * need to do so is described in
				      * the general documentation.
				      *
				      * Return a vector in which for
				      * each old cell the new index is
				      * stored.
				      */
    static
    vector<unsigned int>
    presort_cells (vector<Cell>       &cells,
		   map<Face,FaceData> &faces);
};



#endif
