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
{};


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
    struct Cell : public CellData<dim>
    {
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
    static void track_back (vector<Cell>                               &cells,
			    stack<unsigned int, vector<unsigned int> > &rotation_states,
			    unsigned int                                track_back_to_cell);

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
