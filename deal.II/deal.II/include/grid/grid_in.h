/*----------------------------   grid_in.h     ---------------------------*/
/*      $Id$                 */
#ifndef __grid_in_H
#define __grid_in_H
/*----------------------------   grid_in.h     ---------------------------*/



#include <base/exceptions.h>
#include <base/smartpointer.h>
#include <grid/forward_declarations.h>
#include <iostream>
#include <vector>
#include <string>



/**
 * This class implements an input mechanism for grid data. It allows to
 * read a grid structure into a triangulation object. Future versions
 * will also allow to read data on this grid into vectors.
 *
 * At present, only UCD (unstructured cell data) is supported as input
 * format for grid data. Any numerical data after the block of topological
 * information is ignored.
 *
 * To read grid data, the triangulation to be fed with has to be empty.
 * When giving a file which does not contain the assumed information or
 * which does not keep to the right format, the state of the triangulation
 * will be undefined afterwards. Upon input, only lines in one dimension
 * and line and quads in two dimensions are accepted. All other cell types
 * (e.g. triangles in two dimensions, quads and hexes in 3d) are rejected.
 * The vertex and cell numbering in the UCD file, which
 * need not be consecutively, is lost upon transfer to the triangulation
 * object, since this one needs consecutively numbered elements.
 *
 * Material indicators are accepted to denote the material id of cells and
 * to denote boundary part indication for lines in 2D. Read the according
 * sections in the documentation of the \Ref{Triangulation} class for
 * further details.
 *
 *
 * \subsection{Structure of input grid data}
 * 
 * It is your duty to use a correct numbering of vertices in the cell list,
 * i.e. for lines, you have to first give the vertex with the lower coordinate
 * value, then that with the higher coordinate value. For quadrilaterals in
 * two dimensions, the vertex indices in the #quad# list have to be such that
 * the vertices are numbered in counter-clockwise sense.
 *
 * In two dimensions, another difficulty occurs, which has to do with the sense
 * of a quadrilateral. A quad consists of four lines which have a direction,
 * which is per definitionem as follows:
 * \begin{verbatim}
 *   3-->--2
 *   |     |
 *   ^     ^
 *   |     |
 *   0-->--1
 * \end{verbatim}
 * Now, two adjacent cells must have a vertex numbering such that the direction
 * of the common side is the same. For example, the following two quads
 * \begin{verbatim}
 *   3---4---5
 *   |   |   |
 *   0---1---2
 * \end{verbatim}
 * may be characterised by the vertex numbers (0 1 4 3) and (1 2 5 4), since
 * the middle line would get the direction #1->4# when viewed from both cells.
 * The numbering (0 1 4 3) and (5 4 1 2) would not be allowed, since the left
 * quad would give the common line the direction #1->4#, while the right one
 * would want to use #4->1#, leading to ambiguity. The #Triangulation# object
 * is capable of detecting this special case, which can be eliminated by
 * rotating the indices of the right quad by two. However, it would not
 * know what to do if you gave the vertex indices (4 1 2 5), since then it
 * would have to rotate by one element or three, the decision which to take is
 * not yet implemented.
 *
 * There are more ambiguous cases, where the triangulation may not know what
 * to do at all without the use of very sophisticated algorithms. On such example
 * is the following:
 * \begin{verbatim}
 *   9---10-----11
 *   |   |    / |
 *   6---7---8  |
 *   |   |   |  |
 *   3---4---5  |
 *   |   |    \ |
 *   0---1------2
 * \end{verbatim}
 * Assume that you had numbered the vertices in the cells at the left boundary
 * in a way, that the following line directions are induced:
 * \begin{verbatim}
 *   9->-10-----11
 *   ^   ^    / |
 *   6->-7---8  |
 *   ^   ^   |  |
 *   3->-4---5  |
 *   ^   ^    \ |
 *   0->-1------2
 * \end{verbatim}
 * (This could for example be done by using the indices (0 1 4 3), (3 4 7 6),
 * (6 7 10 9) for the three cells). Now, you will not find a way of giving
 * indices for the right cells, without introducing either ambiguity for
 * one line or other, or without violating that within each cells, there must be
 * one vertex from which both lines are directed away and the opposite one to
 * which both adjacent lines point to.
 *
 * The solution in this case is to renumber one of the three left cells, e.g.
 * by reverting the sense of the line between vertices 7 and 10 by numbering
 * the top left cell by (9 6 7 10).
 *
 * But this is a thing that the triangulation
 * object can't do for you, since it would involve backtracking to cells
 * already created when we find that we can't number the indices of one of
 * the rightmost cells consistently. It is neither clear how to do this
 * backtracking nor whether it can be done with a stopping algorithm, if
 * possible within polynomial time. This kind of numbering must be made
 * upon construction of the coarse grid, unfortunately.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class GridIn
{
  public:
				     /**
				      * Constructor.
				      */
    GridIn ();
    
				     /**
				      * Attach this triangulation
				      * to be fed with the grid data.
				      */
    void attach_triangulation (Triangulation<dim> &tria);

				     /**
				      * Read grid data from an ucd file.
				      * Numerical data is ignored.
				      */
    void read_ucd (istream &);

				     /**
				      * Exception
				      */
    DeclException1 (ExcUnknownIdentifier,
		    string,
		    << "The identifier <" << arg1 << "> as name of a "
		    << "part in an UCD input file is unknown or the "
		    << "respective input routine is not implemented.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoTriangulationSelected);
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidVertexIndex,
		    int, int,
		    << "Trying to access invalid vertex index " << arg2
		    << " while creating cell " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
    
  private:
				     /**
				      * Store address of the triangulation to
				      * be fed with the data read in.
				      */
    SmartPointer<Triangulation<dim> > tria;
};






/*----------------------------   grid_in.h     ---------------------------*/
/* end of #ifndef __grid_in_H */
#endif
/*----------------------------   grid_in.h     ---------------------------*/
