//----------------------------  grid_reordering_internal.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  grid_reordering_internal.h  ---------------------------
#ifndef __deal2__grid_reordering_internal_h
#define __deal2__grid_reordering_internal_h


#include <base/config.h>
#include <grid/tria.h>

#include <map>
#include <vector>





namespace internal
{
/**
 * Implement the algorithm described in the documentation of the
 * GridReordering<2> class.
 *
 * @author Michael Anderson, 2003
 */
  namespace GridReordering2d
  {

                                     /**
                                      * Check whether a given
                                      * arrangement of cells is
                                      * already consisten. If this is
                                      * the case, then we skip the
                                      * reordering pass.
                                      *
                                      * This function works by looping
                                      * over all cells, checking
                                      * whether one of its faces
                                      * already exists in a list of
                                      * edges, and if it already
                                      * exists in reverse order, then
                                      * return @p{false}. If it is not
                                      * already in the list, or in the
                                      * correct direction, then go on
                                      * with the next faces or cell.
                                      */
    bool
    is_consistent (const std::vector<CellData<2> > &cells);
    
    
/**
 * Defines a variety of variables related to the connectivity of a
 * simple quad element. This includes the nodes on each edge, which
 * edges come into each node and what the default deal.II directions
 * are for the quad.
 *
 * @begin{verbatim}
 *       s2
 *       
 *     +-->--+       
 *     |3   2|     
 * s3  ^     ^ s1   
 *     |0   1|     
 *     +-->--+               
 *   
 *       s0           
 * @end{verbatim}
 *
 * @author Michael Anderson, 2003
 */ 
    class ConnectGlobals
    {
      public:
					 /**	
					  * The nodes on each edge in
					  * anti-clockwise order
					  * { {0,1},{1,2},{2,3},{3,0} }
					  */
	static const int EdgeToNode[4][2];

					 /**
					  * The edges comin into each
					  * node, in anti-clockwise
					  * order
					  * { {3,0},{0,1},{1,2},{2,3} }
					  */
	static const int NodeToEdge[4][2];
  
					 /**
					  * The nodes on each edge in
					  * "default direction order".
					  * {{0,1},{1,2},{3,2},{0,3}}
					  */
	static const int DefaultOrientation[4][2];
    };


/**
 * An enriched quad with information about how the mesh fits together
 * so that we can move around the mesh efficiently.
 *
 * @author Michael Anderson, 2003
 */
    class MQuad
    {
      public:
					 /**
					  * v0 - v3 are indexes of the
					  * vertices of the quad, s0 -
					  * s3 are indexes for the
					  * sides of the quad
					  */
	MQuad (const unsigned int  v0,
	       const unsigned int  v1,
	       const unsigned int  v2,
	       const unsigned int  v3,
	       const unsigned int  s0,
	       const unsigned int  s1,
	       const unsigned int  s2,
	       const unsigned int  s3,
	       const CellData<2>  &cd);
	
					 /**
					  * Stores the vertex numbers
					  */
	unsigned int v[4];
					 /**
					  * Stores the side numbers
					  */
	unsigned int side[4]; 

					 /**
					  * Copy of the @p{CellData} object
					  * from which we construct the
					  * data of this object.
					  */
	CellData<2>  original_cell_data;
	
					 /**
					  * Makes an MQuad from the
					  * given CellData and MSide
					  * list.  Is derived from
					  * binary_function to be
					  * usable with STL
					  * containers.
					  *
					  * Also assumes that the
					  * edges listed present in
					  * the CellData are already
					  * present in the elist
					  * vector.
					  */ 
	struct MakeQuad;
    };

/**
 * The enriched side class containing connectivity information.
 * Orientation is from v0 to v1; Initially this should have v0<v1.
 * After global orientation could be either way.
 *
 * @author Michael Anderson, 2003
 */
    struct MSide
    {
					 /**
					  * Constructor.
					  */
	MSide (const unsigned int initv0,
	       const unsigned int initv1);

					 /**
					  * Return whether the sides
					  * are equal, even if their
					  * ends are reversed.
					  */
	bool operator==(const MSide& s2) const;

					 /**
					  * Return the opposite.
					  */
	bool operator!=(const MSide& s2) const;
	
	unsigned int v0;
	unsigned int v1;
	unsigned int Q0;
	unsigned int Q1;

					 /**
					  * Local side numbers on quads 0 and 1.
					  */
	unsigned int lsn0, lsn1;
	bool Oriented;

					 /**
					  * This class makes a MSide have v0<v1
					  */
	struct SideRectify;

					 /**
					  * Provides a side ordering,
					  * s1<s2, without assuming
					  * v0<v1 in either of the
					  * sides.
					  */
	struct SideSortLess;
    };



/**
 * Implement the 2d algorithm for grid reordering described in the
 * documentation of the @ref{GridReordering} class.
 *
 * @author Michael Anderson, 2003
 */   
    class GridReordering
    {
      public:

					 /**
					  * Do the work intended by
					  * this class.
					  */
	void reorient(std::vector<CellData<2> > &quads);
      private:
	
					 /** 
					  * Sets up the internal data
					  * structures so that the we can
					  * do side hopping and face
					  * switching efficiently. This
					  * means we need a whole bunch of
					  * connectivity information
					  */
	void build_graph (const std::vector<CellData<2> > &inquads);

					 /** 
					  * Orient the internal data
					  * into deal.II format The
					  * orientation algorith is as
					  * follows
					  *
					  * 1) Find an unoriented quad
                                          (A)
					  *
					  * 2) Orient an un_oriented
                                          side (s) of (A)
					  *
					  * 3) side hop on (s) of (A)
                                          to get (B)
					  *
					  * 4) if opposite side to (s)
					  * of (B) is unoriented
					  * orient it
					  *
					  * 5) repeat 3) and 4) until
					  * side-hoppong fails (we've
					  * reached a boundary) or (s)
					  * has already been oriented
					  * (we've closed a loop or
					  * unoriented sides).
					  *
					  * 6) Repeat 2), 3) ,4) and
					  * 5) on other unoriented
					  * sides of (A)
					  *
					  * 7) Choose a new unoriented
					  * A.
					  */
	void orient();
  
					 /**
					  * Get the (now correctly
					  * oriented if we've called
					  * orient) quads.
					  */
	void get_quads(std::vector<CellData<2> > &outquads) const;

					 /**
					  * Orient_side(qnum,lsn)
					  * orients the local side lsn
					  * of the quad qnum in the
					  * triangulation. If the side
					  * opposite lsn is oriented
					  * then lsn is oriented to
					  * match it. Otherwise it is
					  * oriented in the "default"
					  * direction for the quad.
					  */
	void orient_side (const unsigned int quadnum,
			  const unsigned int localsidenum);
 
					 /**
					  * Returns true if all sides
					  * of the quad quadnum are
					  * oriented.
					  */
	bool is_fully_oriented_quad (const unsigned int quadnum) const;
  
					 /**
					  * Returns true if the side lsn
					  * of the quad quadnum is
					  * oriented.
					  */
	bool is_oriented_side (const unsigned int quadnum,
			       const unsigned int lsn) const;

					 /**
					  * Returns true is the side is
					  * oriented in the "default"
					  * direction
					  */
	bool is_side_default_oriented (const unsigned int qnum,
				       const unsigned int lsn) const;
  
					 /**
					  * Increases UnOrQLoc from
					  * it's original value to the
					  * next quad with an
					  * unoriented side. Returns
					  * true if there was another
					  * unoriented quad.
					  */
	bool get_unoriented_quad (unsigned int &UnOrQLoc) const;

					 /**
					  * Sets sidenum to the local
					  * sidenumber of an
					  * unoriented side of the
					  * quad quadnum. Returns true
					  * if such a side exists.
					  */
	bool get_unoriented_side (const unsigned int quadnum,
				  unsigned int &sidenum) const;
  
					 /**
					  * side_hop(&qnum, &lsn) has
					  * qnum being the quadnumber
					  * of a quad in the
					  * triangulation, and a local
					  * side number. side_hop then
					  * sets qnum to the
					  * quadnumber across the
					  * other side of the side,
					  * and sets lsn so that
					  * quads[qnum].sides[lsn] is
					  * the same before and after
					  * the call.  if there is no
					  * other quad on the other
					  * side of the current quad,
					  * then side_hop returns
					  * false.
					  */
	bool side_hop (unsigned int &qnum,
		       unsigned int &lsn) const;
  
					 /**
					  * Sets lsn so that it points
					  * to the opposite side of
					  * the current quad (qnum)
					  * that it was originally
					  * pointing to.
					  */
	bool switch_faces (unsigned int &qnum,
			   unsigned int &lsn) const;

					 /**
					  * A list of enriched
					  * sides/edges of the mesh.
					  */
	std::vector<MSide> sides;
					 /**
					  * A list of enriched quads
					  * in the mesh.
					  */
	std::vector<MQuad> mquads;
    };
  }  // namespace GridReordering2d


/**
 * Implement the algorithm described in the documentation of the
 * GridReordering<3> class.
 *
 * @author Michael Anderson, 2003
 */
  namespace GridReordering3d
  {

                                     /**
                                      * During building the
                                      * conectivity information we
                                      * don't need all the heavy duty
                                      * information about edges that
                                      * we will need later. So we can
                                      * save memory and time by using
                                      * a light-weight class for edges.
                                      */
    struct CheapEdge
    {
                                         /**
                                          * The first node
                                          */
        int node0;
                                         /**
                                          * The second node
                                          */
        int node1;
                                         /**
                                          * A simple constructor
                                          */
        CheapEdge(int n0, int n1);
                                         /**
                                          * Need a partial ordering
                                          * for the STL
                                          */
        bool operator< (const CheapEdge & e2) const;
    };

  
    class ElementInfo
    {
      public:
                                         /**
                                          * The numbers of the edges
                                          * coming into node i are
                                          * given by
                                          * edge_to_node[i][k] where
                                          * k=0,1,2.
                                          */
        int edge_to_node[8][3];

                                         /**
                                          * The orientation of edge
                                          * coming into node i is
                                          * given by
                                          * edge_to_node_orient[i][k]
                                          * where k=0,1,2. 1 means the
                                          * given node is the start of
                                          * the edge -1 means the end
                                          * of the edge.
                                          */
        int edge_to_node_orient[8][3];

                                         /**
                                          * nodesonedge[i][0] is the
                                          * start node for edge i.
                                          * nodesonedge[i][1] is the
                                          * end node for edge i.
                                          */
        int nodes_on_edge[12][2];
        int nodes_on_face[6][4];
    };

    class DealElemInfo : public ElementInfo
    {
      public:
        DealElemInfo();
    };



 
                                     /**
                                      * A conectivity and orientation
                                      * aware edge class.
                                      */
    class Edge
    {
      public:
                                         /**
                                          * Simple constructor
                                          */
        Edge (int n0, int n1, int orient=0)
                        :
                        orientation_flag(orient),
                        group(0),
                        num_neighbouring_cubes(0),
                        neighbouring_cubes(NULL)
          {nodes[0]=n0; nodes[1]=n1;};
      
                                         /**
                                          * Simple Destructor
                                          */
        ~Edge();
      
                                         /**
                                          * The IDs for the end nodes
                                          */
        int nodes[2];
                                         /** 
                                          * Whether the edge has been
                                          * oriented (0), points from
                                          * node 0 to node 1 (1), or
                                          * the reverse (-1).
                                          */
        int orientation_flag;

                                         /** 
                                          * Used to determine which
                                          * "sheet" of parallel edges
                                          * the edge falls in when
                                          * oriented. 0 means not yet
                                          * decided.
                                          */
        int group;

        unsigned int num_neighbouring_cubes;
        unsigned int * neighbouring_cubes;
    };

                                     /**
                                      * A connectivity and orientation
                                      * aware cell.
                                      *
                                      * The connectivity of the cell
                                      * is not contained within. (This
                                      * was for flexibility in using
                                      * deal.II's ordering of edges or
                                      * the XDA format etc.) For this
                                      * information we need the
                                      * ElemInfo class.
                                      *
                                      * One thing we do know is that
                                      * the first four edges in the
                                      * edge class are parallel, as
                                      * are the second four, and the
                                      * third four.
                                      *
                                      * TODO: Need to move conectivity information out 
                                      *       of cell and into edge.
                                      */
    class Cell
    {
      public:
                                         /**
                                          * The IDs for each of the edges.
                                          */
        int edges[12];
        
                                         /**
                                          * The IDs for each of the nodes.
                                          */        
        int nodes[8];  

                                         /**
                                          * Which way do the edges
                                          * point.  Whether node 0 of
                                          * the edge is the base of
                                          * the edge in local element
                                          * (1) or node 1 is the base
                                          * (-1).
                                          */
        int local_orientation_flags[12];
        
                                         /**
                                          * An internal flag used to
                                          * determine whether the cell
                                          * is in the queue of cells
                                          * to be oriented in the
                                          * current sheet.
                                          */
        bool waiting_to_be_processed;  

                                         /**
                                          * Copy Constructor
                                          */
        Cell (const Cell &c)
	  {
	    for(int i=0;i<12;++i)
              {
                edges[i]=c.edges[i];
                local_orientation_flags[i]=c.local_orientation_flags[i];
              }
	    for(int i=0;i<8;++i)
              {
                nodes[i]=c.nodes[i];
              }
	    waiting_to_be_processed=c.waiting_to_be_processed;
	  }

                                         /**
                                          * Default Constructor
                                          */
        Cell()
	  {
	    for(int i=0;i<12;++i)
              {
                edges[i]=-1;
                local_orientation_flags[i]=1;
              }
	    for(int i=0;i<8;++i)
              {
                nodes[i]=-1;
              }
	    waiting_to_be_processed=false;
	  }
    };


                                     /**
                                      * This holds all the pieces for
                                      * orientation together.
                                      * 
                                      * Contains lists of nodes, edges
                                      * and cells.  As well as the
                                      * information about how they all
                                      * connect together.
                                      */
    class Mesh
    {
      public:
                                         /**
                                          * Information about how a
                                          * cell is built up from
                                          * nodes and edges.
                                          */
        const ElementInfo & info;
      
                                         /**
                                          * The list of nodes
                                          */
        std::vector< Point<3> > node_list;
                                         /**
                                          * The list of edges 
                                          */
        std::vector<Edge> edge_list;
                                         /**
                                          * The list of cells
                                          */
        std::vector<Cell> cell_list;

                                         /**
                                          * Checks whether every cell
                                          * in the mesh is sensible. By
                                          * calling
                                          * sanity_check(cell_num) on
                                          * every cell.
                                          */
        bool sanity_check() const;
        
                                         /**
                                          * Checks that every node
                                          * matches with its edges. By
                                          * calling
                                          * sanity_check(cell_num,node_num)
                                          * for each node.
                                          */
        bool sanity_check(int cell_num) const;
        
                                         /**
                                          * Checks that each edge
                                          * going into a node is
                                          * correctly setup.
                                          */
        bool sanity_check_node(int cell_num, int i) const;

                                         /**
                                          * Default Constructor
                                          */
        Mesh(const ElementInfo & INFO): info(INFO) { }

                                         /**
                                          * Prints all information
                                          * about the mesh
                                          */
        void dump() const;
        
                                         /**
                                          * Prints all information
                                          * about the cell
                                          */
        void dump_cell(const Cell &c) const;

                                         /**
                                          * Writes edge information to
                                          * a file.
                                          */
        void dump_edges(char const * const fname) const;

      private:
                                         /**
                                          * Unimplemented private copy
                                          * constructor to disable it.
                                          */
        Mesh(const Mesh&);
                                         /**
                                          * Unimplemented private
                                          * assignemnet operator to
                                          * disable it.
                                          */
        Mesh& operator=(const Mesh&);
    };


    class Orienter
    {
      public:
                                         /**
                                          * Constructor.
                                          */
        Orienter()
          {
            for (unsigned int i=0; i<12; ++i)
              edge_orient_array[i]=false;
          }


                                         /**
                                          * The cube we're looking at now.
                                          */
        unsigned int cur_posn;
        
                                         /**
                                          * We have fully oriented all
                                          * cubes before this one.
                                          */
        unsigned int marker_cube;

        std::vector<int> SheetToProcess;

        int cur_edge_group;

        bool edge_orient_array[12];
      
        bool orient_edges (Mesh &m);
        void orient_cubes (Mesh &m);
      
        bool GetNextUnorientedCube (Mesh &m);
        bool is_oriented (const Mesh &m,
                          int cell_num);

        bool OrientEdgesInCurrentCube (Mesh &m);
        bool OrientEdgeSetInCurrentCube (Mesh &m,
                                         int edge_set);
        bool OrientNextUnorientedEdge (Mesh &m);
        bool Consistant (Mesh &m,
                         int cell_num);


        void GetAdjacentCubes (Mesh &m);
        bool GetNextActiveCube (Mesh &m);

        bool CheckCellEdgeGroupConsistancy (const Mesh &m,
                                            const Cell &c) const;
        
        bool CheckCellEdgeGroupConsistancy (const Mesh &m,
                                            const Cell & c,
                                            int egrp) const;

    };


                                     /**
                                      * Creates the connectivity
                                      * information for the mesh m.
                                      */
    void build_mesh (Mesh &m);


    
    
  }  // namespace GridReordering3d
}  // namespace internal



#endif
