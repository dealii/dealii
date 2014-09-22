// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


/**
 * @defgroup grid Grid classes
 *
 * This module groups classes that have to do with the topology and
 * geometry of meshes. A mesh can be thought of a collection of cells;
 * if the mesh has been refined (possibly in an adaptive way), then
 * this collection is grouped into a hierarchy of refinement
 * levels. In addition to cells, the geometric objects that make up a
 * triangulation are the faces of cells (and in 3d the edges of cells)
 * as well as the vertices of the cells. Note that we abuse the word
 * <i>triangulation</i> somewhat, since deal.II only implements
 * triangulations made up of linear, quadrilateral, and hexahedral
 * cells; triangles and tetrahedra are not supported.
 *
 * This collection of cells is managed by the Triangulation class. It holds
 * the relevant data in memory and offers interfaces to query it. Most things
 * you want to do on cells are performed in loops over all cells. For this
 * purpose, the Triangulation class offers the concept of iterators (see @ref
 * Iterators): although implemented differently, they behave like pointers to
 * cells or faces and can be queried for the geometric properties of cells as
 * well as information like neighboring cells or faces of a cell.
 *
 * It is worth noting that the Triangulation class only stores geometry
 * (i.e. the location of vertices and cells) and topology of a mesh
 * (i.e. which cells are neighbors of which other cells, etc). It has nothing
 * to do with finite elements or degrees of freedom that might be defined on a
 * mesh. These functions are performed by the DoFHandler class (see the @ref
 * dofs module) that gets a description of the finite element space and the
 * allocates and manages degrees of freedom on vertices, faces, or cells, as
 * described by the finite element class. This separation makes it possible to
 * have multiple DoFHandler classes work on the same mesh at the same time.
 * 
 * 
 * <h3>Grid generation</h3>
 *
 * There are three ways to create a mesh:
 * <ul>
 * <li> Creation by the GridGenerator class;
 * <li> Reading from a file;
 * <li> Creation by hand.
 * </ul>
 *
 * For the first case, the GridGenerator class provides functions that can
 * generate the simplest and most common geometries automatically. For
 * example, a rectangular (or brick) geometry as well as circles, spheres, or
 * cylinders can be generate with the functions in this class. Most of the
 * tutorial programs use this mechanism.
 *
 * Secondly, it is possible to read in meshes from an input file in a number
 * of different formats using the GridIn class. Using this class, it is
 * possible to read meshes with several 10 or 100 thousand cells, although
 * this is not really recommended: the power of adaptive finite element
 * methods only comes to bear if the initial mesh is as coarse as possible and
 * there is room for a number of adaptive refinement steps. If the initial
 * mesh is too fine already, then one runs out of memory or compute time
 * before adaptive mesh refinement is able to do much good. Nevertheless, the
 * GridIn class can be used in cases of complicated geometries or for
 * comparison or interaction with other programs that compute on meshes that
 * are then exchanged through this class The step-5 tutorial program shows how
 * to use the GridIn class.
 *
 * The third way is to create a mesh by hand, by building a data structure
 * that describes the vertices and cells of a triangulation. This is useful in
 * cases of moderate complexity where a mesh can still be built by hand
 * without resorting to a mesh generator, but where the domain is not one of
 * those already supported by the GridIn class. In this method, the data
 * structure so built is handed to the create_triangulation() function of the
 * Triangulation class. The step-14 tutorial program shows how this can be
 * done.
 *
 *
 * <h3>Grid output</h3>
 *
 * Meshes can be written to output files in a number of different formats. If
 * this involves simulation results obtained on this mesh, then this is done
 * using the DataOut class (described in more detail in the @ref output
 * module). On the other hand, if only the geometry and topology of the mesh
 * is to be written to a file, the GridOut class can do this for you.
 *
 *
 * <h3>Tool classes</h3>
 *
 * The GridTool class offers an assortment of functions that act on grids. For
 * example, this includes moving around nodes, stretching or rotating entire
 * triangulations, computing the diameter of a domain, or subdividing it into
 * chunks of roughly equal size for parallel computations.
 *
 * The GridRefinement class implements a number of mesh refinement algorithms,
 * based on refinement indicators given to its member functions.
 *
 * 
 * <h3>Internal classes</h3>
 *
 * In addition to the above, there are a significant number of classes in this
 * module that are only used in the internal data structures of mesh
 * handling. They are generally in the internal namespace, and not meant for
 * use in application code.
 *
 * 
 * @author Wolfgang Bangerth, 1998-2006
 */
