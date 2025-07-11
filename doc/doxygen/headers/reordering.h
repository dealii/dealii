// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/**
 * @defgroup reordering Grid reordering and cell orientation
 *
 * @brief A group describing how deal.II consistently orients Triangulation
 * objects.
 *
 * @warning The implementation of orientation should be considered an internal
 * detail of the library. Normal users should not need to use the features
 * described in this group: instead, classes like QProjector use orientation
 * information to consistently compute values on faces and lines.
 *
 * The orientation of a line, triangle, or quadrilateral is a permutation of its
 * vertices which does not result in a twisted cell: i.e., orientations are a
 * subset of all possible vertex permutations. Orientations are encoded with
 * types::geometric_orientation.
 *
 * A Triangulation contains both cells and also lower-dimensional objects, such
 * as faces and vertices. While CellAccessor and TriaAccessor objects provide
 * access to all the topological and geometric information of some entity (such
 * as its neighbors, manifold ids, etc), these objects only store indices and a
 * pointer to a Triangulation and look up the data they use from arrays managed
 * by that Triangulation (typically in an
 * internal::TriangulationImplementation::TriaLevel). This is an example of the
 * <a href="https://en.wikipedia.org/wiki/Flyweight_pattern">flyweight
 * pattern</a>. For brevity, we typically write "the line stores" or "the cell
 * stores" but the underlying implementation contains at least one level of
 * indirection, as this data is stored in some way by the Triangulation and is
 * indexed by TriaAccessorBase::index() and TriaAccessorBase::level().
 *
 * In general, each geometric entity only stores the indices of the lower-level
 * entities which bound it (e.g., a cell stores the indices of its faces, but in
 * 3d must query the faces to get the indices of its lines). One exception to
 * this rule is the vertex index cache: for performance reasons cells directly
 * store their vertex indices.
 *
 * Each line and face has a unique index and exists exactly once in a
 * Triangulation. Orientations are defined as the permutation which makes the
 * vertices, defined in the context of the line or face, have the same order as
 * the ones defined in the context of the cell.
 *
 * <h2>Orientation of Lines</h2>
 *
 * In 1D, lines are cells and, since each line appears in the Triangulation
 * exactly once, they do not store any orientation information.
 *
 * A line in 2D is a face which may be shared by two cells, whereas in 3D a line
 * may be shared by an arbitrary number of cells.
 *
 * In 2D, each cell stores the indices of its bounding lines: i.e., each cell
 * stores either 3 (for triangles) or 4 (for quadrilaterals) integers which
 * enumerate those lines. Put another way: while there is no distinct `Line`
 * class in deal.II, each line is uniquely identified by an index (accessed via
 * CellAccessor::line_index()) and lines are, like other geometric entities,
 * implemented with the flyweight pattern and represented by a TriaAccessor.
 * Each line stores and is defined by the vertex indices which bound it and also
 * stores other auxiliary information (such as boundary and manifold ids).
 *
 * In 3D, the line indices of a cell are stored by the faces which bound it:
 * i.e., each cell stores its own face indices (which represent triangles or
 * quadrilaterals) and each face stores its own line indices. Hence, when
 * accessing a cell's lines, after identifying the face which owns that line
 * (via ReferenceCell::standard_line_to_face_and_line_index()) all data lookups
 * proceed in exactly the same way as the 2D case.
 *
 * Each line has both a unique index and a canonical vertex ordering. For
 * example, consider two triangles with vertices `{15, 20, 25}` and `{20, 15,
 * 3}`. The order of the vertices on each cell is defined by the CellData
 * objects passed to Triangulation::create_triangulation(): typically, the order
 * of the vertices is arbitrary aside from the constraint that they form a cell
 * whose mapping to the reference cell has a positive Jacobian (see
 * Triangulation::Triangulation() for more information on whether or not this
 * should be checked). In this example, the first line of the first cell is
 * `{15, 20}` whereas the first line of the second cell is `{20, 15}`. In
 * deal.II, the canonical order of a line's vertices is set by the first cell
 * with that line: i.e., line `0` will have vertices `{15, 20}` since it first
 * appears in the first cell and those are the first two vertices of that cell.
 * Similarly, line `1` is `{20, 25}` and line `2` is `{25, 15}`. This order
 * (i.e., first and second vertex, then second and third, then third and first)
 * is defined by ReferenceCell::line_to_cell_vertices() for each reference cell
 * type.
 *
 * Canonicalization creates an inconsistency because the vertices of the first
 * line on the second cell are reversed. To resolve this inconsistency, each 2D
 * structure (either cells in 2D or faces in 3D) also stores the relative
 * orientations of its bounding lines. In this example, for the first line, the
 * first cell will store numbers::default_geometric_orientation and the second
 * cell will store numbers::reverse_line_orientation. In each case this
 * orientation value encodes the transformation necessary to make the canonical
 * ordering match the cell-local ordering: i.e., the first cell does nothing and
 * the second cell must invert the order.
 *
 * <h2>Orientation of Faces</h2>
 *
 * Unlike lines, which only have two possible orientations, a quadrilateral
 * (i.e., a face of a pyramid, wedge, or hexahedron in 3D) has 8 possible
 * orientations and a triangle (i.e., a face of a tetrahedron, pyramid, or wedge
 * in 3D) has 6. In deal.II, we encode the orientation of a quadrilateral or
 * triangle with three booleans: orientation, rotation, and flip. The default
 * values for these are true, false, and false. These values are typically
 * encoded or decoded from or to a single types::geometric_orientation value
 * (whose exact binary representation is an internal library detail) by the
 * internal::combined_face_orientation() and internal::split_face_orientation()
 * functions.
 *
 * For a quadrilateral, these values correspond to
 * - *orientation* : `true` is the default orientation and `false` means
 *   vertices 1 and 2 are swapped.
 * - *rotation* : all vertices are rotated by 90 degrees counterclockwise.
 * - *flip* : all vertices are rotated by 180 degrees counterclockwise.
 *
 * For a triangle, these values correspond to
 * - *orientation* : `true` is the default orientation and `false` means
 *   vertices 1 and 2 are swapped.
 * - *rotation* : all vertices are rotated by 120 degrees counterclockwise.
 * - *flip* : all vertices are rotated by 240 degrees counterclockwise.
 *
 * Here, 'clockwise' is relative to the vector defined by the cross product of
 * two lines adjacent to the zeroth vertex in their standard orientation (which,
 * e.g., points into the hexahedron for face 0 but out of the hexahedron for
 * face 1).
 *
 * For triangles, to enable indexing from the combined orientation, we do not
 * consider flip-rotate or flip-orient-rotate as those cases are equivalent,
 * respectively, to the identity operation or the orientation = `true` case as
 * flip-rotate is equal to the identity operation. As a consequence, there are
 * only six valid orientations for triangles as faces of tetrahedra, pyramids,
 * or wedges.
 *
 * Like the line case, the stored orientation value defines the way that the
 * vertices of the face must be permuted to match the cell-local ordering. A
 * consequence of this choice is that QProjector uses the inverse orientation
 * (via ReferenceCell::get_inverse_combined_orientation()) to compute the
 * locations of quadrature points, since exactly one of the following
 * possibilities must happen:
 *
 * 1. If we are on the face which defines the canonical ordering of the face
 *    vertices then that face's orientation must be
 *    numbers::default_geometric_orientation, whose inverse is itself (as it is
 *    the identity permutation). Hence, in this case, applying the permutation
 *    to the positions of the vertices will not change the positions of the
 *    quadrature points.
 *
 * 2. If we are on the neighbor's face then, to make vertices match, we must
 *    transform the cell-local vertices so that they match the first cell's
 *    vertex ordering: i.e., the inverse orientation.
 *
 * <h2>Orientations as Permutations of Vertices</h2>
 * As discussed above, TriaAccessor::combined_face_orientation() returns a value
 * less than ReferenceCell::n_face_orientations() which describes how the
 * bounding object (i.e., a face) should be rotated to match the cell-local
 * definition of that object. The encoding of triangles is
 *
 * <div class="threecolumn" style="width: 600px">
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html triangle-0.svg
 *       @image html triangle-1.svg
 *     </div>
 *   </div>
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html triangle-2.svg
 *       @image html triangle-3.svg
 *     </div>
 *   </div>
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html triangle-4.svg
 *       @image html triangle-5.svg
 *     </div>
 *   </div>
 * </div>
 *
 * Similarly, quadrilaterals are encoded as
 *
 * <div class="fourcolumn" style="width: 800px">
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html quadrilateral-0.svg
 *       @image html quadrilateral-1.svg
 *     </div>
 *   </div>
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html quadrilateral-2.svg
 *       @image html quadrilateral-3.svg
 *     </div>
 *   </div>
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html quadrilateral-4.svg
 *       @image html quadrilateral-5.svg
 *     </div>
 *   </div>
 *   <div class="parent">
 *     <div class="img" align="center">
 *       @image html quadrilateral-6.svg
 *       @image html quadrilateral-7.svg
 *     </div>
 *   </div>
 * </div>
 *
 * The numbers on each edge are the line numbers with the two bounding vertices
 * as subscripts: e.g., $0_0^1$ is the first line of a triangle and when it
 * appears as $0_1^0$ then the line is in the reversed orientation (i.e.,
 * numbers::reverse_line_orientation).
 *
 * <h2>Orientation of Quadrilateral Meshes</h2>
 *
 * Purely quadrilateral meshes are a special case, since deal.II will (with the
 * exception of faces which are neighbors across periodic boundaries)
 * consistently orient purely quadrilateral meshes. Hence, in this case, the
 * orientation of all lines will be numbers::default_geometric_orientation. See
 * @cite AABB17 for more information on this algorithm.
 *
 * For example, in two dimensions, a quad consists of four lines which have a
 * direction, which is by definition as follows:
 * @verbatim
 *   2-->--3
 *   |     |
 *   ^     ^
 *   |     |
 *   0-->--1
 * @endverbatim
 * Now, two adjacent cells must have a vertex numbering such that the
 * direction of the common side is the same. For example, the following two
 * quads
 * @verbatim
 *   3---4---5
 *   |   |   |
 *   0---1---2
 * @endverbatim
 * may be characterized by the vertex numbers {0, 1, 3, 4} and {1, 2, 4, 5},
 * since the middle line would be {1, 4} when viewed from both cells. The
 * numbering {0, 1, 3, 4} and {5, 4, 2, 1} would not be allowed, since the
 * left quad would give the common line the direction {1, 4} whereas
 * the right one would use {4, 1}, leading to an inconsistency.
 *
 * As a sidenote, we remark that if one adopts the idea that having directions
 * of faces is useful, then the orientation of the four faces of a cell as
 * shown above is almost necessary. In particular, it is not possible to
 * orient them such that they represent a (counter-)clockwise sense, since
 * then we couldn't already find a valid orientation of the following patch of
 * three cells:
 * @verbatim
 *       o
 *     /   \
 *   o       o
 *   | \   / |
 *   |   o   |
 *   |   |   |
 *   o---o---o
 * @endverbatim
 * (The reader is asked to try to find a conforming choice of line directions;
 * it will soon be obvious that there can't exists such a thing, even if we
 * allow that there might be cells with clockwise and counterclockwise
 * orientation of the lines at the same time.)
 *
 * One might argue that the definition of unique directions for faces and
 * edges, and the definition of directions relative to the cells they bound,
 * is a misfeature of deal.II. In fact, it makes reading in grids created by
 * mesh generators rather difficult, as they usually don't follow these
 * conventions when generating their output. On the other hand, there are good
 * reasons to introduce such conventions, as they can make programming much
 * simpler in many cases, leading to an increase in speed of some computations
 * as one can avoid expensive checks in many places because the orientation of
 * faces is known by assumption that it is guaranteed by the triangulation.
 *
 * As a last question for this first section: is it guaranteed that such
 * orientations of faces always exist for a given subdivision of a domain into
 * cells? The linear complexity algorithm described below for 2d also proves
 * that the answer is yes for 2d. For 3d, the answer is no (which also
 * underlines that using such orientations might be an -- unfortunately
 * uncurable -- misfeature of deal.II). A simple counter-example in 3d
 * illustrates this: take a string of 3d cells and bend it together to a
 * torus. Since opposing lines in a cell need to have the same direction,
 * there is a simple ordering for them, for example all lines radially
 * outward, tangentially clockwise, and axially upward. However, if before
 * joining the two ends of the string of cells, the string is twisted by 180
 * degrees, then no such orientation is possible any more, as can easily be
 * checked. In effect, some meshes could not be used in deal.II. In order to
 * overcome this problem, the <code>face_rotation</code>,
 * <code>face_flip</code> and <code>line_orientation</code> flags have been
 * introduced. With these, it is possible to treat all purely hexahedral
 * meshes. However, in order to reduce the effect of possible bugs, it should
 * still be tried to reorder a grid. Only if this procedure fails, the
 * original connectivity information should be used.
 *
 * <h3>Examples of problems</h3>
 *
 * As noted, reordering the vertex lists of cells such that the resulting grid
 * is not a trivial problem. In particular, it is often not sufficient to only
 * look at the neighborhood of a cell that cannot be added to a set of other
 * cells without violating the requirements stated above. We will show two
 * examples where this is obvious.
 *
 * The first such example is the following, which we will call the ``four
 * cells at the end'' because of the four cells that close of the right end of
 * a row of three vertical cells each (in the following picture we only show
 * one such column of three cells at the left, but we will indicate what
 * happens if we prolong this list):
 * @verbatim
 *   9---10-----11
 *   |   |    / |
 *   6---7---8  |
 *   |   |   |  |
 *   3---4---5  |
 *   |   |    \ |
 *   0---1------2
 * @endverbatim
 * Assume that you had numbered the vertices in the cells at the left boundary
 * in a way, that the following line directions are induced:
 * @verbatim
 *   9->-10-----11
 *   ^   ^    / |
 *   6->-7---8  |
 *   ^   ^   |  |
 *   3->-4---5  |
 *   ^   ^    \ |
 *   0->-1------2
 * @endverbatim
 * (This could for example be done by using the indices {0, 1, 3, 4},
 * {3, 4, 6, 7}, {6, 7, 9, 10} for the three cells). Now, you will
 * not find a way of giving indices for the right cells, without introducing
 * either ambiguity for one line or other, or without violating that within
 * each cells, there must be one vertex from which both lines are directed
 * away and the opposite one to which both adjacent lines point to.
 *
 * The solution in this case is to renumber one of the three left cells, e.g.
 * by reverting the sense of the line between vertices 7 and 10 by numbering
 * the top left cell by {9, 6, 10, 7}:
 * @verbatim
 *   9->-10-----11
 *   v   v    / |
 *   6->-7---8  |
 *   ^   ^   |  |
 *   3->-4---5  |
 *   ^   ^    \ |
 *   0->-1------2
 * @endverbatim
 *
 * The point here is the following: assume we wanted to prolong the grid to
 * the left like this:
 * @verbatim
 *   o---o---o---o---o------o
 *   |   |   |   |   |    / |
 *   o---o---o---o---o---o  |
 *   |   |   |   |   |   |  |
 *   o---o---o---o---o---o  |
 *   |   |   |   |   |    \ |
 *   o---o---o---o---o------o
 * @endverbatim
 * Then we run into the same problem as above if we order the cells at the
 * left uniformly, thus forcing us to revert the ordering of one cell (the one
 * which we could order as {9, 6, 7, 10} above). However, since opposite
 * lines have to have the same direction, this in turn would force us to
 * rotate the cell left of it, and then the one left to that, and so on until
 * we reach the left end of the grid. This is therefore an example we have
 * to track back right until the first column of three cells to find a
 * consistent ordering, if we had initially ordered them uniformly.
 *
 * As a second example, consider the following simple grid, where the order in
 * which the cells are numbered is important:
 * @verbatim
 *   3-----2-----o-----o ... o-----7-----6
 *   |     |     |     |     |     |     |
 *   |  0  |  N  | N-1 | ... |  2  |  1  |
 *   |     |     |     |     |     |     |
 *   0-----1-----o-----o ... o-----4-----5
 * @endverbatim
 * We have here only indicated the numbers of the vertices that are relevant.
 * Assume that the user had given the cells 0 and 1 by the vertex indices
 * {0, 1, 3, 2} and {6, 7, 5, 4}. Then, if we follow this orientation,
 * the grid after creating the lines for these two cells would look like this:
 * @verbatim
 *   3-->--2-----o-----o ... o-----7--<--6
 *   |     |     |     |     |     |     |
 *   ^  0  ^  N  | N-1 | ... |  2  v  1  v
 *   |     |     |     |     |     |     |
 *   0-->--1-----o-----o ... o-----4--<--5
 * @endverbatim
 * Now, since opposite lines must point in the same direction, we can only add
 * the cells 2 through N-1 to cells 1 such that all vertical lines point down.
 * Then, however, we cannot add cell N in any direction, as it would have two
 * opposite lines that do not point in the same direction. We would have to
 * rotate either cell 0 or 1 in order to be able to add all the other cells
 * such that the requirements of deal.II triangulations are met.
 *
 * These two examples demonstrate that if we have added a certain number of
 * cells in some orientation of faces and can't add the next one without
 * introducing faces that had already been added in another direction, then it
 * might not be sufficient to only rotate cells in the neighborhood of the
 * cell that we failed to add. It might be necessary to go back a long way and
 * rotate cells that have been entered long ago.
 *
 *
 * <h3>Solution</h3>
 *
 * From the examples above, it is obvious that if we encounter a cell that
 * cannot be added to the cells which have already been entered, we can not
 * usually point to a cell that is the culprit and that must be entered in a
 * different orientation. Furthermore, even if we knew which cell, there might
 * be large number of cells that would then cease to fit into the grid and
 * which we would have to find a different orientation as well (in the second
 * example above, if we rotated cell 1, then we would have to rotate the cells
 * 1 through N-1 as well).
 *
 * A brute force approach to this problem is the following: if cell N can't be
 * added, then try to rotate cell N-1. If we can't rotate cell N-1 any more,
 * then try to rotate cell N-2 and try to add cell N with all orientations of
 * cell N-1. And so on. Algorithmically, we can visualize this by a tree
 * structure, where node N has as many children as there are possible
 * orientations of node N+1 (in two space dimensions, there are four
 * orientations in which each cell can be constructed from its four vertices;
 * for example, if the vertex indices are {0 1 3 2}, then the four
 * possibilities would be {0, 1, 3, 2}, {1, 3, 2, 0}, {3, 2, 0, 1},
 * and {2, 0, 1, 3}. When adding one cell after the other, we
 * traverse this tree in a depth-first (pre-order) fashion. When we encounter
 * that one path from the root (cell 0) to a leaf (the last cell) is not
 * allowed (i.e. that the orientations of the cells which are encoded in the
 * path through the tree does not lead to a valid triangulation), we have to
 * track back and try another path through the tree.
 *
 * In practice, of course, we do not follow each path to a final node and then
 * find out whether a path leads to a valid triangulation, but rather use an
 * inductive argument: if for all previously added cells the triangulation is
 * a valid one, then we can find out whether a path through the tree can yield
 * a valid triangulation by checking whether entering the present cell would
 * introduce any faces that have a nonunique direction; if that is so, then we
 * can stop following all paths below this point and track back immediately.
 *
 * Nevertheless, it is already obvious that the tree has $4^N$ leaves
 * in two space dimensions, since each of the $N$ cells can be added in four
 * orientations. Most of these nodes can be discarded rapidly, since firstly
 * the orientation of the first cell is irrelevant, and secondly if we add one
 * cell that has a neighbor that has already been added, then there are
 * already only two possible orientations left, so the total number of checks
 * we have to make until we find a valid way is significantly smaller than
 * $4^N$. However, the algorithm is still exponential in time and
 * linear in memory (we only have to store the information for the present
 * path in form of a stack of orientations of cells that have already been
 * added).
 *
 * In fact, the two examples above show that the exponential estimate is not a
 * pessimistic one: we indeed have to track back to one of the very first cells
 * there to find a way to add all cells in a consistent fashion.
 *
 * This discouraging situation is greatly improved by the fact that we have an
 * alternative algorithm for 2d that is always linear in runtime (discovered
 * and implemented by Michael Anderson of TICAM, University of Texas, in
 * 2003), and that for 3d we can find an algorithm that in practice is usually
 * only roughly linear in time and memory. We will describe these algorithms
 * in the following. A full description and theoretical analysis is given in
 * @cite AABB17 .
 *
 *
 * <h3>The 2d linear complexity algorithm</h3>
 *
 * The algorithm uses the fact that opposite faces of a cell need to have the
 * same orientation. So you start with one arbitrary line, choose an
 * orientation. Then the orientation of the opposite face is already fixed.
 * Then go to the two cells across the two faces we have fixed: for them, one
 * face is fixed, so we can also fix the opposite face. Go on with doing so.
 * Eventually, we have done this for a string of cells. Then take one of the
 * non-fixed faces of a cell which has already two fixed faces and do all this
 * again.
 *
 * In more detail, the algorithm is best illustrated using an example. We
 * consider the mesh below:
 * @verbatim
 *   9------10-------11
 *   |      |        /|
 *   |      |       / |
 *   |      |      /  |
 *   6------7-----8   |
 *   |      |     |   |
 *   |      |     |   |
 *   |      |     |   |
 *   3------4-----5   |
 *   |      |      \  |
 *   |      |       \ |
 *   |      |        \|
 *   0------1---------2
 * @endverbatim
 * First a cell is chosen ( (0,1,3,4) in this case). A single side of the cell
 * is oriented arbitrarily (3->4). This choice of orientation is then
 * propagated through the mesh, across sides and elements. (0->1), (6->7) and
 * (9->10). The involves edge-hopping and face hopping, giving a path through
 * the mesh shown in dots.
 * @verbatim
 *   9-->--10-------11
 *   |  .  |        /|
 *   |  .  |       / |
 *   |  .  |      /  |
 *   6-->--7-----8   |
 *   |  .  |     |   |
 *   |  .  |     |   |
 *   |  .  |     |   |
 *   3-->--4-----5   |
 *   |  .  |      \  |
 *   |  X  |       \ |
 *   |  .  |        \|
 *   0-->--1---------2
 * @endverbatim
 * This is then repeated for the other sides of the chosen element, orienting
 * more sides of the mesh.
 * @verbatim
 *   9-->--10-------11
 *   |     |        /|
 *   v.....v.......V |
 *   |     |      /. |
 *   6-->--7-----8 . |
 *   |     |     | . |
 *   |     |     | . |
 *   |     |     | . |
 *   3-->--4-----5 . |
 *   |     |      \. |
 *   ^..X..^.......^ |
 *   |     |        \|
 *   0-->--1---------2
 * @endverbatim
 * Once an element has been completely oriented it need not be considered
 * further. These elements are filled with o's in the diagrams. We then move
 * to the next element.
 * @verbatim
 *   9-->--10->-----11
 *   | ooo |  .     /|
 *   v ooo v  .    V |
 *   | ooo |  .   /  |
 *   6-->--7-->--8   |
 *   |     |  .  |   |
 *   |     |  .  |   |
 *   |     |  .  |   |
 *   3-->--4-->--5   |
 *   | ooo |  .   \  |
 *   ^ ooo ^  X    ^ |
 *   | ooo |  .     \|
 *   0-->--1-->------2
 * @endverbatim
 * Repeating this gives
 * @verbatim
 *   9-->--10->-----11
 *   | ooo | oooooo /|
 *   v ooo v ooooo V |
 *   | ooo | oooo /  |
 *   6-->--7-->--8   |
 *   |     |     |   |
 *   ^.....^..X..^...^
 *   |     |     |   |
 *   3-->--4-->--5   |
 *   | ooo | oooo \  |
 *   ^ ooo ^ ooooo ^ |
 *   | ooo | oooooo \|
 *   0-->--1-->------2
 * @endverbatim
 * and the final oriented mesh is
 * @verbatim
 *   9-->--10->-----11
 *   |     |        /|
 *   v     v       V |
 *   |     |      /  |
 *   6-->--7-->--8   |
 *   |     |     |   |
 *   ^     ^     ^   ^
 *   |     |     |   |
 *   3-->--4-->--5   |
 *   |     |      \  |
 *   ^     ^       ^ |
 *   |     |        \|
 *   0-->--1-->-------2
 * @endverbatim
 * It is obvious that this algorithm has linear run-time, since it only ever
 * touches each face exactly once.
 *
 * The algorithm just described in the two-dimensional case is
 * implemented for both 2d and (in generalized form) for 3d in this
 * class. The 3d case uses sheets instead of strings of cells to work
 * on. If a grid is orientable, then the algorithm is able to do its
 * work in linear time; if it is not orientable, then it aborts in
 * linear time as well.
 *
 * Both algorithms are described in the paper "On orienting edges of
 * unstructured two- and three-dimensional meshes", R. Agelek,
 * M. Anderson, W.  Bangerth, W. L. Barth, ACM Transactions on
 * Mathematical Software, vol. 44, article 5, 2017. A preprint is
 * available as <a href="https://arxiv.org/abs/1512.02137">arxiv
 * 1512.02137</a>.
 *
 *
 * <h3>For the curious</h3>
 *
 * Prior to the implementation of the algorithms described above (originally
 * implemented by Michael Anderson in 2002, and re-implemented by Wolfgang
 * Bangerth in 2016 based on the work in @cite AABB17),
 * we used a branch-and-cut algorithm initially
 * implemented in 2000 by Wolfgang Bangerth. Although it is no longer used,
 * here is how it works, and why it doesn't always work for large meshes since
 * its run-time can be exponential in bad cases.
 *
 * The first observation is that although there are counterexamples, problems
 * are usually local. For example, in the second example mentioned above, if
 * we had numbered the cells in a way that neighboring cells have similar cell
 * numbers, then the amount of backtracking needed is greatly reduced.
 * Therefore, in the implementation of the algorithm, the first step is to
 * renumber the cells in a Cuthill-McKee fashion: start with the cell with the
 * least number of neighbors and assign to it the cell number zero. Then find
 * all neighbors of this cell and assign to them consecutive further numbers.
 * Then find their neighbors that have not yet been numbered and assign to
 * them numbers, and so on. Graphically, this represents finding zones of
 * cells consecutively further away from the initial cells and number them in
 * this front-marching way. This already greatly improves locality of problems
 * and consequently reduced the necessary amount of backtracking.
 *
 * The second point is that we can use some methods to prune the tree, which
 * usually lead to a valid orientation of all cells very quickly.
 *
 * The first such method is based on the observation that if we fail to insert
 * one cell with number N, then this may not be due to cell N-1 unless N-1 is
 * a direct neighbor of N. The reason is obvious: the chosen orientation of
 * cell M could only affect the possibilities to add cell N if either it were
 * a direct neighbor or if there were a sequence of cells that were added
 * after M and that connected cells M and N. Clearly, for M=N-1, the latter
 * cannot be the case. Conversely, if we fail to add cell N, then it is not
 * necessary to track back to cell N-1, but we can track back to the neighbor
 * of N with the largest cell index and which has already been added.
 *
 * Unfortunately, this method can fail to yield a valid path through the tree
 * if not applied with care. Consider the following situation, initially
 * extracted from a mesh of 950 cells generated automatically by the program
 * BAMG (this program usually generates meshes that are quite badly balanced,
 * often have many -- sometimes 10 or more -- neighbors of one vertex, and
 * exposed several problems in the initial algorithm; note also that the
 * example is in 2d where we now have the much better algorithm described
 * above, but the same observations also apply to 3d):
 * @verbatim
 * 13----------14----15
 * | \         |     |
 * |  \    4   |  5  |
 * |   \       |     |
 * |    12-----10----11
 * |     |     |     |
 * |     |     |  7  |
 * |     |     |     |
 * |  3  |     8-----9
 * |     |     |     |
 * |     |     |  6  |
 * |     |     |     |
 * 4-----5-----6-----7
 * |     |     |     |
 * |  2  |  1  |  0  |
 * |     |     |     |
 * 0-----1-----2-----3
 * @endverbatim
 * Note that there is a hole in the middle. Assume now that the user described
 * the first cell 0 by the vertex numbers {2, 3, 6, 7}, and cell 5 by
 * {15, 14, 11, 10}, and assume that cells 1, 2, 3, and 4 are numbered
 * such that 5 can be added in initial rotation. All other cells are numbered
 * in the usual way, i.e. starting at the bottom left and counting
 * counterclockwise. Given this description of cells, the algorithm will start
 * with cell zero and add one cell after the other, up until the sixth one.
 * Then the situation will be the following:
 * @verbatim
 * 13----->---14--<--15
 * | \         |     |
 * |  >    4   v  5  v
 * |   \       |     |
 * |    12->--10--<--11
 * |     |     |     |
 * ^     |     |  7  |
 * |     |     |     |
 * |  3  ^     8-->--9
 * |     |     |     |
 * |     |     ^  6  ^
 * |     |     |     |
 * 4-->--5-->--6-->--7
 * |     |     |     |
 * ^  2  ^  1  ^  0  ^
 * |     |     |     |
 * 0-->--1-->--2-->--3
 * @endverbatim
 * Coming now to cell 7, we see that the two opposite lines at its top and
 * bottom have different directions; we will therefore find no orientation of
 * cell 7 in which it can be added without violation of the consistency of the
 * triangulation. According to the rule stated above, we track back to the
 * neighbor with greatest index, which is cell 6, but since its bottom line is
 * to the right, its top line must be to the right as well, so we won't be
 * able to find an orientation of cell 6 such that 7 will fit into the
 * triangulation. Then, if we have finished all possible orientations of cell
 * 6, we track back to the neighbor of 6 with the largest index and which has
 * been added already. This would be cell 0. However, we know that the
 * orientation of cell 0 can't be important, so we conclude that there is no
 * possible way to orient all the lines of the given cells such that they
 * satisfy the requirements of deal.II triangulations. We know that this can't
 * be, so it results in an exception be thrown.
 *
 * The bottom line of this example is that when we looked at all possible
 * orientations of cell 6, we couldn't find one such that cell 7 could be
 * added, and then decided to track back to cell 0. We did not even attempt to
 * turn cell 5, after which it would be simple to add cell 7. Thus, the
 * algorithm described above has to be modified: we are only allowed to track
 * back to that neighbor that has already been added, with the largest cell
 * index, if we fail to add a cell in any orientation. If we track back
 * further because we have exhausted all possible orientations but could add
 * the cell (i.e. we track back since another cell, further down the road
 * couldn't be added, irrespective of the orientation of the cell which we are
 * presently considering), then we are not allowed to track back to one of its
 * neighbors, but have to track back only one cell index.
 *
 * The second method to prune the tree is that usually we cannot add a new
 * cell since the orientation of one of its neighbors that have already been
 * added is wrong. Thus, if we may try to rotate one of the neighbors (of
 * course making sure that rotating that neighbor does not violate the
 * consistency of the triangulation) in order to allow the present cell to be
 * added.
 *
 * While the first method could be explained in terms of backtracking in the
 * tree of orientations more than one step at once, turning a neighbor means
 * jumping to a totally different place in the tree. For both methods, one can
 * find arguments that they will never miss a path that is valid and only skip
 * paths that are invalid anyway.
 *
 * These two methods have proven extremely efficient. We have been able to
 * read very large grids (several ten thousands of cells) without the need to
 * track back much. In particular, the time to find an ordering of the cells
 * was found to be mostly linear in the number of cells, and the time to
 * reorder them is usually much smaller (for example by one order of
 * magnitude) than the time needed to read the data from a file, and also to
 * actually generate the triangulation from this data using the
 * Triangulation::create_triangulation() function.
 *
 * @ingroup grid
 */
