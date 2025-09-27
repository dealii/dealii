// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_grid_cell_data_h
#define dealii_grid_cell_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <deal.II/grid/reference_cell.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * The CellData class (and the related SubCellData class) is used to
 * provide a comprehensive, but minimal, description of the cells when
 * creating a triangulation via Triangulation::create_triangulation().
 * Specifically, each CellData object -- describing one cell in a
 * triangulation -- has member variables for indices of the $2^d$ vertices
 * (the actual coordinates of the vertices are described in a separate
 * vector passed to Triangulation::create_triangulation(), so the CellData
 * object only needs to store indices into that vector), the material
 * id of the cell that can be used in applications to describe which
 * part of the domain a cell belongs to (see
 * @ref GlossMaterialId "the glossary entry on material ids"),
 * and a manifold id that is used to describe the geometry object
 * that is responsible for this cell (see
 * @ref GlossManifoldIndicator "the glossary entry on manifold ids")
 * to describe the manifold this object belongs to.
 *
 * This structure is also used to represent data for faces and edges when used
 * as a member of the SubCellData class. In this case, the template argument
 * @p structdim of an object will be less than the dimension @p dim of the
 * triangulation. If this is so, then #vertices array represents the indices of
 * the vertices of one face or edge of one of the cells passed to
 * Triangulation::create_triangulation(). Furthermore, for faces the
 * material id has no meaning, and the @p material_id field is reused
 * to store a @p boundary_id instead to designate which part of the boundary
 * the face or edge belongs to (see
 * @ref GlossBoundaryIndicator "the glossary entry on boundary ids").
 *
 * An example showing how this class can be used is in the
 * <code>create_coarse_grid()</code> function of step-14. There are also
 * many more use cases in the implementation of the functions of the
 * GridGenerator namespace.
 *
 * @ingroup grid
 */
template <int structdim>
struct CellData
{
  /**
   * Indices of the vertices of this cell. These indices correspond
   * to entries in the vector of vertex locations passed to
   * Triangulation::create_triangulation().
   *
   * By default, the constructor of this class initializes this variable to
   * have as many entries as it takes to describe a hypercube cell (i.e., a
   * ReferenceCells::Line, ReferenceCells::Quadrilateral, or
   * ReferenceCells::Hexahedron). This is historical and dates back to the time
   * where deal.II could only deal with these kinds of cells. If you want an
   * object of the current type to describe, for example, a triangle or
   * tetrahedron, then you either have to call this constructor with an explicit
   * argument different from the default value, or manually resize the
   * `vertices` member variable after construction.
   *
   * The kind of cell described by the current object is then determined by
   * calling ReferenceCell::n_vertices_to_type() on the number of vertices
   * described by this array.
   */
  std::vector<unsigned int> vertices;

  /**
   * Material or boundary indicator of this cell.
   * This field is a union that stores <i>either</i> a boundary or
   * a material id, depending on whether the current object is used
   * to describe a cell (in a vector of CellData objects) or a
   * face or edge (as part of a SubCellData object).
   */
  union
  {
    /**
     * The material id of the cell being described. See the documentation
     * of the CellData class for examples of how to use this field.
     *
     * This variable can only be used if the current object is used to
     * describe a cell, i.e., if @p structdim equals the dimension
     * @p dim of a triangulation.
     */
    types::material_id material_id;

    /**
     * The boundary id of a face or edge being described. See the documentation
     * of the CellData class for examples of how to use this field.
     *
     * This variable can only be used if the current object is used to
     * describe a face or edge, i.e., if @p structdim is less than the dimension
     * @p dim of a triangulation. In this case, the CellData object this
     * variable belongs to will be part of a SubCellData object.
     */
    types::boundary_id boundary_id;
  };

  /**
   * Manifold identifier of this object. This identifier should be used to
   * identify the manifold to which this object belongs, and from which this
   * object will collect information on how to add points upon refinement.
   *
   * See the documentation of the CellData class for examples of how to use
   * this field.
   */
  types::manifold_id manifold_id;

  /**
   * Default constructor. Sets the member variables to the following values:
   *
   * - vertex indices to invalid values
   * - boundary or material id zero (the default for boundary or material ids)
   * - manifold id to numbers::flat_manifold_id.
   *
   * By default, the constructor initializes the `vertices` member variable to
   * have as many entries as it takes to describe a hypercube cell (i.e., a
   * ReferenceCells::Line, ReferenceCells::Quadrilateral, or
   * ReferenceCells::Hexahedron). This is historical and dates back to the time
   * where deal.II could only deal with these kinds of cells. If you want an
   * object of the current type to describe, for example, a triangle or
   * tetrahedron, then you either have to call this constructor with an explicit
   * argument different from the default value, or manually resize the
   * `vertices` member variable after construction.
   */
  CellData(const unsigned int n_vertices =
             ReferenceCells::get_hypercube<structdim>().n_vertices());

  /**
   * Comparison operator.
   */
  bool
  operator==(const CellData<structdim> &other) const;

  /**
   * Read or write the data of this object to or from a stream for the
   * purpose of serialization using the [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  static_assert(structdim > 0,
                "The class CellData can only be used for structdim>0.");
};



/**
 * The SubCellData class is used to describe information about faces and
 * edges at the boundary of a mesh when creating a triangulation via
 * Triangulation::create_triangulation(). It contains member variables
 * that describe boundary edges and boundary quads.
 *
 * The class has no template argument and is used both in the description
 * of boundary edges in 2d (in which case the contents of the
 * @p boundary_quads member variable are ignored), as well as in the
 * description of boundary edges and faces in 3d (in which case both the
 * @p boundary_lines and @p boundary_quads members may be used). It is also
 * used as the argument to Triangulation::create_triangulation() in 1d,
 * where the contents of objects of the current type are simply ignored.
 *
 * By default, Triangulation::create_triangulation() simply assigns
 * default boundary indicators and manifold indicators to edges and
 * quads at the boundary of the mesh. (See the glossary entries on
 * @ref GlossBoundaryIndicator "boundary ids"
 * and
 * @ref GlossManifoldIndicator "manifold ids"
 * for more information on what they represent.) As a consequence,
 * it is not <i>necessary</i> to explicitly describe the properties
 * of boundary objects. In all cases, these properties can also be
 * set at a later time, once the triangulation has already been
 * created. On the other hand, it is sometimes convenient to describe
 * boundary indicators or manifold ids at the time of creation. In
 * these cases, the current class can be used by filling the
 * @p boundary_lines and @p boundary_quads vectors with
 * CellData<1> and CellData<2> objects that correspond to boundary
 * edges and quads for which properties other than the default
 * values should be used.
 *
 * Each entry in the @p boundary_lines and @p boundary_quads vectors
 * then needs to correspond to an edge or quad of the cells that
 * are described by the vector of CellData objects passed to
 * Triangulation::create_triangulation(). I.e., the vertex indices
 * stored in each entry need to correspond to an edge or face
 * of the triangulation that has the same set of vertex indices,
 * and in the same order. For these boundary edges or quads, one can
 * then set either or both the CellData::boundary_id and
 * CellData::manifold_id.
 *
 * There are also use cases where one may want to set the manifold id
 * of an <i>interior</i> edge or face. Such faces, identified by
 * their vertex indices, may also appear in the
 * @p boundary_lines and @p boundary_quads vectors (despite the names of
 * these member variables). However, it is then obviously not allowed
 * to set a boundary id (because the object is not actually part of
 * the boundary). As a consequence, to be valid, the CellData::boundary_id
 * of interior edges or faces needs to equal
 * numbers::internal_face_boundary_id.
 *
 * @ingroup grid
 */
struct SubCellData
{
  /**
   * A vector of CellData<1> objects that describe boundary and manifold
   * information for edges of 2d or 3d triangulations. For 2d triangulations,
   * edges (lines) are of course the faces of the cells.
   *
   * This vector must not be used in the creation of 1d triangulations.
   */
  std::vector<CellData<1>> boundary_lines;

  /**
   * A vector of CellData<2> objects that describe boundary and manifold
   * information for triangles and quads of 3d triangulations. The name of
   * the variable is historical and dates back to a time when deal.II only
   * supported hexahedral meshes in 3d, where then all boundary faces were
   * necessarily quadrilaterals. However, the variable is now also used to
   * describe boundary triangles for tetrahedral and mixed meshes.
   *
   * Whether an element in this array describes a boundary triangle or a
   * boundary quad is determined by how many elements the `vertices`
   * member variable of `CellData<2>` stores.
   *
   * This vector must not be used in the creation of 1d or 2d triangulations.
   */
  std::vector<CellData<2>> boundary_quads;

  /**
   * Determine whether the member variables above which may not be used in a
   * given dimension are really empty. In other words, this function returns
   * whether
   * both @p boundary_lines and @p boundary_quads are empty vectors
   * when @p dim equals one, and whether the @p boundary_quads
   * vector is empty when @p dim equals two.
   */
  bool
  check_consistency(const unsigned int dim) const;
};


template <int structdim>
template <class Archive>
void
CellData<structdim>::serialize(Archive &ar, const unsigned int /*version*/)
{
  ar &vertices;
  ar &material_id;
  ar &boundary_id;
  ar &manifold_id;
}

DEAL_II_NAMESPACE_CLOSE

#endif
