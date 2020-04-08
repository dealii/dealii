// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_grid_construction_utilities_h
#define dealii_grid_construction_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria.h>


DEAL_II_NAMESPACE_OPEN

/*------------------------------------------------------------------------*/

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
   */
  unsigned int vertices[GeometryInfo<structdim>::vertices_per_cell];

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
   * - manifold id to numbers::flat_manifold_id
   */
  CellData();

  /**
   * Comparison operator.
   */
  bool
  operator==(const CellData<structdim> &other) const;

  /**
   * Boost serialization function
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
   * information for edges of 2d or 3d triangulations.
   *
   * This vector may not be used in the creation of 1d triangulations.
   */
  std::vector<CellData<1>> boundary_lines;

  /**
   * A vector of CellData<2> objects that describe boundary and manifold
   * information for quads of 3d triangulations.
   *
   * This vector may not be used in the creation of 1d or 2d triangulations.
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

/**
 * A namespace dedicated to the struct Description, which can be used in
 * Triangulation::create_triangulation().
 */
namespace TriangulationDescription
{
  /**
   * Configuration flags for Triangulations.
   * Settings can be combined using bitwise OR.
   *
   * @author Peter Munch, 2019
   */
  enum Settings
  {
    /**
     * Default settings, other options are disabled.
     */
    default_setting = 0x0,
    /**
     * This flag needs to be set to use the geometric multigrid
     * functionality. This option requires additional computation and
     * communication.
     */
    construct_multigrid_hierarchy = 0x1
  };

  /**
   * Information needed for each locally relevant cell, stored in
   * Description and used during construction of a
   * Triangulation. This struct stores
   * the cell id, the subdomain_id and the level_subdomain_id as well as
   * information related to manifold_id and boundary_id.
   *
   * @note Similarly to dealii::CellData, this structure stores information
   * about a cell. However, in contrast to dealii::CellData, it also stores
   * a unique id, partitioning information, and information related to cell
   * faces and edges.
   *
   * @author Peter Munch, 2019
   */
  template <int dim>
  struct CellData
  {
    /**
     * Boost serialization function
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int /*version*/);

    /**
     * Comparison operator.
     */
    bool
    operator==(const CellData<dim> &other) const;

    /**
     * Unique CellID of the cell.
     */
    CellId::binary_type id;

    /**
     * subdomain_id of the cell.
     */
    types::subdomain_id subdomain_id;

    /**
     * level_subdomain_id of the cell.
     */
    types::subdomain_id level_subdomain_id;

    /**
     * Manifold id of the cell.
     */
    types::manifold_id manifold_id;

    /**
     * Manifold id of all lines of the cell.
     *
     * @note Only used for 2D and 3D.
     */
    std::array<types::manifold_id, GeometryInfo<dim>::lines_per_cell>
      manifold_line_ids;

    /**
     * Manifold id of all face quads of the cell.
     *
     * @note Only used for 3D.
     */
    std::array<types::manifold_id,
               dim == 1 ? 1 : GeometryInfo<3>::quads_per_cell>
      manifold_quad_ids;

    /**
     * List of face number and boundary id of all non-internal faces of the
     * cell.
     */
    std::vector<std::pair<unsigned int, types::boundary_id>> boundary_ids;
  };

  /**
   * Data used in Triangulation::create_triangulation().
   *
   * @author Peter Munch, 2019
   */
  template <int dim, int spacedim>
  struct Description
  {
    /**
     * Boost serialization function
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int /*version*/);

    /**
     * Comparison operator.
     */
    bool
    operator==(const Description<dim, spacedim> &other) const;

    /**
     * Cells of the locally-relevant coarse-grid triangulation.
     */
    std::vector<dealii::CellData<dim>> coarse_cells;

    /**
     * Vertices of the locally-relevant coarse-grid triangulation.
     */
    std::vector<Point<spacedim>> coarse_cell_vertices;

    /**
     * List that for each locally-relevant coarse cell provides the
     * corresponding global @ref GlossCoarseCellId.
     */
    std::vector<types::coarse_cell_id> coarse_cell_index_to_coarse_cell_id;

    /**
     * CellData for each locally relevant cell on each level. cell_infos[i]
     * contains the CellData for each locally relevant cell on the ith
     * level.
     */
    std::vector<std::vector<CellData<dim>>> cell_infos;

    /**
     * The MPI communicator used to create this struct. It will be compared
     * to the communicator inside of the Triangulation
     * and an assert is thrown if they do not match.
     *
     * @note Please note this is necessary since the communicator inside of
     * parallel::TriangulationBase is const and cannot be changed after the
     * constructor has been called.
     *
     */
    MPI_Comm comm;

    /**
     * Properties to be use in the construction of the triangulation.
     */
    Settings settings;
  };


  /**
   * A namespace for TriangulationDescription::Description utility functions.
   *
   * @ingroup TriangulationDescription
   */
  namespace Utilities
  {
    /**
     * Construct TriangulationDescription::Description from a given
     * partitioned triangulation `tria` and a specified process.
     * The input triangulation can be either
     * a serial triangulation of type dealii::Triangulation which has been
     * colored (subdomain_id and/or level_subdomain_id has been set) or a
     * distributed triangulation of type
     * dealii::parallel::distributed::Triangulation, where the partitioning is
     * adopted unaltered.
     *
     * @param tria Partitioned input triangulation.
     * @param comm MPI_Communicator to be used. In the case
     *   of dealii::parallel::distributed::Triangulation, the communicators have
     * to match.
     * @param construct_multilevel_hierarchy Signal if the multigrid levels
     *        should be constructed.
     * @param my_rank_in Construct Description for the specified rank (only
     *   working for serial triangulations that have been partitioned by
     *   functions like GridToold::partition_triangulation()).
     * @return Description to be used to set up a Triangulation.
     *
     * @note Multilevel hierarchies are supported if it is enabled in
     *   parallel::fullydistributed::Triangulation.
     *
     * @note Hanging nodes in the input triangulation are supported. However,
     *   to be able to use this feature in the case of
     *   parallel::fullydistributed::Triangulation, the user has to enable
     *   multilevel hierarchy support in
     *   parallel::fullydistributed::Triangulation.
     *
     * @author Peter Munch, 2019
     */
    template <int dim, int spacedim = dim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const dealii::Triangulation<dim, spacedim> &tria,
      const MPI_Comm                              comm,
      const bool         construct_multilevel_hierarchy = false,
      const unsigned int my_rank_in = numbers::invalid_unsigned_int);


    /**
     * Construct a TriangulationDescription::Description. In contrast
     * to the function above, this function is also responsible for creating
     * a serial triangulation and for its partitioning (by calling the
     * provided `std::function` objects). Internally only selected processes (
     * every n-th/each root of a group of size group_size) create a serial
     * triangulation and the TriangulationDescription::Description for all
     * processes in its group, which is communicated.
     *
     * @note A reasonable group size is the size of a NUMA domain or the
     * size of a compute node.
     *
     * @param serial_grid_generator A function which creates a serial triangulation.
     * @param serial_grid_partitioner A function which can partition a serial
     *   triangulation, i.e., sets the sudomain_ids of the active cells.
     *   The function takes as the first argument a serial triangulation,
     *   as the second argument the MPI communicator, and as the third
     *   argument the group size.
     * @param comm MPI communicator.
     * @param group_size The size of each group.
     * @param construct_multilevel_hierarchy Construct multigrid levels.
     * @return Description to be used to set up a Triangulation.
     *
     * @author Peter Munch, 2019
     */
    template <int dim, int spacedim = dim>
    Description<dim, spacedim>
    create_description_from_triangulation_in_groups(
      const std::function<void(dealii::Triangulation<dim, spacedim> &)>
        &                                            serial_grid_generator,
      const std::function<void(dealii::Triangulation<dim, spacedim> &,
                               const MPI_Comm,
                               const unsigned int)> &serial_grid_partitioner,
      const MPI_Comm                                 comm,
      const int                                      group_size = 1,
      const bool construct_multilevel_hierarchy                 = false);

  } // namespace Utilities



  template <int dim>
  template <class Archive>
  void
  CellData<dim>::serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &id;
    ar &subdomain_id;
    ar &level_subdomain_id;
    ar &manifold_id;
    if (dim >= 2)
      ar &manifold_line_ids;
    if (dim >= 3)
      ar &manifold_quad_ids;
    ar &boundary_ids;
  }


  template <int dim, int spacedim>
  template <class Archive>
  void
  Description<dim, spacedim>::serialize(Archive &ar,
                                        const unsigned int /*version*/)
  {
    ar &coarse_cells;
    ar &coarse_cell_vertices;
    ar &coarse_cell_index_to_coarse_cell_id;
    ar &cell_infos;
    ar &settings;
  }



  template <int dim>
  bool
  CellData<dim>::operator==(const CellData<dim> &other) const
  {
    if (this->id != other.id)
      return false;
    if (this->subdomain_id != other.subdomain_id)
      return false;
    if (this->level_subdomain_id != other.level_subdomain_id)
      return false;
    if (this->manifold_id != other.manifold_id)
      return false;
    if (dim >= 2 && this->manifold_line_ids != other.manifold_line_ids)
      return false;
    if (dim >= 3 && this->manifold_quad_ids != other.manifold_quad_ids)
      return false;
    if (this->boundary_ids != other.boundary_ids)
      return false;

    return true;
  }



  template <int dim, int spacedim>
  bool
  Description<dim, spacedim>::
  operator==(const Description<dim, spacedim> &other) const
  {
    if (this->coarse_cells != other.coarse_cells)
      return false;
    if (this->coarse_cell_vertices != other.coarse_cell_vertices)
      return false;
    if (this->coarse_cell_index_to_coarse_cell_id !=
        other.coarse_cell_index_to_coarse_cell_id)
      return false;
    if (this->cell_infos != other.cell_infos)
      return false;
    if (this->settings != other.settings)
      return false;

    return true;
  }
} // namespace TriangulationDescription


DEAL_II_NAMESPACE_CLOSE

#endif
