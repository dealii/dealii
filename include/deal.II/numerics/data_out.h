// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_data_out_h
#define dealii_data_out_h



#include <deal.II/base/config.h>

#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/numerics/data_out_dof_data.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DataOutImplementation
  {
    /**
     * A derived class for use in the DataOut class. This is a class for the
     * AdditionalData kind of data structure discussed in the documentation of
     * the WorkStream context.
     */
    template <int dim, int spacedim>
    struct ParallelData : public ParallelDataBase<dim, spacedim>
    {
      ParallelData(
        const unsigned int               n_datasets,
        const unsigned int               n_subdivisions,
        const std::vector<unsigned int> &n_postprocessor_outputs,
        const dealii::hp::MappingCollection<dim, spacedim> &mapping,
        const std::vector<
          std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                                                     &finite_elements,
        const UpdateFlags                             update_flags,
        const std::vector<std::vector<unsigned int>> &cell_to_patch_index_map);

      std::vector<Point<spacedim>> patch_evaluation_points;

      const std::vector<std::vector<unsigned int>> *cell_to_patch_index_map;
    };
  } // namespace DataOutImplementation
} // namespace internal



/**
 * This class is the main class to provide output of data described by finite
 * element fields defined on a collection of cells.
 *
 * This class is an actual implementation of the functionality proposed by the
 * DataOut_DoFData class. It offers a function build_patches() that generates
 * the data to be written in some graphics format. Most of the interface and
 * an example of its use is described in the documentation of this base class.
 *
 * The only thing this class offers is the function build_patches() which
 * loops over all cells of the triangulation stored by the
 * attach_dof_handler() function of the base class (with the exception of
 * cells of parallel::distributed::Triangulation objects that are not owned by
 * the current processor) and converts the data on these to actual patches
 * which are the objects that are later output by the functions of the base
 * classes. You can give a parameter to the function which determines how many
 * subdivisions in each coordinate direction are to be performed, i.e. of how
 * many subcells each patch shall consist. The default is one, but you may want
 * to choose a higher number for higher order elements, for example two for
 * quadratic elements, three for cubic elements, and so on. (See
 * step-11 for an example.) The purpose
 * of this parameter is because most graphics programs do not allow to specify
 * higher order polynomial functions in the file formats: only data at
 * vertices can be plotted and is then shown as a bilinear interpolation
 * within the interior of cells. This may be insufficient if you have higher
 * order finite elements, and the only way to achieve better output is to
 * subdivide each cell of the mesh into several cells for graphical output. Of
 * course, what you get to see is still a bilinear interpolation on each cell
 * of the output (where these cells are not subdivisions of the cells of the
 * triangulation in use) due to the same limitations in output formats, but at
 * least a bilinear interpolation of a higher order polynomial on a finer
 * mesh.
 *
 * Note that after having called build_patches() once, you can call one or
 * more of the write() functions of DataOutInterface. You can therefore output
 * the same data in more than one format without having to rebuild the
 * patches.
 *
 *
 * <h3>User interface information</h3>
 *
 * The base classes of this class, DataOutBase, DataOutInterface and
 * DataOut_DoFData offer several interfaces of their own. Refer to the
 * DataOutBase class's documentation for a discussion of the different output
 * formats presently supported, DataOutInterface for ways of selecting which
 * format to use upon output at run-time and without the need to adapt your
 * program when new formats become available, as well as for flags to
 * determine aspects of output. The DataOut_DoFData() class's documentation
 * has an example of using nodal data to generate output.
 *
 *
 * <h3>Extensions</h3>
 *
 * By default, this class produces patches for all active cells. Sometimes,
 * this is not what you want, maybe because there are simply too many (and too
 * small to be seen individually) or because you only want to see a certain
 * region of the domain (for example only in the fluid part of the domain in
 * step-46), or for some other reason.
 *
 * For this, internally build_patches() does not generate the sequence of cells
 * to be converted into patches itself, but relies on the two private
 * std::function objects first_cell_function() and next_cell_function(). By
 * default, they return the first active cell, and the next active cell,
 * respectively. But this can be changed using the set_cell_selection() function
 * that allows you to replace this behavior. What set_cell_selection() wants to
 * know is how you want to pick out the first cell on which output should be
 * generated, and how given one cell on which output is generated you want to
 * pick the next cell.
 *
 * This may,
 * for example, include only cells that are in parts of a domain (e.g., if you
 * don't care about the solution elsewhere, think for example a buffer region
 * in which you attenuate outgoing waves in the Perfectly Matched Layer
 * method) or if you don't want output to be generated at all levels of an
 * adaptively refined mesh because this creates too much data (in this case,
 * the set of cells returned by your implementations of the `first_cell` and
 * `next_cell` arguments to set_cell_selection() will include
 * non-active cells, and DataOut::build_patches()
 * will simply take interpolated values of the solution instead of the exact
 * values on these cells children for output).
 *
 * @ingroup output
 */
template <int dim, int spacedim = dim>
class DataOut : public DataOut_DoFData<dim, dim, spacedim, spacedim>
{
public:
  /**
   * Typedef to the iterator type of the dof handler class under
   * consideration.
   */
  using cell_iterator =
    typename DataOut_DoFData<dim, dim, spacedim, spacedim>::cell_iterator;

  /**
   * The type of the function object returning the first cell as used in
   * set_cell_selection().
   */
  using FirstCellFunctionType =
    typename std::function<cell_iterator(const Triangulation<dim, spacedim> &)>;

  /**
   * The type of the function object returning the next cell as used in
   * set_cell_selection().
   */
  using NextCellFunctionType =
    typename std::function<cell_iterator(const Triangulation<dim, spacedim> &,
                                         const cell_iterator &)>;

  /**
   * Enumeration describing the part of the domain in which cells
   * should be written with curved boundaries. In reality, no file
   * format we are aware of really supports curved boundaries, but
   * this can be emulated by plotting edges as a sequence of straight
   * lines (and faces in 3d as a collection of bilinear patches) if
   * DataOut::build_patches() is called with a number of subdivisions
   * greater than 1.
   *
   * The elements of this enumeration then describe for which cells
   * DataOut::build_patches() will query the manifold or boundary
   * description for curved geometries.
   */
  enum CurvedCellRegion
  {
    /**
     * The geometry or boundary description will never be queried for
     * curved geometries. This means that even if you have more than
     * one subdivision per cell (see DataOut::build_patches() for what
     * exactly this means) and even if the geometry really is curved,
     * each cell will still be subdivided as if it was just a bi- or
     * trilinear cell.
     */
    no_curved_cells,

    /**
     * The geometry or boundary description will be queried for curved
     * geometries for cells located at the boundary, i.e., for cells
     * that have at least one face at the boundary. This is sufficient
     * if you have not attached a manifold description to the
     * interiors of cells but only to faces at the boundary.
     */
    curved_boundary,

    /**
     * The geometry description will be queried for all cells and all
     * faces, whether they are at the boundary or not. This option is
     * appropriate if you have attached a manifold object to cells
     * (not only to boundary faces).
     */
    curved_inner_cells
  };

  /**
   * Constructor.
   */
  DataOut();

  /**
   * This is the central function of this class since it builds the list of
   * patches to be written by the low-level functions of the base class. A
   * patch is, in essence, some intermediate representation of the data on
   * each cell of a triangulation and DoFHandler object that can then be used
   * to write files in some format that is readable by visualization programs.
   *
   * You can find an overview of the use of this function in the general
   * documentation of this class. An example is also provided in the
   * documentation of this class's base class DataOut_DoFData.
   *
   * @param n_subdivisions A parameter that determines how many "patches" this
   * function will build out of every cell. If you do not specify this value
   * in calling, or provide the default value zero, then this is interpreted
   * as DataOutInterface::default_subdivisions which most of the time will be
   * equal to one (unless you have set it to something else). The purpose of
   * this parameter is to subdivide each cell of the mesh into $2\times 2,
   * 3\times 3, \ldots$ "patches" in 2d, and $2\times 2\times 2, 3\times
   * 3\times 3, \ldots$ (if passed the value 2, 3, etc) where each patch
   * represents the data from a regular subdivision of the cell into equal
   * parts. Most of the times, this is not necessary and outputting one patch
   * per cell is exactly what you want to plot the solution. That said, the
   * data we write into files for visualization can only represent (bi-,
   * tri)linear data on each cell, and most visualization programs can in fact
   * only visualize this kind of data. That's good enough if you work with
   * (bi-, tri)linear finite elements, in which case what you get to see is
   * exactly what has been computed. On the other hand, if you work with (bi-,
   * tri)quadratic elements, then what is written into the output file is just
   * a (bi-, tri)linear interpolation onto the current mesh, i.e., only the
   * values at the vertices. If this is not good enough, you can, for example,
   * specify @p n_subdivisions equal to 2 to plot the solution on a
   * once-refined mesh, or if set to 3, on a mesh where each cell is represented
   * by 3-by-3 patches. On each of these smaller patches, given the limitations
   * of output formats, the data is still linearly interpolated, but a linear
   * interpolation of quadratic data on a finer mesh is still a better
   * representation of the actual quadratic surface than on the original mesh.
   * In other words, using this parameter can not help you plot the solution
   * exactly, but it can get you closer if you use finite elements of higher
   * polynomial degree.
   *
   * @note Specifying `n_subdivisions>1` is useful when using higher order
   *   finite elements, but in general it does not actually result in the
   *   visualization showing higher order polynomial surfaces -- rather, you
   *   just get a (bi-, tri-)linear interpolation of that higher order
   *   surface on a finer mesh. However, when outputting the solution in the
   *   VTK and VTU file formats via DataOutInterface::write_vtk() or
   *   DataOutInterface::write_vtu() (where DataOutInterface is a base
   *   class of the current class) as we often do in the tutorials,
   *   you can provide a set of flags via the DataOutBase::VtkFlags
   *   structure that includes the
   *   DataOutBase::VtkFlags::write_higher_order_cells flag. When set, the
   *   subdivisions produced by this function will be interpreted as
   *   support points for a higher order polynomial that will then actually
   *   be visualized as such. This is shown in step-11, for example. It
   *   is worth noting, however, that this requires a
   *   sufficiently new version of one of the VTK-based visualization
   *   programs.
   */
  virtual void
  build_patches(const unsigned int n_subdivisions = 0);

  /**
   * Same as above, except that the additional first parameter defines a
   * mapping that is to be used in the generation of output. If
   * <tt>n_subdivisions>1</tt>, the points interior of subdivided patches
   * which originate from cells at the boundary of the domain can be computed
   * using the mapping, i.e., a higher order mapping leads to a representation
   * of a curved boundary by using more subdivisions. Some mappings like
   * MappingQEulerian result in curved cells in the interior of the domain.
   * The same is true if you have attached a manifold description to
   * the cells of a triangulation (see
   * @ref manifold "Manifolds"
   * for more information). However, there is no easy way to query the mapping
   * or manifold whether it really does lead to curved cells.
   * Thus the last argument @p curved_region takes one of three values
   * resulting in no curved cells at all, curved cells at the boundary
   * (default) or curved cells in the whole domain. For more information
   * about these three options, see the CurvedCellRegion enum's
   * description.
   *
   * Even for non-curved cells, the mapping argument can be used for
   * Eulerian mappings (see class MappingQ1Eulerian) where a mapping is used
   * not only to determine the position of points interior to a cell, but also
   * of the vertices.  It offers an opportunity to watch the solution on a
   * deformed triangulation on which the computation was actually carried out,
   * even if the mesh is internally stored in its undeformed configuration and
   * the deformation is only tracked by an additional vector that holds the
   * deformation of each vertex.
   */
  virtual void
  build_patches(const Mapping<dim, spacedim> &mapping,
                const unsigned int            n_subdivisions = 0,
                const CurvedCellRegion        curved_region  = curved_boundary);

  /**
   * Same as above, but for hp::MappingCollection.
   */
  virtual void
  build_patches(const hp::MappingCollection<dim, spacedim> &mapping,
                const unsigned int                          n_subdivisions = 0,
                const CurvedCellRegion curved_region = curved_boundary);

  /**
   * A function that allows selecting for which cells output should be
   * generated. This function takes two arguments, both `std::function`
   * objects that can be used what the first cell on which output is
   * generated is supposed to be, and what given one cell the next
   * function is supposed to be. Through these function objects,
   * it is possible to select a subset of cells on which output should
   * be produced (e.g., only selecting those cells that belong to a
   * part of the domain -- say, the fluid domain in a code such as step-46),
   * or to completely change *where* output is produced (e.g., to produce
   * output on non-active cells of a multigrid hierarchy or if the finest
   * level of a mesh is so fine that generating graphical output would lead
   * to an overwhelming amount of data).
   *
   * @param[in] first_cell A function object that takes as argument the
   *   triangulation this class works on and that should return the first cell
   *   on which output should be generated.
   * @param[in] next_cell A function object that takes as arguments the
   *   triangulation as well as the last cell
   *   on which output was generated, and that should return the next
   *   cell on which output should be generated. If there is no next
   *   cell, i.e., if the input argument to the `next_cell` function object
   *   is the last cell on which output is to be generated, then `next_cell`
   *   must return `triangulation.end()`.
   *
   * These function objects are not difficult to write, but also not immediately
   * obvious. As a consequence, there is a second variation of this function
   * that takes a IteratorFilter argument and generates the corresponding
   * functions itself.
   *
   * @note This function is also called in the constructor of this class,
   *   where the default behavior is set. By default, this class will select all
   *   @ref GlossLocallyOwnedCell "locally owned"
   *   and
   *   @ref GlossActive "active"
   *   cells for output.
   *
   * @note If you have cell data (in contrast to nodal, or dof, data) such as
   *   error indicators, then you must make sure that the `first_cell` and
   *   `next_cell` function objects only walk over active cells, since cell data
   *   cannot be interpolated to a coarser cell. If you do have cell data and
   *   use this pair of functions and they return a non-active cell, then an
   *   exception will be thrown.
   */
  void
  set_cell_selection(
    const std::function<cell_iterator(const Triangulation<dim, spacedim> &)>
                                                              &first_cell,
    const std::function<cell_iterator(const Triangulation<dim, spacedim> &,
                                      const cell_iterator &)> &next_cell);

  /**
   * A variation of the previous function that selects a subset of all
   * cells for output based on the filter encoded in the FilteredIterator
   * object given as argument. A typical way to generate the argument
   * is via the make_filtered_iterator() function.
   *
   * Alternatively, since FilteredIterator objects can be created from
   * just a predicate (i.e., a function object that returns a `bool`), it is
   * possible to call this function with just a lambda function, which will then
   * automatically be converted to a FilteredIterator object. For example, the
   * following piece of code works:
   * @code
   *   DataOut<dim> data_out;
   *   data_out.set_cell_selection(
   *          [](const typename Triangulation<dim>::cell_iterator &cell) {
   *              return (cell->is_active() && cell->subdomain_id() == 0);
   *          });
   * @endcode
   * In this case, the lambda function selects all of those cells that are
   * @ref GlossActive "active"
   * and whose subdomain id is zero. These will then be the only cells on
   * which output is generated.
   *
   * @note Not all filters will result in subsets of cells for which
   *   output can actually be generated. For example, if you are working
   *   on parallel meshes where data is only available on some cells,
   *   then you better make sure that your `filtered_iterator` only
   *   loops over the
   *   @ref GlossLocallyOwnedCell "locally owned"
   *   cells; likewise, in most cases you will probably only want to work on
   *   @ref GlossActive "active"
   *   cells since this is where the solution actually lives. In particular,
   *   if you have added vectors that represent data defined on cells
   *   (instead of nodal data), then you can not generate output on non-active
   *   cells and your iterator filter should reflect this.
   */
  void
  set_cell_selection(const FilteredIterator<cell_iterator> &filtered_iterator);

  /**
   * Return the two function objects that are in use for determining the first
   * and the next cell as set by set_cell_selection().
   */
  std::pair<FirstCellFunctionType, NextCellFunctionType>
  get_cell_selection() const;

private:
  /**
   * A function object that is used to select what the first cell is going to
   * be on which to generate graphical output. See the set_cell_selection()
   * function for more information.
   */
  std::function<cell_iterator(const Triangulation<dim, spacedim> &)>
    first_cell_function;

  /**
   * A function object that is used to select what the next cell is going to
   * be on which to generate graphical output, given a previous cell. See
   * the set_cell_selection() function for more information.
   */
  std::function<cell_iterator(const Triangulation<dim, spacedim> &,
                              const cell_iterator &)>
    next_cell_function;

  /**
   * Build one patch. This function is called in a WorkStream context.
   *
   * The first argument here is the iterator, the second the scratch data
   * object. All following are tied to particular values when calling
   * WorkStream::run(). The function does not take a CopyData object but
   * rather allocates one on its own stack for memory access efficiency
   * reasons.
   */
  void
  build_one_patch(
    const std::pair<cell_iterator, unsigned int> *cell_and_index,
    internal::DataOutImplementation::ParallelData<dim, spacedim> &scratch_data,
    const unsigned int     n_subdivisions,
    const CurvedCellRegion curved_cell_region);
};


DEAL_II_NAMESPACE_CLOSE

#endif
